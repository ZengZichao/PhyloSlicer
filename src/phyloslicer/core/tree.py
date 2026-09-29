"""TreeArray: an array-based ultrametric-tree representation optimised for time slicing.

Node indexing convention (ape-like)
-----------------------------------
* tips are numbered ``0 .. n_tips - 1`` (index equals row in ``tip_labels``);
* internal nodes are numbered ``n_tips .. n_tips + n_internal - 1``;
* the root always has index ``n_tips``.

All depth quantities ("node depths") are measured from the root in branch-length
units: root depth is 0 and (for ultrametric trees) every tip has depth
``root_age``.  Ages in units of time before present are ``root_age - depth``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

import numpy as np

__all__ = ["TreeArray", "SlicedTree"]

#: Default *relative* tolerance of :attr:`TreeArray.is_ultrametric`.
#:
#: The test is ``spread(tip depths) <= ULTRAMETRIC_RTOL * max(tip depths)``,
#: i.e. both sides are in branch-length units and scale with the tree, so the
#: verdict is invariant to the unit the user happens to store branch lengths in
#: (Myr vs yr).  A fixed absolute tolerance would not be (see the review, A5).
ULTRAMETRIC_RTOL = 1e-8


@dataclass(eq=False)
class TreeArray:
    """Array representation of a (typically ultrametric) phylogenetic tree.

    Attributes
    ----------
    tip_labels : list of str
        Tip labels; the column index of a species in a presence matrix equals
        its position in this list.
    edges : (E, 2) int64 ndarray
        Each row is ``[parent, child]`` in the node numbering described above.
    lengths : (E,) float64 ndarray
        Branch lengths.
    root_edge : float
        Length of the stem above the root node, if the input Newick carried
        one.  It is **not** part of :attr:`edges`/:attr:`lengths` and therefore
        excluded from :meth:`total_pd`, matching ``ape::keep.tip(root.edge = 0)``;
        it is kept only so :meth:`to_newick` can restore it and so the excluded
        amount stays auditable.
    """

    tip_labels: list[str]
    edges: np.ndarray
    lengths: np.ndarray
    root_edge: float = 0.0
    _depths: np.ndarray | None = field(default=None, repr=False)
    _tip_intervals: np.ndarray | None = field(default=None, repr=False)
    _lengths_hash: int | None = field(default=None, repr=False)
    _edges_hash: tuple | None = field(default=None, repr=False)

    # ------------------------------------------------------------------ #
    # basic properties
    # ------------------------------------------------------------------ #
    @property
    def n_tips(self) -> int:
        return len(self.tip_labels)

    @property
    def n_edges(self) -> int:
        return len(self.lengths)

    @property
    def root(self) -> int:
        return self.n_tips

    def __post_init__(self) -> None:
        self.edges = np.asarray(self.edges, dtype=np.int64)
        self.lengths = np.asarray(self.lengths, dtype=np.float64)
        if self.edges.ndim != 2 or self.edges.shape[1] != 2:
            raise ValueError("edges must be a 2-column [parent, child] array")
        if len(self.lengths) != self.edges.shape[0]:
            raise ValueError("lengths must have one value per edge")
        if not np.isfinite(self.lengths).all():
            # NaN slips past every `<`/`>` guard, and an infinite branch makes
            # `spread <= rtol * scale` true because both sides are inf - so
            # without this a garbage tree is silently accepted *and* certified
            # as ultrametric, and only the window boundaries come back NaN.
            bad = int(np.count_nonzero(~np.isfinite(self.lengths)))
            raise ValueError(
                f"branch lengths must be finite numbers: {bad} NaN or Inf value(s) "
                "found; drop or repair those edges before building the tree"
            )
        self.root_edge = float(self.root_edge)
        if not np.isfinite(self.root_edge) or self.root_edge < 0.0:
            raise ValueError(
                f"root_edge must be a finite non-negative number, got {self.root_edge}"
            )
        if self.n_tips < 2:
            raise ValueError("trees need at least 2 tips")
        if (self.lengths < 0).any():
            raise ValueError("negative branch lengths are not supported")
        if len(set(self.tip_labels)) != self.n_tips:
            raise ValueError("tip labels must be unique")
        children = self.edges[:, 1]
        # the tip-interval representation requires tips to be numbered
        # 0..n_tips-1 with exactly one incoming edge each, and the root to be
        # id n_tips with no parent edge (internal ids above n_tips are free)
        if (children < 0).any() or (children == self.n_tips).any():
            raise ValueError(
                "invalid node numbering: child ids must be non-negative, and the "
                f"root id {self.n_tips} must not appear as a child"
            )
        tip_in = np.bincount(children[children < self.n_tips], minlength=self.n_tips)
        if (tip_in != 1).any():
            raise ValueError(
                "invalid node numbering: every tip 0..n_tips-1 needs exactly one "
                "incoming edge and internal node ids must be >= n_tips"
            )
        if not (self.edges[:, 0] == self.n_tips).any():
            raise ValueError(f"node {self.n_tips} (the root) has no outgoing edges")
        self._require_single_parent_and_acyclic()
        self._canonicalize_tip_order()
        visited = self._reachable_from_root()
        unknown = {int(v) for v in np.unique(self.edges)} - visited
        if unknown:
            raise ValueError(f"tree contains nodes unreachable from the root: {sorted(unknown)}")

    def _require_single_parent_and_acyclic(self) -> None:
        """Explicitly enforce the two structural invariants the kernels rely on.

        Every node may have at most one incoming edge (single-parent, i.e. a
        tree rather than a network), and the resulting graph must be acyclic.
        Both are assumed by :meth:`parent_of`, :meth:`node_depths` and
        :meth:`tip_intervals`, which each resolve a node to *one* parent.
        Previously the guarantee was only indirect (a duplicated tip edge trips
        the tip-in-degree check, a duplicated internal edge usually leaves the
        second parent unreachable), so the invariant could in principle be
        satisfied by accident rather than by construction.
        """
        children = self.edges[:, 1]
        n_nodes = int(self.edges.max()) + 1
        in_degree = np.bincount(children, minlength=n_nodes)
        duplicated = np.flatnonzero(in_degree > 1)
        if duplicated.size:
            raise ValueError(
                f"each node must have at most one parent edge; nodes with several "
                f"incoming edges: {sorted(int(v) for v in duplicated)} (reticulate "
                "input is not supported)"
            )
        # with unique parent edges a cycle can only be detected by walking
        # parent pointers; an iterative three-colour walk over the parent
        # function is O(n_nodes)
        parent = np.full(n_nodes, -1, dtype=np.int64)
        parent[children] = self.edges[:, 0]
        state = np.zeros(n_nodes, dtype=np.int8)  # 0 unseen, 1 on path, 2 done
        for start in range(n_nodes):
            if state[start]:
                continue
            path: list[int] = []
            node = start
            while node >= 0 and state[node] == 0:
                state[node] = 1
                path.append(node)
                node = int(parent[node])
            if node >= 0 and state[node] == 1:
                raise ValueError(
                    f"the edge list contains a cycle through node {node}; "
                    "a phylogenetic tree must be acyclic"
                )
            for done in path:
                state[done] = 2

    def _canonicalize_tip_order(self) -> None:
        """Renumber tips in DFS (preorder) order from the root.

        The tip-interval representation used across the package is only valid
        when every clade's tips form a contiguous index range, which requires
        DFS tip numbering.  Inputs built by hand or by simulators may use any
        order, so it is normalised here once at construction.
        """
        children = self.children_of()
        n_tips = self.n_tips
        dfs_tips: list[int] = []
        stack = [self.root]
        while stack:
            node = stack.pop()
            kids = children.get(node, [])
            if kids:
                stack.extend(reversed(kids))
            else:
                dfs_tips.append(node)
        if dfs_tips == list(range(n_tips)):
            return  # already canonical
        old_to_new = np.empty(int(self.edges.max()) + 1, dtype=np.int64)
        old_to_new[:n_tips] = -1
        for new_idx, old_idx in enumerate(dfs_tips):
            old_to_new[old_idx] = new_idx
        # internal ids are preserved: [n_tips, n_tips+1, ...] == arange
        old_to_new[n_tips:] = np.arange(n_tips, len(old_to_new))
        self.tip_labels = [self.tip_labels[old] for old in dfs_tips]
        remapped = self.edges.copy()
        child = remapped[:, 1]
        is_tip = child < n_tips
        remapped[is_tip, 1] = old_to_new[child[is_tip]]
        self.edges = remapped
        self._depths = None
        self._tip_intervals = None

    # ------------------------------------------------------------------ #
    # derived structure
    # ------------------------------------------------------------------ #
    def children_of(self) -> dict[int, list[int]]:
        """Map parent node -> ordered list of child nodes."""
        out: dict[int, list[int]] = {}
        for p, c in self.edges.tolist():
            out.setdefault(p, []).append(c)
        return out

    def _reachable_from_root(self) -> set[int]:
        children = self.children_of()
        seen: set[int] = {self.root}
        stack = [self.root]
        while stack:
            for child in children.get(stack.pop(), []):
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        return seen

    def parent_of(self) -> dict[int, int]:
        """Map child node -> parent node.

        ``__post_init__`` rejects any node with more than one incoming edge, so
        the dict comprehension below cannot silently overwrite a duplicate
        parent (it used to be the last-write-wins accident described in the
        review); it is exactly the inverse of the single-parent edge list.
        """
        return {c: p for p, c in self.edges.tolist()}

    def _check_cache(self) -> None:
        """Invalidate cached depths/intervals if edges or lengths were mutated."""
        edges_sig = (self.edges.shape, self.edges.tobytes())
        lengths_sig = self.lengths.tobytes() if self.lengths is not None else None
        if self._edges_hash is not None and self._edges_hash != edges_sig:
            self._depths = None
            self._tip_intervals = None
        if self._lengths_hash is not None and self._lengths_hash != lengths_sig:
            self._depths = None
        self._edges_hash = edges_sig
        self._lengths_hash = lengths_sig

    def node_depths(self) -> np.ndarray:
        """(n_nodes,) depth-from-root of every node (branch-length units)."""
        self._check_cache()
        if self._depths is None:
            n_nodes = int(self.edges.max()) + 1
            depths = np.zeros(n_nodes)
            parent = self.parent_of()
            edge_of_child = {int(c): i for i, c in enumerate(self.edges[:, 1].tolist())}
            children = self.children_of()
            stack = [self.root]
            while stack:
                node = stack.pop()
                for child in children.get(node, []):
                    depths[child] = depths[node] + self.lengths[edge_of_child[child]]
                    stack.append(child)
            missing = [n for n in range(n_nodes) if n != self.root and n not in parent]
            if missing:
                raise ValueError("tree contains disconnected nodes")
            self._depths = depths
        return self._depths

    def tip_intervals(self) -> np.ndarray:
        """(E, 2) half-open tip-index intervals ``[lo, hi)`` covered by each edge.

        Edge ``e`` lies on the root-to-tip path of tip ``t`` iff
        ``lo_e <= t < hi_e``; equivalently, the tips descending from the child
        node of ``e`` form the contiguous range ``tip_labels[lo:hi]``.  This
        single structure powers the fast site-edge membership construction
        used by every index.
        """
        self._check_cache()
        if self._tip_intervals is None:
            children = self.children_of()
            n_nodes = int(self.edges.max()) + 1
            lo = np.zeros(n_nodes, dtype=np.int64)
            hi = np.zeros(n_nodes, dtype=np.int64)
            # iterative post-order traversal
            order: list[int] = []
            stack = [self.root]
            visited: set[int] = set()
            while stack:
                node = stack[-1]
                kids = children.get(node, [])
                if node not in visited and kids:
                    stack.extend(reversed(kids))
                    visited.add(node)
                    continue
                stack.pop()
                if kids:
                    lo[node] = min(lo[c] for c in kids)
                    hi[node] = max(hi[c] for c in kids)
                else:  # tip
                    lo[node] = node
                    hi[node] = node + 1
                order.append(node)
            # note: post-order ensures children are finished before parents,
            # but the min/max above already read correct child values because
            # the loop pops children before parents.
            intervals = np.zeros((self.n_edges, 2), dtype=np.int64)
            edge_of_child = {int(c): i for i, c in enumerate(self.edges[:, 1].tolist())}
            for child, edge_i in edge_of_child.items():
                intervals[edge_i] = (lo[child], hi[child])
            self._tip_intervals = intervals
        return self._tip_intervals

    @property
    def root_age(self) -> float:
        """Age of the tree = the deepest tip depth, in branch-length units.

        For ultrametric trees every tip shares this depth; for
        non-ultrametric input it is the maximum tip depth (and a warning is
        expected from :meth:`check_ultrametric` when slicing).
        """
        return float(self.node_depths()[: self.n_tips].max())

    def tip_ages(self) -> np.ndarray:
        """Age (before present) of each tip = root_age - tip depth."""
        return self.root_age - self.node_depths()[: self.n_tips]

    @property
    def is_ultrametric(self) -> bool:
        """Whether all tips sit at the same depth from the root.

        The judgement is *relative to the size of the tree*: the tip-depth
        spread must not exceed ``ULTRAMETRIC_RTOL * root_age``.  Because both
        sides carry branch-length units, rescaling the tree (Myr -> yr) can
        never flip the verdict, which is the property an absolute tolerance
        lacks.  An earlier revision compared tip **ages** (= ``root_age`` -
        tip depth, hence identically 0 for an ultrametric tree) against their
        median, so ``rtol`` multiplied a quantity that is structurally zero and
        the check silently degenerated into a fixed ``1e-9`` absolute test.
        """
        return self.is_ultrametric_within(ULTRAMETRIC_RTOL)

    def is_ultrametric_within(self, rtol: float) -> bool:
        """:attr:`is_ultrametric` with an explicit relative tolerance."""
        depths = self.node_depths()[: self.n_tips]
        scale = float(depths.max())
        spread = float(depths.max() - depths.min())
        if not np.isfinite(scale) or not np.isfinite(spread):
            # inf/inf would otherwise satisfy `spread <= rtol * scale` and
            # certify a tree with an infinite branch as ultrametric
            return False
        if scale == 0.0:
            return spread == 0.0
        return bool(spread <= rtol * scale)

    def total_pd(self) -> float:
        """Total phylogenetic diversity (sum of all branch lengths)."""
        return float(self.lengths.sum())

    def check_ultrametric(self, allow_non_ultrametric: bool = False) -> bool:
        """Guard the ultrametricity assumption every slice rests on.

        Mirrors ``treesliceR::nodes_config``, which calls ``stop()`` on a
        non-ultrametric tree: the default here is a hard :class:`ValueError`.
        ``allow_non_ultrametric=True`` downgrades it to a warning and returns
        ``False`` so the flag can be recorded as an audit column by the
        callers; it returns ``True`` when the tree is ultrametric.

        Slicing coordinates (``root_age``, ``depth_for_pd``, ``windows_time``,
        :attr:`SliceStack.ages`) are all derived from "every tip ends at the
        same depth", so on a non-ultrametric tree the numbers stay computable
        but lose their meaning as "diversity before time *t*".
        """
        if self.is_ultrametric:
            return True
        if allow_non_ultrametric:
            import warnings

            warnings.warn(
                "the input tree is not ultrametric; slicing semantics assume "
                "tip ages are equal - results may be meaningless "
                "(pass allow_non_ultrametric=False to make this an error)",
                stacklevel=2,
            )
            return False
        raise ValueError(
            "the inputted phylogenetic tree is not ultrametric; time slicing "
            "requires all tips at the same depth. Rescale or re-ultrametrise "
            "the tree, or pass allow_non_ultrametric=True to proceed with a "
            "warning and an `ultrametric=False` audit flag."
        )

    def check_non_degenerate(self) -> None:
        """Reject trees whose branch lengths carry no temporal signal.

        A tree parsed without any ``:length`` field has ``root_age == 0`` and
        ``total_pd == 0``; every window over it is empty, so slicing would
        "succeed" while returning all-zero slices and all-zero PD profiles.
        Failing loudly is the documented alternative to that silent result.
        """
        if self.n_edges == 0:
            raise ValueError("the tree has no edges")
        if self.root_age <= 0.0:
            raise ValueError(
                "the tree has zero depth (every branch length is 0), so no "
                "time window can be placed in it; supply branch lengths or "
                "rescale before slicing"
            )
        if self.total_pd() <= 0.0:
            raise ValueError(
                "the tree stores no phylogenetic diversity (total PD is 0); "
                "PD-based criteria are undefined on it"
            )

    # ------------------------------------------------------------------ #
    # transforms
    # ------------------------------------------------------------------ #
    def subarray(self, tip_labels: Sequence[str]) -> "TreeArray":
        """Prune the tree down to a subset of tips (rooted at their MRCA).

        The root→MRCA stem is discarded, matching ``ape::keep.tip`` with
        ``root.edge = 0``.

        The returned tree's tips are in this tree's canonical (DFS) order,
        **not** in the order they were requested: ``subarray(["D", "C"])``
        yields ``["C", "D"]``.  Reorder any accompanying matrix by the result's
        :attr:`tip_labels` if you need a specific column order.
        """
        # one dict build instead of ``list.index`` per requested label: the
        # latter made this call itself O(n_kept * n_tips)
        positions = {label: i for i, label in enumerate(self.tip_labels)}
        try:
            keep = [positions[t] for t in tip_labels]
        except KeyError as exc:
            raise KeyError(f"tip label not present in the tree: {exc.args[0]!r}") from None
        if len(set(keep)) != len(keep):
            raise ValueError("duplicate tip labels in subset")
        if len(keep) < 2:
            raise ValueError(
                f"a subtree needs at least 2 tips but only {len(keep)} "
                f"({list(tip_labels)}) were requested"
            )
        return _prune_to_tips(self, np.asarray(sorted(keep), dtype=np.int64))

    def with_lengths(self, new_lengths: np.ndarray) -> "SlicedTree":
        """Return a :class:`SlicedTree` sharing this topology with new lengths."""
        return SlicedTree(parent=self, lengths=np.asarray(new_lengths, dtype=np.float64))

    def to_newick(self) -> str:
        from ..io.newick import write_newick

        return write_newick(self)

    # ------------------------------------------------------------------ #
    # constructors / adapters
    # ------------------------------------------------------------------ #
    @classmethod
    def from_newick(cls, source) -> "TreeArray":
        from ..io.newick import read_newick

        return read_newick(source)

    @classmethod
    def from_nexus(cls, path) -> list["TreeArray"]:
        from ..io.newick import read_nexus

        return read_nexus(path)

    @classmethod
    def from_dendropy(cls, obj) -> "TreeArray":
        from ..io.adapters import from_dendropy

        return from_dendropy(obj)

    @classmethod
    def from_biopython(cls, obj) -> "TreeArray":
        from ..io.adapters import from_biopython

        return from_biopython(obj)

    def to_dendropy(self):
        from ..io.adapters import to_dendropy

        return to_dendropy(self)

    def to_biopython(self):
        from ..io.adapters import to_biopython

        return to_biopython(self)

    # convenience
    def __len__(self) -> int:
        return self.n_tips

    def __eq__(self, other: object) -> bool:
        """Structural equality: same tip order, same topology, same branch lengths.

        ``TreeArray`` is declared ``eq=False`` because the dataclass-generated
        ``__eq__`` compares the ``edges``/``lengths`` ndarrays with ``==`` and
        then casts the element-wise result to ``bool``, which raises
        ``ValueError: The truth value of an array ... is ambiguous`` for every
        real tree - so ``a == b``, ``a in list_of_trees`` and ``assert got ==
        expected`` all failed with a message pointing at numpy, not here.

        The comparison goes through the Newick rendering rather than the raw
        rows, so two trees that store the same clades in a different edge order
        (which the parser and a hand-built array routinely do) still compare
        equal, while any difference in labels, topology, lengths or root stem
        does not.
        """
        if not isinstance(other, TreeArray):
            return NotImplemented
        if self.tip_labels != other.tip_labels or self.n_edges != other.n_edges:
            return False
        try:
            return self.to_newick() == other.to_newick()
        except ValueError:  # a structurally odd tree falls back to its arrays
            return bool(
                np.array_equal(self.edges, other.edges)
                and np.array_equal(self.lengths, other.lengths)
                and self.root_edge == other.root_edge
            )

    def __hash__(self) -> int:
        """Identity hash (a TreeArray is mutable: caches and in-place edits)."""
        return object.__hash__(self)

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"TreeArray(n_tips={self.n_tips}, n_edges={self.n_edges}, "
            f"root_age={self.root_age:.4g}, ultrametric={self.is_ultrametric})"
        )


@dataclass
class SlicedTree:
    """A tree slice: same topology as the parent tree, truncated branch lengths.

    Zero-length edges are retained so the slice can always be traced back to
    the parent tree; use :meth:`collapse` (or ``dropNodes=True`` in the
    compatibility layer) to remove them.
    """

    parent: TreeArray
    lengths: np.ndarray
    _collapsed: TreeArray | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if len(self.lengths) != self.parent.n_edges:
            raise ValueError("length vector must match the number of parent edges")

    @property
    def tree(self) -> TreeArray:
        """The (uncollapsed) sliced tree as a TreeArray."""
        return TreeArray(
            tip_labels=list(self.parent.tip_labels),
            edges=self.parent.edges.copy(),
            lengths=self.lengths.copy(),
        )

    @property
    def pd(self) -> float:
        """Phylogenetic diversity stored in this slice."""
        return float(self.lengths.sum())

    def collapse(self) -> TreeArray:
        """Collapse zero-length internal edges and re-number nodes."""
        if self._collapsed is None:
            self._collapsed = _collapse_zero_internal(
                self.parent.tip_labels, self.parent.edges, self.lengths
            )
        return self._collapsed

    def branch_length_hash(self) -> str:
        """Order-insensitive hash of the positive branch lengths (validation)."""
        import hashlib

        positive = np.sort(self.lengths[self.lengths > 0])
        payload = ",".join(f"{v:.10f}" for v in positive)
        return hashlib.sha256(payload.encode()).hexdigest()


# ---------------------------------------------------------------------- #
# structural helpers
# ---------------------------------------------------------------------- #
def _prune_to_tips(tree: TreeArray, tip_idx: np.ndarray) -> TreeArray:
    """Return the subtree spanning ``tip_idx``, rooted at their MRCA.

    The root→MRCA stem is **discarded** (matching ``ape::keep.tip`` with
    ``root.edge = 0``).  Unary internal nodes can remain when a sibling
    clade was pruned away; they are harmless for every kernel (PD,
    intervals, depths) and are removed later by :meth:`SlicedTree.collapse`
    if a compact tree is needed.

    Complexity is ``O(E log n_kept + E)``: the spanning-edge set comes out of
    two vectorised binary searches over the sorted kept-tip indices (no
    ``(E, n_kept)`` dense boolean mask, which was the memory peak), the MRCA
    out of one masked argmax over node depths (no per-node ``list.index``
    scan), and the stem removal out of an array-driven parent walk (no linear
    re-scan of the edge list per step).
    """
    if len(tip_idx) < 2:
        raise ValueError("a subtree needs at least 2 tips")
    edges = tree.edges
    children = edges[:, 1]
    n_tips = tree.n_tips
    n_nodes = int(edges.max()) + 1
    intervals = tree.tip_intervals()
    lo, hi = intervals[:, 0], intervals[:, 1]

    # ---- edges on some root -> kept-tip path --------------------------- #
    # child(node) of edge e owns the contiguous tip range [lo_e, hi_e), so e is
    # in the spanning subtree iff at least one kept tip index falls in it.
    tip_sorted = np.sort(np.asarray(tip_idx, dtype=np.int64))
    left = np.searchsorted(tip_sorted, lo, side="left")
    right = np.searchsorted(tip_sorted, hi, side="left")
    keep_edge = left < right

    # ---- identify the MRCA and drop the root→MRCA stem ------------------ #
    # With tips numbered in DFS order the MRCA is the deepest internal node
    # whose tip interval covers [min(kept), max(kept)]; every edge on the walk
    # from the MRCA up to the root is a stem edge and must be dropped.
    depths = tree.node_depths()
    internal_child = children >= n_tips
    covers_all = keep_edge & internal_child & (lo <= tip_sorted[0]) & (tip_sorted[-1] < hi)
    mrca = tree.root
    if covers_all.any():
        candidates = np.flatnonzero(covers_all)
        child_depths = depths[children[candidates]]
        mrca = int(children[candidates[int(np.argmax(child_depths))]])
    if mrca != tree.root:
        parent = np.full(n_nodes, -1, dtype=np.int64)
        edge_of_child = np.full(n_nodes, -1, dtype=np.int64)
        parent[children] = edges[:, 0]
        edge_of_child[children] = np.arange(children.size, dtype=np.int64)
        node = mrca
        while node != tree.root and node >= 0:
            e = int(edge_of_child[node])
            if e < 0:  # pragma: no cover - unreachable for valid trees
                break
            keep_edge[e] = False
            node = int(parent[node])

    kept = np.flatnonzero(keep_edge)
    sub_edges = edges[kept]

    # ---- renumber ------------------------------------------------------ #
    # the subtree root is the internal node that never appears as a child
    # among the kept edges; it receives id new_n_tips (the root convention).
    sub_children = sub_edges[:, 1]
    sub_internal_children = np.unique(sub_children[sub_children >= n_tips])
    old_internal = np.union1d(np.unique(sub_edges[:, 0]), sub_internal_children)
    subtree_roots = np.setdiff1d(old_internal, sub_internal_children)
    if subtree_roots.size != 1:
        raise ValueError("kept tips do not span a single subtree")
    new_n_tips = tip_sorted.size
    ordered = np.concatenate(([subtree_roots[0]], np.setdiff1d(old_internal, subtree_roots)))

    mapping = np.full(n_nodes, -1, dtype=np.int64)
    mapping[tip_sorted] = np.arange(new_n_tips, dtype=np.int64)
    mapping[ordered] = new_n_tips + np.arange(ordered.size, dtype=np.int64)
    new_edges = mapping[sub_edges]

    tip_labels = [tree.tip_labels[t] for t in tip_sorted.tolist()]
    return TreeArray(tip_labels=tip_labels, edges=new_edges, lengths=tree.lengths[kept].copy())


def _collapse_zero_internal(
    tip_labels: list[str], edges: np.ndarray, lengths: np.ndarray
) -> TreeArray:
    """Remove zero-length internal-node hops, re-attaching children upstream."""
    n_tips = len(tip_labels)
    child_to_parent = {c: p for p, c in edges.tolist()}
    removed = {
        c for (p, c), seg in zip(edges.tolist(), lengths.tolist()) if seg == 0.0 and c >= n_tips
    }

    def nearest_alive(node: int) -> int:
        cur = node
        while cur in removed:
            cur = child_to_parent[cur]
        return cur

    new_edges: list[tuple[int, int]] = []
    new_lengths: list[float] = []
    for (p, c), seg in zip(edges.tolist(), lengths.tolist()):
        if seg == 0.0 and c >= n_tips:
            continue
        new_edges.append((nearest_alive(p), c))
        new_lengths.append(seg)

    # a removed root hop would re-root the tree on the nearest alive ancestor;
    # this cannot happen for slices of an ultrametric tree (root keeps a
    # positive share of the deepest window) but guard anyway.
    old_internal = sorted(set(p for p, _ in new_edges) | {c for _, c in new_edges if c >= n_tips})
    internal_map = {old: n_tips + i for i, old in enumerate(old_internal)}
    out_edges = np.empty((len(new_edges), 2), dtype=np.int64)
    for i, (p, c) in enumerate(new_edges):
        out_edges[i, 0] = internal_map.get(p, p)
        out_edges[i, 1] = internal_map.get(c, c)
    return TreeArray(tip_labels=list(tip_labels), edges=out_edges, lengths=np.asarray(new_lengths))


def iter_edges(tree: TreeArray) -> Iterable[tuple[int, int, float]]:
    for (p, c), seg in zip(tree.edges.tolist(), tree.lengths.tolist()):
        yield p, c, float(seg)
