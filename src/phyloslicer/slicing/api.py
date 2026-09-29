"""Slicing API: rootward / tipward / interval / pieces / prune + SliceStack."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.tree import SlicedTree, TreeArray
from .windows import contribution_matrix, make_windows

__all__ = [
    "slice_rootward",
    "slice_tipward",
    "slice_interval",
    "slice_pieces",
    "prune_tips",
    "SliceStack",
]


# ---------------------------------------------------------------------- #
# single slices
# ---------------------------------------------------------------------- #
def _depth_from_time(tree: TreeArray, time: float, criterion: str) -> float:
    """Convert a user threshold into a depth-from-root cut position.

    * ``criterion="time"``: ``time`` is the window width in Ma before present;
      the cut depth is ``T - time`` (a rootward cut keeps ``[T - time, T]``, a
      tipward cut keeps ``[0, T - time]``).
    * ``criterion="pd"``: for rootward cuts ``time`` is the *root-side* PD
      amount to remove; for tipward cuts it is the *root-side* PD amount to
      keep (identical to treesliceR's behaviour, which converts the PD
      threshold into a depth measured from the root in both cases).

    In both cases the return value is a **depth** in ``[0, T]``; it is never a
    PD amount.  ``depth_for_pd`` already handles the two endpoints exactly
    (``0 -> 0``, ``total_pd -> T``), so no endpoint special-casing is needed
    (and an earlier revision that added one inverted all four endpoint
    behaviours - is what this branch prevents.)
    """
    T = tree.root_age
    if not np.isfinite(time):
        # both comparisons below are False for NaN, so without this a NaN
        # threshold passed straight through and yielded an all-NaN slice
        raise ValueError(f"the threshold must be a finite number, got {time!r}")
    if criterion == "time":
        if time < 0:
            raise ValueError("time must be non-negative")
        if time > T:
            raise ValueError(
                f"the threshold inputted ({time:g}) is bigger than the available "
                f"by the phylogeny ({T:g})"
            )
        # window in depth coordinates is [T - time, T] (rootward) or
        # [0, T - time] (tipward); time = 0 therefore yields an empty rootward
        # slice and the full tipward tree, exactly as the general formula does
        return T - time
    if criterion == "pd":
        total = tree.total_pd()
        if total <= 0:
            raise ValueError(
                "the tree stores no phylogenetic diversity (total PD is 0), so "
                "a PD threshold cannot be converted into a cut depth"
            )
        if time < 0 or time > total:
            raise ValueError(
                f"the threshold inputted ({time:g}) lies outside the available "
                f"PD range [0, {total:g}]"
            )
        from .windows import depth_for_pd

        return depth_for_pd(tree, time)
    raise ValueError("criterion must be 'time' or 'pd'")


def _keep_window(tree: TreeArray, lo: float, hi: float, collapse: bool):
    lengths = contribution_matrix(tree, np.array([[lo, hi]]))[:, 0]
    sliced = tree.with_lengths(lengths)
    return sliced.collapse() if collapse else sliced


def slice_rootward(
    tree: TreeArray,
    time: float,
    criterion: str = "time",
    collapse: bool = False,
    allow_non_ultrametric: bool = False,
) -> SlicedTree | TreeArray:
    """Keep the youngest portion of the tree (the crown slice of width ``time``).

    Mirrors treesliceR's ``squeeze_root``: with ``criterion="time"`` the result
    stores the branches younger than ``time`` Ma (the sliced tree's root starts
    at ``time`` Ma); with ``criterion="pd"`` the root-side PD amount ``time``
    is removed instead.

    ``time=0`` yields an empty slice and ``time=T`` (or ``time=total_pd`` for
    ``criterion="pd"``) the whole tree, matching ``squeeze_root``.

    Raises :class:`ValueError` for a non-ultrametric tree, as ``nodes_config()``
    does in R; pass ``allow_non_ultrametric=True`` for the legacy warn-and-go
    behaviour.
    """
    tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    tree.check_non_degenerate()
    x = _depth_from_time(tree, time, criterion)
    return _keep_window(tree, x, tree.root_age, collapse)


def slice_tipward(
    tree: TreeArray,
    time: float,
    criterion: str = "time",
    collapse: bool = False,
    allow_non_ultrametric: bool = False,
) -> SlicedTree | TreeArray:
    """Keep the oldest portion of the tree (drop the youngest ``time`` Ma).

    Mirrors treesliceR's ``squeeze_tips``; ``criterion="pd"`` keeps the
    root-side PD amount ``time``.  ``time=0`` yields the whole tree and
    ``time=T`` (``total_pd`` for PD) an empty slice.
    """
    tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    tree.check_non_degenerate()
    x = _depth_from_time(tree, time, criterion)
    return _keep_window(tree, 0.0, x, collapse)


def slice_interval(
    tree: TreeArray,
    start: float,
    stop: float,
    invert: bool = False,
    criterion: str = "time",
    collapse: bool = False,
    allow_non_ultrametric: bool = False,
) -> SlicedTree | TreeArray | tuple:
    """Extract the temporal interval [stop, start] Ma (start > stop for "time").

    With ``invert=True`` returns the complementary ``(root_piece, tip_piece)``
    tuple: everything older than ``start`` and everything younger than ``stop``.
    Note that treesliceR's ``squeeze_int(invert = TRUE)`` returns the same two
    trees in the opposite order; :func:`phyloslicer.compat.squeeze_int` follows
    R, this function follows the documented ``(older, younger)`` convention.
    """
    tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    tree.check_non_degenerate()
    if not np.isfinite(start) or not np.isfinite(stop):
        raise ValueError(f"start and stop must be finite numbers, got {start!r} and {stop!r}")
    if criterion == "time":
        if start <= stop:
            raise ValueError("the thresholds set in arguments [start] and [stop] are incompatible")
        if stop < 0 or start > tree.root_age:
            raise ValueError(
                f"the interval [{stop:g}, {start:g}] Ma lies outside the phylogeny "
                f"(root age {tree.root_age:g})"
            )
    elif criterion == "pd":
        total = tree.total_pd()
        if total <= 0:
            raise ValueError("the tree stores no phylogenetic diversity (total PD is 0)")
        if start >= stop:
            raise ValueError("the thresholds set in arguments [start] and [stop] are incompatible")
        if start < 0 or stop > total:
            raise ValueError(
                f"the interval [{start:g}, {stop:g}] PD lies outside the phylogeny "
                f"(total PD {total:g})"
            )
    else:
        raise ValueError("criterion must be 'time' or 'pd'")

    if invert:
        older = slice_tipward(
            tree,
            start,
            criterion=criterion,
            collapse=collapse,
            allow_non_ultrametric=allow_non_ultrametric,
        )
        younger = slice_rootward(
            tree,
            stop,
            criterion=criterion,
            collapse=collapse,
            allow_non_ultrametric=allow_non_ultrametric,
        )
        return older, younger
    # interior interval [lo, hi] in depth coordinates
    if criterion == "time":
        lo = tree.root_age - start
        hi = tree.root_age - stop
    else:  # "pd": start < stop are root-side PD amounts
        from .windows import depth_for_pd

        lo = depth_for_pd(tree, start)
        hi = depth_for_pd(tree, stop)
    return _keep_window(tree, lo, hi, collapse)


# ---------------------------------------------------------------------- #
# SliceStack
# ---------------------------------------------------------------------- #
@dataclass
class SliceStack:
    """A stack of k slices of one tree, without copying the tree.

    Attributes
    ----------
    tree : TreeArray
        The full tree being sliced.
    depth_windows : (k, 2) ndarray
        Slice windows ``[lo, hi]`` in depth-from-root coordinates (root -> tip).
    criterion : str
        "time" or "pd".
    ultrametric : bool
        Audit flag: whether the sliced tree passed the ultrametricity check.
        It is ``False`` only when the check was explicitly waived with
        ``allow_non_ultrametric=True``.
    """

    tree: TreeArray
    depth_windows: np.ndarray
    criterion: str = "time"
    ultrametric: bool = True
    _contrib: np.ndarray | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.depth_windows = np.asarray(self.depth_windows, dtype=np.float64)

    # -- geometry ------------------------------------------------------ #
    @property
    def n_slices(self) -> int:
        return len(self.depth_windows)

    @property
    def windows(self) -> np.ndarray:
        """(k, 2) windows in Ma before present, ``[older, younger]``, root -> tip."""
        T = self.tree.root_age
        return np.column_stack([T - self.depth_windows[:, 0], T - self.depth_windows[:, 1]])

    @property
    def ages(self) -> np.ndarray:
        """(k,) older boundary of each slice in Ma before present (root -> tip)."""
        return self.windows[:, 0]

    @property
    def widths(self) -> np.ndarray:
        """(k,) width of each slice **in depth-from-root units** (root -> tip).

        These are always depth/height differences, so for ``criterion="time"``
        they are the (equal) durations in Ma, while for ``criterion="pd"`` the
        windows are equal in *accumulated PD* and therefore unequal here.  The
        PD held by each window is ``contribution_matrix().sum(axis=0)``.
        """
        return self.depth_windows[:, 1] - self.depth_windows[:, 0]

    def contribution_matrix(self) -> np.ndarray:
        """(E, k) edge x slice contribution matrix (cached)."""
        if self._contrib is None:
            self._contrib = contribution_matrix(self.tree, self.depth_windows)
        return self._contrib

    def contribution(self, tip_subset=None) -> np.ndarray:
        """(k,) branch-length contribution of a set of tips per slice.

        ``tip_subset`` is a list of tip labels or indices; ``None`` uses all
        tips.  This is the core abstraction replacing treesliceR's k full-tree
        copies: one matrix multiply instead of k tree objects.
        """
        C = self.contribution_matrix()
        if tip_subset is None:
            return C.sum(axis=0)
        positions = {label: i for i, label in enumerate(self.tree.tip_labels)}
        idx = np.asarray(
            [positions[t] if isinstance(t, str) else int(t) for t in tip_subset],
            dtype=np.int64,
        )
        intervals = self.tree.tip_intervals()
        lo, hi = intervals[:, 0], intervals[:, 1]
        edge_mask = ((lo[None, :] <= idx[:, None]) & (idx[:, None] < hi[None, :])).any(axis=0)
        return edge_mask.astype(float) @ C

    def to_sliced_trees(self, collapse: bool = False) -> list[SlicedTree | TreeArray]:
        """Materialise the k slices as tree objects (only when really needed)."""
        C = self.contribution_matrix()
        out = [self.tree.with_lengths(C[:, j]) for j in range(self.n_slices)]
        if collapse:
            return [s.collapse() for s in out]
        return out

    # -- convenience --------------------------------------------------- #
    def __len__(self) -> int:
        return self.n_slices

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        return (
            f"SliceStack(n_slices={self.n_slices}, criterion={self.criterion!r}, "
            f"n_tips={self.tree.n_tips}, ultrametric={self.ultrametric})"
        )


def slice_pieces(
    tree: TreeArray,
    n: int | None = None,
    width: float | None = None,
    criterion: str = "time",
    allow_non_ultrametric: bool = False,
) -> SliceStack:
    """Cut a tree into multiple slices (mirrors treesliceR's ``phylo_pieces``).

    Exactly one of ``n`` (number of slices) or ``width`` (slice width in Ma or
    PD units; rounded treesliceR-style to the nearest achievable value) must
    be given - replacing R's easily-confused ``n`` + ``method`` pair.

    A non-ultrametric tree is rejected (as ``nodes_config()`` does in R) unless
    ``allow_non_ultrametric=True``; the resulting :class:`SliceStack` carries
    the verdict in its ``ultrametric`` attribute.
    """
    ultrametric = tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    tree.check_non_degenerate()
    windows = make_windows(tree, n=n, width=width, criterion=criterion)
    return SliceStack(
        tree=tree, depth_windows=windows, criterion=criterion, ultrametric=bool(ultrametric)
    )


# ---------------------------------------------------------------------- #
# tip pruning
# ---------------------------------------------------------------------- #
def prune_tips(
    tree: TreeArray,
    threshold: float | np.ndarray,
    quantiles: bool = False,
    side: str = "after",
    allow_non_ultrametric: bool = False,
) -> TreeArray | list[TreeArray]:
    """Prune tips by their terminal-branch length (mirrors treesliceR's ``prune_tips``).

    Parameters
    ----------
    threshold : float or array
        Terminal-branch-length threshold(s); with ``quantiles=True`` the
        value(s) are interpreted as quantile probabilities in [0, 1] of the
        observed terminal branch lengths.
    side : {"after", "before"}
        ``"after"`` keeps tips whose terminal branch is >= the threshold
        (treesliceR ``method=1``); ``"before"`` drops them (``method=2``).
        Any other value raises; R's ``prune_tips`` returns ``NULL`` for a
        ``method`` outside {1, 2}.

    Raises
    ------
    ValueError
        when a threshold would leave fewer than two tips.  R returns ``NULL``
        there; this raises, because a ``None`` reaching a downstream ``.tip_labels``
        is a harder bug to trace than a named threshold.
    """
    if side not in ("after", "before"):
        raise ValueError("side must be 'after' or 'before'")
    # R's prune_tips routes through nodes_config(), which stops on a
    # non-ultrametric tree; the terminal-branch threshold is only interpretable
    # as a "splitting time" when every tip sits at the same depth.
    tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    tree.check_non_degenerate()

    child = tree.edges[:, 1]
    tip_edge_mask = child < tree.n_tips
    if tip_edge_mask.sum() != tree.n_tips:
        # every tip must have exactly one incoming (terminal) edge
        raise ValueError("unexpected tree shape: tips without a unique terminal edge")
    tip_branch = np.full(tree.n_tips, np.nan)
    tip_branch[child[tip_edge_mask]] = tree.lengths[tip_edge_mask]

    thresholds = np.atleast_1d(np.asarray(threshold, dtype=np.float64))
    if quantiles:
        if (thresholds < 0).any() or (thresholds > 1).any():
            raise ValueError("quantile values must fall within the range 0 to 1")
        thresholds = np.quantile(tip_branch, thresholds)

    def one(thr: float) -> TreeArray:
        """Prune for one threshold, or explain why the cut is impossible.

        An infeasible threshold used to come back as a bare ``None`` - the same
        silent ``NULL`` treesliceR is criticised for, and which this package's
        own compat layer documents as "worse for a scripted migration than an
        error".  ``ValueError`` now names the threshold and the range that
        would work.
        """
        if side == "after":
            keep = np.flatnonzero(tip_branch >= thr)
            why = (
                f"threshold {thr:g} leaves {keep.size} tip(s) with a terminal "
                f"branch of at least that length (the longest terminal branch in "
                f"the tree is {tip_branch.max():g})"
            )
        else:
            keep = np.flatnonzero(tip_branch < thr)
            why = (
                f"threshold {thr:g} prunes every tip below it, leaving {keep.size} "
                f"(the shortest terminal branch in the tree is {tip_branch.min():g})"
            )
        if keep.size < 2:
            raise ValueError(
                f"{why}; a pruned tree still needs at least 2 tips. "
                f"Use a threshold in ({tip_branch.min():g}, {tip_branch.max():g}) "
                "for this tree, or pass quantiles=True to express it as a quantile"
            )
        return tree.subarray([tree.tip_labels[i] for i in keep])

    results = [one(float(t)) for t in thresholds]
    scalar = thresholds.size == 1
    if scalar:
        return results[0]
    return results
