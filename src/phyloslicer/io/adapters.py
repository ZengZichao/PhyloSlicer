"""Adapters for DendroPy and Biopython tree objects (optional dependencies)."""

from __future__ import annotations

import numpy as np

from ..core.tree import TreeArray

__all__ = [
    "from_dendropy",
    "to_dendropy",
    "from_biopython",
    "to_biopython",
]


def from_dendropy(obj) -> TreeArray:
    """Convert a ``dendropy.Tree`` to a :class:`TreeArray`.

    Tip order is the tree's own leaf order.  The labels come from the taxa
    *attached to those leaves*, not from ``obj.taxon_namespace``: posterior
    pools, ``split_partitions()`` and NEXUS files with a wider ``TAXA`` block
    all make several trees share one namespace, and taking the namespace as the
    tip set used to inflate ``n_tips`` so that TreeArray's structural check
    rejected a perfectly good tree with "every tip needs exactly one incoming
    edge".
    """
    import warnings

    try:
        import dendropy  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "install dendropy to use this adapter: pip install phyloslicer[io]"
        ) from exc

    leaves = list(obj.leaf_node_iter())
    if len(leaves) < 2:
        raise ValueError(f"a sliceable tree needs at least 2 tips, this one has {len(leaves)}")
    labels: list[str] = []
    tip_index: dict[int, int] = {}
    for i, node in enumerate(leaves):
        taxon = node.taxon
        if taxon is None or not taxon.label:
            # a leaf with no taxon at all used to raise a bare KeyError: None
            labels.append(f"tip_{i + 1}")
        else:
            labels.append(taxon.label)
        tip_index[id(node)] = i
    used = {id(n.taxon) for n in leaves if n.taxon is not None}
    unused = [t.label for t in obj.taxon_namespace if id(t) not in used]
    if unused:
        warnings.warn(
            f"the DendroPy taxon namespace carries {len(unused)} taxon(a) this tree "
            f"does not use and they are not imported as tips: "
            f"{', '.join(map(str, unused[:5]))}{'...' if len(unused) > 5 else ''}",
            stacklevel=2,
        )
    edges: list[tuple[int, int]] = []
    lengths: list[float] = []
    internal_ids: dict[int, int] = {}
    n_tips = len(labels)

    def node_id(node) -> int:
        if node.is_leaf():
            return tip_index[id(node)]
        if id(node) not in internal_ids:
            internal_ids[id(node)] = n_tips + len(internal_ids)
        return internal_ids[id(node)]

    # root has no parent edge; internal root id must be n_tips per convention
    root_id = n_tips
    if id(obj.seed_node) in internal_ids:
        raise ValueError("unexpected seed node assignment")
    internal_ids[id(obj.seed_node)] = root_id

    for node in obj.preorder_node_iter():
        if node is obj.seed_node:
            continue
        child = node_id(node)
        parent = node_id(node.parent_node)
        edges.append((parent, child))
        lengths.append(float(node.edge_length) if node.edge_length is not None else 0.0)

    return TreeArray(
        tip_labels=labels,
        edges=np.asarray(edges, dtype=np.int64),
        lengths=np.asarray(lengths, dtype=np.float64),
    )


def to_dendropy(tree: TreeArray):
    """Convert a :class:`TreeArray` to a ``dendropy.Tree``."""
    import dendropy

    dnd = dendropy.Tree()
    dnd.taxon_namespace = dendropy.TaxonNamespace(tree.tip_labels)
    taxa = {name: dnd.taxon_namespace.get_taxon(name) for name in tree.tip_labels}

    nodes: dict[int, object] = {}
    children = tree.children_of()

    def build(node_id: int):
        from dendropy import Node

        if node_id not in nodes:
            nodes[node_id] = Node()
        node = nodes[node_id]
        if node_id < tree.n_tips:
            node.taxon = taxa[tree.tip_labels[node_id]]
            node.label = tree.tip_labels[node_id]
        for child in children.get(node_id, []):
            child_node = build(child)
            child_node.edge_length = float(tree.lengths[(tree.edges[:, 1] == child)][0])
            node.add_child(child_node)
        return node

    dnd.seed_node = build(tree.root)
    dnd.label = "phyloslicer_tree"
    return dnd


def from_biopython(obj) -> TreeArray:
    """Convert a ``Bio.Phylo`` tree to a :class:`TreeArray`."""
    tips = list(obj.get_terminals())
    labels = [t.name or f"tip_{i + 1}" for i, t in enumerate(tips)]
    tip_index = {id(t): i for i, t in enumerate(tips)}
    n_tips = len(tips)

    edges: list[tuple[int, int]] = []
    lengths: list[float] = []
    internal_ids: dict[int, int] = {id(obj.root): n_tips}
    counter = 0

    parent_of: dict[int, int] = {}
    for clade in obj.find_clades(order="preorder"):
        for child in clade.clades:
            parent_of[id(child)] = id(clade)

    for clade in obj.find_clades(order="preorder"):
        if clade is obj.root:
            continue
        if clade.is_terminal():
            cid = tip_index[id(clade)]
        else:
            counter += 1
            cid = n_tips + counter
            internal_ids[id(clade)] = cid
        pid = internal_ids[parent_of[id(clade)]]
        edges.append((pid, cid))
        lengths.append(float(clade.branch_length) if clade.branch_length is not None else 0.0)

    return TreeArray(
        tip_labels=labels,
        edges=np.asarray(edges, dtype=np.int64),
        lengths=np.asarray(lengths, dtype=np.float64),
    )


def to_biopython(tree: TreeArray):
    from Bio.Phylo.Newick import Clade, Tree

    children = tree.children_of()

    def build(node_id: int, length: float | None = None) -> Clade:
        if node_id < tree.n_tips:
            clade = Clade(branch_length=length, name=tree.tip_labels[node_id])
            return clade
        clade = Clade(branch_length=length)
        for child in children.get(node_id, []):
            child_len = float(tree.lengths[(tree.edges[:, 1] == child)][0])
            clade.clades.append(build(child, child_len))
        return clade

    root = build(tree.root, None)
    return Tree(root=root, rooted=True)
