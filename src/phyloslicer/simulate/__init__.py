"""Tree and community simulation for tests, benchmarks and tutorials."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..core.tree import TreeArray
from ..io.newick import write_newick

__all__ = ["yule_tree", "birth_death_tree", "make_dataset", "to_newick_file"]


def yule_tree(
    n: int, age: float = 1.0, seed: int | None = None, tip_prefix: str = "sp_"
) -> TreeArray:
    """Ultrametric tree of exactly ``n`` tips via the backward-Yule (coalescent) process."""
    if n < 2:
        raise ValueError("n must be >= 2")
    rng = np.random.default_rng(seed)
    heights = [0.0] * n
    active = list(range(n))
    edges: list[tuple[int, int]] = []
    lengths: list[float] = []
    node_id = n
    t = 0.0
    while len(active) > 1:
        rate = len(active) * (len(active) - 1) / 2.0
        t += rng.exponential(1.0 / rate)
        i, j = rng.choice(len(active), size=2, replace=False)
        a, b = active[int(i)], active[int(j)]
        for child, h in ((a, heights[a]), (b, heights[b])):
            edges.append((node_id, child))
            lengths.append(t - h)
        active = [x for k, x in enumerate(active) if k not in (int(i), int(j))]
        active.append(node_id)
        heights.append(t)
        node_id += 1
    total = t
    scale = age / total if total > 0 else 1.0
    # internal nodes were created n, n+1, ..., with the final MRCA last; the
    # TreeArray convention puts the root at id n_tips, so renumber in reverse
    # creation order (MRCA -> n, first coalescence -> 2n-2).
    old_internal = list(range(n, 2 * n - 1))
    new_id = {old: n + (old_internal[-1] - old) for old in old_internal}
    edges = [(new_id.get(p, p), new_id.get(c, c)) for p, c in edges]
    return TreeArray(
        tip_labels=[f"{tip_prefix}{i + 1}" for i in range(n)],
        edges=np.asarray(edges, dtype=np.int64),
        lengths=np.asarray(lengths) * scale,
    )


class _Node:
    __slots__ = ("t0", "parent", "children")

    def __init__(self, t0: float, parent: "_Node | None"):
        self.t0 = t0
        self.parent = parent
        self.children: list[_Node] = []


def birth_death_tree(
    n: int,
    age: float = 1.0,
    birth_rate: float = 1.0,
    death_rate: float = 0.2,
    seed: int | None = None,
    max_attempts: int = 500,
    max_lineages: int = 200_000,
    tip_prefix: str = "sp_",
) -> TreeArray:
    """Ultrametric tree of exactly ``n`` tips under a forward birth-death process.

    Lineages are born at ``birth_rate`` and die at ``death_rate`` starting from
    one lineage; the simulation stops the moment the extant count reaches ``n``,
    dead side lineages are pruned, and the reconstructed tree of the survivors
    is rescaled so its depth is ``age``.  The shape therefore comes from the
    rates and the depth from ``age``, exactly as in :func:`yule_tree`.

    ``max_attempts`` bounds the retries after an attempt that goes extinct
    before ever reaching ``n``; ``max_lineages`` ceilings the live-lineage count
    inside one attempt, so no parameter choice can turn the call into an
    unbounded computation.

    Raises
    ------
    ValueError
        for ``n < 2``, non-positive ``age`` or ``birth_rate``, or a negative
        ``death_rate``.
    RuntimeError
        when no attempt reached ``n`` lineages within ``max_attempts``.  The
        usual cause is ``death_rate >= birth_rate``, which makes extinction
        likely before the count ever grows.
    """
    if n < 2:
        raise ValueError("n must be >= 2")
    if age <= 0.0:
        raise ValueError("age must be positive")
    if birth_rate <= 0.0:
        raise ValueError("birth_rate must be positive")
    if death_rate < 0.0:
        raise ValueError("death_rate must be non-negative")
    if death_rate >= birth_rate:
        # a critical or subcritical process goes extinct with probability 1
        # (mu > lambda) or creeps up only diffusively (mu == lambda), so
        # reaching n may take many attempts - worth saying before the caller
        # waits out a RuntimeError
        import warnings

        warnings.warn(
            f"death_rate={death_rate:g} is not below birth_rate={birth_rate:g}: the "
            f"lineage count tends to shrink, so reaching {n} extant tips may take "
            f"many attempts (raise max_attempts or lower death_rate)",
            stacklevel=2,
        )
    rng = np.random.default_rng(seed)
    for _ in range(max_attempts):
        built = _one_bd_simulation(n, age, birth_rate, death_rate, rng, tip_prefix, max_lineages)
        if built is not None:
            edges, lengths, labels = built
            return TreeArray(tip_labels=labels, edges=edges, lengths=lengths)
    raise RuntimeError(
        f"could not grow a birth-death tree to {n} extant tips within "
        f"{max_attempts} attempts at birth_rate={birth_rate:g}, "
        f"death_rate={death_rate:g}; with death_rate at or above birth_rate the "
        f"process usually dies out first - lower death_rate, raise "
        f"max_attempts, or call yule_tree(n, age) for pure birth"
    )


def _one_bd_simulation(n, age, lam, mu, rng, tip_prefix="sp_", max_lineages=200_000):
    root = _Node(0.0, None)
    alive = [root]
    t = 0.0
    # Grow until the extant count first reaches n.  Each birth adds exactly one
    # lineage, so the loop ends with len(alive) == n after O(n) events.  Running
    # to a fixed `age` and then *rejecting* every attempt whose count differed
    # from n - the previous design - never terminated usefully for a
    # supercritical rate: the expected count at that age is e**((lambda-mu)*age),
    # astronomically above n, so n >= 20 was unattainable while one attempt
    # itself grew without bound.
    while len(alive) < n:
        n_alive = len(alive)
        if n_alive > max_lineages:
            return None
        total_rate = (lam + mu) * n_alive
        if total_rate <= 0:
            return None
        t += rng.exponential(1.0 / total_rate)
        parent = alive[int(rng.integers(n_alive))]
        if rng.random() < lam / (lam + mu):  # birth
            child = _Node(t, parent)
            parent.children.append(child)
            alive.append(child)
        else:  # death - the lineage leaves the extant pool
            k = int(rng.integers(n_alive))
            alive.pop(k)
        if not alive:
            return None  # this attempt went extinct; the caller retries
    horizon = t
    survivors = alive
    if len(survivors) != n:
        return None

    survivor_ids = {id(node) for node in survivors}

    # retain only ancestors of survivors (the MRCA may be younger than the root)
    retained: set[int] = set()
    for node in survivors:
        cur = node
        while cur is not None:
            retained.add(id(cur))
            cur = cur.parent

    retained_nodes: list[_Node] = []

    def gather(node: _Node) -> None:
        if id(node) in retained:
            retained_nodes.append(node)
        for child in node.children:
            gather(child)

    gather(root)
    has_retained_child = {
        id(node): any(id(c) in retained for c in node.children) for node in retained_nodes
    }

    # A lineage alive at the horizon that also split earlier is BOTH an extant
    # species (its continuing branch) and the attachment point of its
    # daughters; represent the continuing branch as a synthetic terminal
    # hanging on the split node with length (horizon - split time).
    internals = [node for node in retained_nodes if has_retained_child[id(node)]]
    internals.sort(key=lambda nd: 0 if nd.parent is None else 1)  # root first
    synthetic_tips = [node for node in internals if id(node) in survivor_ids]
    leaf_nodes = [node for node in retained_nodes if not has_retained_child[id(node)]]

    n_tips = len(leaf_nodes) + len(synthetic_tips)
    if n_tips != n:
        return None

    internal_ids = {id(node): n_tips + i for i, node in enumerate(internals)}
    tip_ids: dict = {}
    tip_counter = 0
    for node in leaf_nodes:
        tip_ids[id(node)] = tip_counter
        tip_counter += 1
    for node in synthetic_tips:
        tip_ids[("synth", id(node))] = tip_counter
        tip_counter += 1

    def nearest_retained_ancestor(node: _Node) -> _Node | None:
        cur = node.parent
        while cur is not None and id(cur) not in retained:
            cur = cur.parent
        return cur  # None only for the root itself

    edges: list[tuple[int, int]] = []
    lengths: list[float] = []

    def emit(node: _Node) -> None:
        parent = nearest_retained_ancestor(node)
        if id(node) in internal_ids:  # internal node: its incoming edge uses the internal id
            if parent is not None:
                edges.append((internal_ids[id(parent)], internal_ids[id(node)]))
                lengths.append(node.t0 - parent.t0)
            for child in node.children:
                if id(child) in retained:
                    emit(child)
        else:  # terminal lineage: branch runs from the parent down to the horizon
            my_id = tip_ids[id(node)]
            if parent is not None:
                edges.append((internal_ids[id(parent)], my_id))
                lengths.append(horizon - parent.t0)

    emit(root)
    for node in synthetic_tips:
        edges.append((internal_ids[id(node)], tip_ids[("synth", id(node))]))
        lengths.append(horizon - node.t0)

    labels = [f"{tip_prefix}{i + 1}" for i in range(n_tips)]
    out = np.asarray(lengths, dtype=np.float64)
    if horizon > 0:
        out = out * (age / horizon)  # the rates give the shape, `age` gives the depth
    return (
        np.asarray(edges, dtype=np.int64).reshape(-1, 2),
        out,
        labels,
    )


def make_dataset(
    n_sites: int = 100,
    n_species: int = 50,
    structure: str = "random",
    seed: int | None = None,
    grid_side: int | None = None,
    tip_labels: list[str] | None = None,
):
    """Simulate a presence/absence matrix with optional spatial clustering.

    Returns ``(mat: DataFrame sites x species, coords: DataFrame with x, y)``.
    ``structure="clustered"`` builds species ranges as random disks so nearby
    sites share species (spatial autocorrelation); ``"random"`` draws
    Bernoulli presences.  When ``tip_labels`` is given the matrix columns
    are set to those labels (useful for direct alignment with a simulated
    tree's ``tip_labels``); otherwise columns are ``"1".."n_species"``.
    """
    if structure not in ("random", "clustered"):
        raise ValueError("structure must be 'random' or 'clustered'")
    rng = np.random.default_rng(seed)
    if grid_side is None:
        grid_side = int(np.ceil(np.sqrt(n_sites)))
    xs = np.arange(grid_side, dtype=float)
    gx, gy = np.meshgrid(xs, xs)
    coords = pd.DataFrame({"x": gx.ravel()[:n_sites], "y": gy.ravel()[:n_sites]})

    if structure == "random":
        mat = (rng.random((n_sites, n_species)) < 0.3).astype(float)
    else:
        mat = np.zeros((n_sites, n_species))
        for s in range(n_species):
            cx, cy = rng.uniform(0, grid_side - 1, size=2)
            radius = rng.uniform(1.0, max(2.0, grid_side / 3.0))
            inside = (coords["x"] - cx) ** 2 + (coords["y"] - cy) ** 2 <= radius**2
            mat[inside.to_numpy(), s] = 1.0
        empty = mat.sum(axis=0) == 0
        if empty.any():
            site = rng.integers(0, n_sites, size=int(empty.sum()))
            mat[site, np.flatnonzero(empty)] = 1.0
    if tip_labels is not None:
        if len(tip_labels) != n_species:
            raise ValueError("tip_labels length must equal n_species")
        columns = list(tip_labels)
    else:
        columns = [str(i + 1) for i in range(n_species)]
    mat_df = pd.DataFrame(mat, columns=columns)
    return mat_df, coords


def to_newick_file(tree: TreeArray, path) -> None:
    with open(path, "w") as fh:
        fh.write(write_newick(tree) + "\n")
