"""Interval-intersection kernels: the mathematical core of PhyloSlicer.

Every edge of the tree is the time interval ``[b_e, d_e]`` in depth-from-root
coordinates (root = 0, tips = ``root_age`` T).  A slice is a window ``[lo, hi]``
in the same coordinates, and the branch length contributed by edge ``e`` to
window ``j`` is simply the length of the intersection

    C[e, j] = max(0, min(d_e, hi_j) - max(b_e, lo_j))

which can be computed for all edges and all windows at once with array
broadcasting - no tree is ever copied (contrast treesliceR's
``phylo_pieces``, which materialises k full tree copies).
"""

from __future__ import annotations

import numpy as np

from ..core.tree import TreeArray

__all__ = [
    "contribution_matrix",
    "pd_depth_curve",
    "depth_for_pd",
    "depths_for_pd",
    "windows_time",
    "windows_pd",
    "make_windows",
]


def contribution_matrix(tree: TreeArray, windows_depth: np.ndarray) -> np.ndarray:
    """(E, k) matrix of edge contributions to each window.

    Parameters
    ----------
    tree : TreeArray
    windows_depth : (k, 2) ndarray
        Windows ``[lo, hi]`` in depth-from-root coordinates, root -> tip order.
    """
    windows_depth = np.asarray(windows_depth, dtype=np.float64)
    if windows_depth.ndim != 2 or windows_depth.shape[1] != 2:
        raise ValueError("windows must be a (k, 2) [lo, hi] array")
    if (windows_depth[:, 1] < windows_depth[:, 0]).any():
        raise ValueError("window hi must be >= lo")
    depths = tree.node_depths()
    b = depths[tree.edges[:, 0]]
    d = depths[tree.edges[:, 1]]
    lo = windows_depth[None, :, 0]
    hi = windows_depth[None, :, 1]
    return np.clip(np.minimum(d[:, None], hi) - np.maximum(b[:, None], lo), 0.0, None)


def pd_depth_curve(tree: TreeArray) -> tuple[np.ndarray, np.ndarray]:
    """Cumulative-PD-vs-depth function of the tree.

    Returns
    -------
    (breaks, cum_pd) : (m,) breakpoints (incl. 0 and root_age) and the
        cumulative PD accumulated from the root to each breakpoint.

    The breakpoint set is every node depth (both ends of every edge), not just
    the parents' depths: on a non-ultrametric tree a tip can end above
    ``root_age``, and without that depth as a boundary the branch would be
    counted as alive over a segment it does not cover, so ``cum_pd`` stopped
    agreeing with :meth:`TreeArray.total_pd`.  For ultrametric trees the extra
    depths are already present (every child depth is some parent depth or
    ``root_age``), so the curve is unchanged there.

    Note on the active-branch count (documented deviation from treesliceR):
    ``treesliceR::squeeze_root`` and ``phylo_pieces`` build the equal-PD curve
    with ``nBranch = 2:(length(unique(nodes$YearBegin)) + 1)``, i.e. they assume
    a strictly binary tree whose node depths are all distinct, so on a
    multifurcating tree R's cumulative PD does not even reach ``total_pd`` at
    the root age and its "equal PD" slices are not equal.  This implementation
    counts the branches genuinely alive on each depth segment, which makes
    ``cum_pd(root_age) == tree.total_pd()`` hold exactly for any topology.
    The deviation is deliberate, correct, and *not* cross-validated against R
    (every committed reference scenario uses ``criterion="my"``); see
    ``validation/NOTES.md`` and the difference tables in the user manuals.
    """
    depths = tree.node_depths()
    b = depths[tree.edges[:, 0]]
    d = depths[tree.edges[:, 1]]
    T = tree.root_age
    breaks = np.unique(np.concatenate(([0.0], b, d, [T])))
    mids = (breaks[:-1] + breaks[1:]) / 2.0
    # number of branches alive on each segment (binary trees: 2, 3, ..., S)
    active = ((b[None, :] <= mids[:, None]) & (mids[:, None] < d[None, :])).sum(axis=1)
    seg_pd = active * np.diff(breaks)
    cum_pd = np.concatenate(([0.0], np.cumsum(seg_pd)))
    return breaks, cum_pd


def depths_for_pd(tree: TreeArray, target_pd) -> np.ndarray:
    """Vectorised inverse of :func:`pd_depth_curve`.

    One curve is built for the whole batch, which is why :func:`windows_pd`
    does not call :func:`depth_for_pd` once per window: with ``k`` slices on an
    ``E``-edge tree that recomputed the curve ``k`` times (measured: 4.9 s
    versus 0.38 s for a 100-window, 5000-tip ``criterion='pd'`` run).
    """
    breaks, cum_pd = pd_depth_curve(tree)
    total = cum_pd[-1]
    targets = np.atleast_1d(np.asarray(target_pd, dtype=np.float64))
    # the slack is relative to the tree's own size, like every other tolerance
    # here: `total_pd()` sums the edge lengths while `cum_pd[-1]` accumulates
    # segment areas, so on a 5000-edge tree the two differ by ~1e-9 in absolute
    # terms and a fixed 1e-12 would reject the caller's own total PD
    slack = 1e-9 * max(1.0, float(total))
    if np.any(targets < 0.0) or np.any(targets > total + slack):
        raise ValueError(
            f"target PD {targets.min():.6g}..{targets.max():.6g} outside the "
            f"available range [0, {total:.6g}]"
        )
    targets = np.minimum(targets, total)
    atol = 1e-12 * max(1.0, float(total))
    j = np.searchsorted(cum_pd, targets, side="left")
    exact = np.isclose(cum_pd[np.minimum(j, cum_pd.size - 1)], targets, rtol=0, atol=atol)
    jj = np.clip(j, 1, cum_pd.size - 1)  # keeps the interpolation branch in range
    seg_pd = cum_pd[jj] - cum_pd[jj - 1]
    seg_len = breaks[jj] - breaks[jj - 1]
    interp = np.where(
        seg_pd > 0.0, breaks[jj - 1] + (targets - cum_pd[jj - 1]) * seg_len / seg_pd, breaks[jj]
    )
    return np.where(exact, breaks[np.clip(j, 0, breaks.size - 1)].astype(float), interp)


def depth_for_pd(tree: TreeArray, target_pd: float) -> float:
    """Depth-from-root at which the cumulative PD from the root reaches ``target_pd``."""
    return float(depths_for_pd(tree, target_pd)[0])


def _round_pieces_like_r(total: float, width: float) -> tuple[int, float]:
    """treesliceR-compatible rounding of a requested width to (n_slices, actual_width)."""
    if width <= 0:
        raise ValueError("width must be positive")
    if width > total:
        raise ValueError("requested width exceeds the available total")
    expected = total / width
    floor_w = total / np.floor(expected)
    ceil_w = total / np.ceil(expected)
    if width == total:
        return 1, float(total)
    if abs(width - floor_w) < abs(width - ceil_w) and floor_w != total:
        return int(np.floor(expected)), float(floor_w)
    return int(np.ceil(expected)), float(ceil_w)


def windows_time(tree: TreeArray, n: int | None = None, width: float | None = None) -> np.ndarray:
    """Equal-time windows in depth coordinates (root -> tip)."""
    T = tree.root_age
    if n is not None:
        if width is not None:
            raise ValueError("pass either n or width, not both")
        k = int(n)
        if k < 1:
            raise ValueError("n must be a positive integer")
        bounds = np.linspace(0.0, T, k + 1)
    elif width is not None:
        k, _w = _round_pieces_like_r(T, float(width))
        bounds = np.linspace(0.0, T, k + 1)
    else:
        raise ValueError("pass either n or width")
    return np.column_stack([bounds[:-1], bounds[1:]])


def windows_pd(tree: TreeArray, n: int | None = None, width: float | None = None) -> np.ndarray:
    """Equal-PD windows (from the root side) in depth coordinates (root -> tip)."""
    total = tree.total_pd()
    if n is not None:
        if width is not None:
            raise ValueError("pass either n or width, not both")
        k = int(n)
        if k < 1:
            raise ValueError("n must be a positive integer")
    elif width is not None:
        k, _ = _round_pieces_like_r(total, float(width))
    else:
        raise ValueError("pass either n or width")
    thresholds = np.linspace(0.0, total, k + 1)
    depths = depths_for_pd(tree, thresholds)
    return np.column_stack([depths[:-1], depths[1:]])


def make_windows(
    tree: TreeArray,
    n: int | None = None,
    width: float | None = None,
    criterion: str = "time",
) -> np.ndarray:
    """Dispatch window construction by criterion ("time" or "pd")."""
    if criterion == "time":
        if tree.root_age <= 0.0:
            raise ValueError(
                "the tree has zero depth (every branch length is 0), so no time "
                "window can be placed in it"
            )
        return windows_time(tree, n=n, width=width)
    if criterion == "pd":
        if tree.total_pd() <= 0.0:
            raise ValueError(
                "the tree stores no phylogenetic diversity (total PD is 0), so it "
                "cannot be cut into equal-PD windows"
            )
        return windows_pd(tree, n=n, width=width)
    raise ValueError("criterion must be 'time' or 'pd'")
