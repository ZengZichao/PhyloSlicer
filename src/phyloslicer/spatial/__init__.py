"""Spatial helpers: adjacency construction and rate mapping.

``geopandas`` is optional - everything works with plain DataFrames carrying
rectangle bounds (``minx/miny/maxx/maxy``) or centroid coordinates (``x/y``).
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

__all__ = [
    "adjacency_from_grid",
    "adjacency_from_coords",
    "adjacency_from_bounds",
    "adjacency_criterion",
    "resolve_adjacency_tol",
    "to_adjacency_dict",
    "from_adjacency_dict",
    "rate_map",
]


def _cell_scale(bounds: np.ndarray) -> float:
    """Median cell edge length - the natural length unit of the grid."""
    width = bounds[:, 2] - bounds[:, 0]
    height = bounds[:, 3] - bounds[:, 1]
    edges = np.concatenate([width, height])
    positive = edges[edges > 0]
    if positive.size == 0:
        return 1.0  # degenerate point grid: fall back to a unit tolerance scale
    return float(np.median(positive))


def resolve_adjacency_tol(
    bounds: np.ndarray, rel_tol: float = 1e-9, tol: float | None = None
) -> float:
    """Contact tolerance for ``bounds``: ``tol`` if given, else relative to cells.

    ``adjacency_from_bounds`` / :func:`adjacency_from_grid` used to compare raw
    coordinate differences against a fixed absolute ``1e-9``.  On a
    degree-unit 0.5 deg grid that judged a 1e-10 deg (~11 um) genuine gap as
    "touching" while rejecting 1e-8 residuals from a projection or
    trigonometric transform as "not touching" - silently changing the beta
    neighbourhood set, and hence every CpB value.  The default now scales with
    the median cell edge, so the criterion is invariant to the coordinate unit.
    """
    if tol is not None:
        if not np.isfinite(tol) or tol < 0:
            raise ValueError("tol must be a finite non-negative number")
        return float(tol)
    if not np.isfinite(rel_tol) or rel_tol < 0:
        raise ValueError("rel_tol must be a finite non-negative number")
    bounds = np.asarray(bounds, dtype=float)
    return float(rel_tol) * _cell_scale(bounds)


def adjacency_from_bounds(
    bounds: np.ndarray,
    method: str = "queen",
    rel_tol: float = 1e-9,
    tol: float | None = None,
) -> np.ndarray:
    """Adjacency of axis-aligned rectangles given (n, 4) [minx, miny, maxx, maxy].

    ``queen`` includes corner contacts; ``rook`` requires a shared edge of
    positive length.  Two rectangles "touch" when their intervals meet to within
    ``tol``, which defaults to ``rel_tol * median(cell edge length)``; see
    :func:`resolve_adjacency_tol`.  Use :func:`adjacency_criterion` to record
    the criterion actually applied (numpy arrays cannot carry metadata).
    """
    if method not in ("rook", "queen"):
        raise ValueError("method must be 'rook' or 'queen'")
    bounds = np.asarray(bounds, dtype=float)
    if bounds.ndim != 2 or bounds.shape[1] != 4:
        raise ValueError("bounds must be an (n, 4) [minx, miny, maxx, maxy] array")
    if bounds.shape[0] == 0:
        raise ValueError("bounds is empty; at least one site is required")
    if (bounds[:, 2] < bounds[:, 0]).any() or (bounds[:, 3] < bounds[:, 1]).any():
        raise ValueError("bounds rows need minx <= maxx and miny <= maxy")
    tol = resolve_adjacency_tol(bounds, rel_tol, tol)
    minx, miny, maxx, maxy = bounds[:, 0], bounds[:, 1], bounds[:, 2], bounds[:, 3]
    ox = np.minimum(maxx[:, None], maxx[None, :]) - np.maximum(minx[:, None], minx[None, :])
    oy = np.minimum(maxy[:, None], maxy[None, :]) - np.maximum(miny[:, None], miny[None, :])
    overlap_x = np.maximum(ox, 0.0)
    overlap_y = np.maximum(oy, 0.0)
    touch_x = np.abs(ox) <= tol  # x-intervals meet at a vertical line
    touch_y = np.abs(oy) <= tol
    edge_contact = ((overlap_x > tol) & touch_y) | ((overlap_y > tol) & touch_x)
    corner_contact = touch_x & touch_y
    overlap_contact = (overlap_x > tol) & (overlap_y > tol)
    if method == "queen":
        adj = edge_contact | corner_contact | overlap_contact
    else:
        adj = edge_contact | overlap_contact
    np.fill_diagonal(adj, True)  # focal sites include themselves (treesliceR convention)
    return adj.astype(np.int8)


def adjacency_criterion(
    bounds: np.ndarray, method: str = "queen", rel_tol: float = 1e-9, tol: float | None = None
) -> dict:
    """Audit summary of the criterion :func:`adjacency_from_bounds` would apply.

    Returns the resolved tolerance, the cell scale it was derived from, the
    method, and the resulting number of adjacency edges - enough to reproduce or
    challenge a neighbourhood set after the fact.
    """
    bounds = np.asarray(bounds, dtype=float)
    applied = resolve_adjacency_tol(bounds, rel_tol, tol)
    adj = adjacency_from_bounds(bounds, method=method, tol=applied)
    n = len(bounds)
    return {
        "method": method,
        "tol": applied,
        "rel_tol": float(rel_tol),
        "cell_scale": _cell_scale(bounds),
        "n_sites": n,
        "n_neighbour_pairs": int((adj.sum() - n) // 2),
        "mean_neighbours_per_focal": float(adj.sum(axis=1).mean()),
    }


def adjacency_from_grid(
    gdf: pd.DataFrame,
    method: str = "queen",
    id_col: str | None = None,
    rel_tol: float = 1e-9,
    tol: float | None = None,
) -> np.ndarray:
    """Adjacency matrix from a grid of polygons.

    With geopandas installed, true polygon semantics are used: ``queen``
    counts any shared boundary contact (edges or corners), ``rook`` only
    contacts along a boundary segment of positive length.  Otherwise the
    rectangle bounds columns (``minx/miny/maxx/maxy``) are used with the same
    semantics.

    ``tol`` / ``rel_tol`` set the numeric contact tolerance; without an explicit
    ``tol`` it is scaled by the median cell edge length (see
    :func:`resolve_adjacency_tol`), so the criterion does not depend on whether
    the grid is measured in degrees, metres or kilometres.
    """
    if method not in ("rook", "queen"):
        raise ValueError("method must be 'rook' or 'queen'")
    if hasattr(gdf, "geometry") and getattr(gdf, "geometry", None) is not None:
        try:
            import geopandas as gpd  # noqa: F401

            n = len(gdf)
            if n == 0:
                raise ValueError("the grid is empty")
            geoms = list(gdf.geometry.values)
            # median cell edge, so the contact tolerance is relative to the
            # grid's own coordinate scale (see resolve_adjacency_tol)
            rectangles = np.asarray([list(g.bounds) for g in geoms], dtype=float)
            applied = resolve_adjacency_tol(rectangles, rel_tol, tol)
            adj = np.zeros((n, n), dtype=np.int8)
            for i in range(n):
                if method == "queen":
                    # tolerate numeric noise at shared edges/corners
                    adj[i] = np.asarray(geoms[i].buffer(applied).intersects(geoms), dtype=np.int8)
                else:
                    # rook: the shared boundary must be a segment, not a point
                    boundary = geoms[i].boundary
                    adj[i] = np.asarray(
                        [boundary.intersection(g).length > applied for g in geoms],
                        dtype=np.int8,
                    )
            np.fill_diagonal(adj, 1)
            return adj
        except ImportError:
            warnings.warn("geopandas not available; falling back to bounds-based adjacency")
    for col in ("minx", "miny", "maxx", "maxy"):
        if col not in gdf.columns:
            raise ValueError(
                "grid must either be a GeoDataFrame or provide bounds columns minx/miny/maxx/maxy"
            )
    bounds = gdf[["minx", "miny", "maxx", "maxy"]].to_numpy(dtype=float)
    return adjacency_from_bounds(bounds, method=method, rel_tol=rel_tol, tol=tol)


def adjacency_from_coords(
    x: np.ndarray,
    y: np.ndarray,
    method: str = "knn",
    k: int = 8,
    r: float | None = None,
) -> np.ndarray:
    """Adjacency from coordinates: ``knn`` (k nearest neighbours) or ``radius``."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 1:
        raise ValueError("x and y must be 1-D arrays of equal length")
    n = x.size
    d = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :])
    np.fill_diagonal(d, np.inf)
    if method == "knn":
        if k < 1 or k >= n:
            raise ValueError(f"k must lie in [1, {n - 1}]")
        idx = np.argsort(d, axis=1)[:, :k]
        adj = np.zeros((n, n), dtype=np.int8)
        rows = np.repeat(np.arange(n), k)
        adj[rows, idx.ravel()] = 1
        adj = np.maximum(adj, adj.T)  # symmetrise
    elif method == "radius":
        if r is None or r <= 0:
            raise ValueError("radius method requires r > 0")
        adj = (d <= r).astype(np.int8)
    else:
        raise ValueError("method must be 'knn' or 'radius'")
    np.fill_diagonal(adj, 1)
    return adj


def to_adjacency_dict(adj: np.ndarray, site_ids=None) -> dict:
    """Sparse dict form (site -> [neighbour sites]) for interoperability.

    Returns a plain ``{str(site_id): [neighbour site_ids]}`` dict; keys are
    always stringified so JSON / YAML serialisation works.  Use
    :func:`from_adjacency_dict` with the same ``site_ids`` to rebuild the
    dense matrix (the original id order is recovered via the ``site_ids``
    list, not by ``int()``-parsing the keys).

    Round-trip contract, which the previous docstring left implicit: the
    diagonal is dropped here and unconditionally restored by
    :func:`from_adjacency_dict`, and the reconstruction always sets both
    directions of a listed neighbour pair when ``symmetrise=True``.  A
    round-trip is therefore only exact for matrices that are reflexive
    (ones on the diagonal) and symmetric; an asymmetric input comes back
    symmetrised, which changes the focal neighbourhoods.  The functions check
    and say so rather than quietly altering the result.
    """
    adj = np.asarray(adj)
    if adj.ndim != 2 or adj.shape[0] != adj.shape[1]:
        raise ValueError("adj must be a square (n, n) matrix")
    if not np.array_equal(adj, adj.T):
        warnings.warn(
            "the adjacency matrix is asymmetric: from_adjacency_dict() mirrors "
            "every listed link by default (symmetrise=True), so the rebuilt "
            "matrix is NOT identical to this one and the focal neighbourhoods "
            "change.  Pass symmetrise=False to keep it exactly as written.",
            stacklevel=2,
        )
    n = adj.shape[0]
    ids = list(range(n)) if site_ids is None else list(site_ids)
    if len(ids) != n:
        raise ValueError(f"site_ids has {len(ids)} entries but adj is {n} x {n}")
    return {str(ids[i]): [ids[j] for j in np.flatnonzero(adj[i] == 1) if j != i] for i in range(n)}


def from_adjacency_dict(
    d: dict,
    site_ids=None,
    n_sites: int | None = None,
    symmetrise: bool = True,
) -> np.ndarray:
    """Rebuild a dense adjacency matrix from :func:`to_adjacency_dict` output.

    When ``site_ids`` is provided, keys/values are looked up in the id list
    (supporting non-integer ids); otherwise keys are parsed as ``int``.
    Unknown ids raise :class:`KeyError` with the offending id, the number of
    known ids and up to five examples, instead of the bare
    ``KeyError('some_string')`` the dict comprehension used to leak.

    The diagonal is always filled (a focal neighbourhood includes the focal
    site, treesliceR convention).  ``symmetrise=True`` (default) mirrors every
    listed neighbour link, because the "focal row" convention only makes sense
    for a symmetric neighbourhood relation; pass ``False`` to keep the matrix
    exactly as written.
    """
    keys = list(d.keys())
    if site_ids is not None:
        id_list = list(site_ids)
        str_to_idx = {str(v): i for i, v in enumerate(id_list)}
        n = n_sites if n_sites is not None else len(id_list)
        adj = np.zeros((n, n), dtype=np.int8)

        def lookup(token, *, where):
            key = str(token)
            if key not in str_to_idx:
                raise KeyError(
                    f"adjacency {where} {token!r} is not among the {len(id_list)} "
                    f"site_ids (first ids: {[str(v) for v in id_list[:5]]}"
                    f"{'...' if len(id_list) > 5 else ''}); pass the site_ids that "
                    "were used with to_adjacency_dict, or n_sites if they are ints"
                )
            idx = str_to_idx[key]
            if idx >= n:
                raise KeyError(
                    f"adjacency {where} {token!r} resolves to index {idx}, outside "
                    f"the requested {n} x {n} matrix"
                )
            return idx

        for key, neigh in d.items():
            i = lookup(key, where="key")
            adj[i, i] = 1
            for nb in neigh:
                j = lookup(nb, where="neighbour")
                adj[i, j] = 1
                if symmetrise:
                    adj[j, i] = 1
        return adj
    n = n_sites if n_sites is not None else len(keys)

    def parse(token, *, where):
        try:
            idx = int(token)
        except (TypeError, ValueError) as exc:
            raise KeyError(
                f"adjacency {where} {token!r} is not an integer site index; pass the "
                "site_ids used with to_adjacency_dict for non-integer identifiers"
            ) from exc
        if not 0 <= idx < n:
            raise KeyError(
                f"adjacency {where} {token!r} indexes site {idx}, outside the "
                f"{n} x {n} matrix implied by n_sites / the number of keys"
            )
        return idx

    adj = np.zeros((n, n), dtype=np.int8)
    for key, neigh in d.items():
        i = parse(key, where="key")
        adj[i, i] = 1
        for j in neigh:
            jj = parse(j, where="neighbour")
            adj[i, jj] = 1
            if symmetrise:
                adj[jj, i] = 1
    return adj


def rate_map(
    rates_df: pd.DataFrame,
    grid: pd.DataFrame,
    rate: str = "CpD",
    quantiles: bool = False,
    palette: str = "viridis",
    ax=None,
):
    """Choropleth of a rate column over a grid (mirrors ``CpR_graph(map=...)``).

    ``grid`` must carry either polygon geometry (GeoDataFrame) or bounds
    columns (``minx/miny/maxx/maxy``); the rate column values from
    ``rates_df`` are plotted in grid row order (not by grid index value).
    """
    import matplotlib.pyplot as plt
    from matplotlib import colormaps

    def _get_cmap(name: str):
        # plt.get_cmap is deprecated/removed in newer matplotlib; the
        # colormaps registry is the supported interface (matplotlib >= 3.6)
        try:
            return colormaps[name]
        except KeyError as exc:
            raise ValueError(f"unknown palette {name!r}") from exc

    if rate not in rates_df.columns:
        raise ValueError(f"rate column {rate!r} not found in rates_df")
    values = rates_df[rate].to_numpy(dtype=float)
    if len(values) != len(grid):
        raise ValueError("rates_df and grid must have the same number of rows")

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    if hasattr(grid, "geometry") and getattr(grid, "geometry", None) is not None:
        try:
            cmap = _get_cmap(palette)
            if quantiles:
                levels = np.nanquantile(values, np.linspace(0, 1, 6))
                binned = np.digitize(values, levels[1:-1])
                vrange = (0, max(1, binned.max()))
            else:
                vrange = (np.nanmin(values), np.nanmax(values))
                binned = None
            norm = plt.Normalize(*vrange)
            for poly, val, binv in zip(
                grid.geometry, values, binned if binned is not None else values
            ):
                color = cmap(norm(binv)) if binned is not None else cmap(norm(val))
                if not np.isfinite(val):
                    color = (0.85, 0.85, 0.85, 1.0)
                # Polygon -> exterior ring; MultiPolygon -> every part
                parts = (
                    [poly]
                    if hasattr(poly, "exterior") and poly.exterior is not None
                    else list(getattr(poly, "geoms", []))
                )
                for part in parts:
                    if getattr(part, "exterior", None) is None:
                        continue
                    xs, ys = part.exterior.xy
                    ax.fill(xs, ys, color=color, linewidth=0.2, edgecolor="0.4")
        except ImportError:
            raise
    else:
        if not {"minx", "miny", "maxx", "maxy"}.issubset(grid.columns):
            raise ValueError("grid must be a GeoDataFrame or provide bounds columns")
        cmap = _get_cmap(palette)
        finite = values[np.isfinite(values)]
        vmin, vmax = (np.nanmin(finite), np.nanmax(finite)) if finite.size else (0, 1)
        if quantiles:
            edges = np.nanquantile(finite, np.linspace(0, 1, 6))
            codes = np.searchsorted(edges[1:-1], values)
            norm = plt.Normalize(0, max(1, codes.max()))
        else:
            codes = values
            norm = plt.Normalize(vmin, vmax)
        for (_, row), val, code in zip(grid.iterrows(), values, codes):
            color = cmap(norm(code)) if np.isfinite(val) else (0.85, 0.85, 0.85, 1.0)
            ax.fill(
                [row["minx"], row["maxx"], row["maxx"], row["minx"]],
                [row["miny"], row["miny"], row["maxy"], row["maxy"]],
                color=color,
                linewidth=0.2,
                edgecolor="0.4",
            )
    sm = plt.cm.ScalarMappable(cmap=_get_cmap(palette), norm=norm)
    fig = ax.figure
    fig.colorbar(sm, ax=ax, shrink=0.7, label=rate)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    return ax
