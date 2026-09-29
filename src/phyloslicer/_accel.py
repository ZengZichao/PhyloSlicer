"""Optional acceleration layer (``pip install phyloslicer[accel]``).

The batched coarse grid search of :mod:`phyloslicer.rates.fitting` is the only
hot loop that is compiled.  Both back-ends compute *the same* per-site SSE
matrix of shape ``(n_sites, n_grid)`` and return it, so the caller's
downstream logic (``argmin``, bracketing, bisection, Newton polish) is a
single code path that cannot diverge between an accelerated and a plain
install.  Importing this module never fails, and it never requires numba.

Historical note: an earlier version returned ``(best_rate, best_sse)`` from
the compiled branch while the numpy branch returned a full SSE grid, so any
install that could import numba raised ``AxisError`` one level up.  The
back-ends are now required to be shape- and value-identical, which
``tests/test_rates.py`` asserts whenever numba is importable.
"""

from __future__ import annotations

import numpy as np

try:  # pragma: no cover - depends on environment
    from numba import njit, prange

    HAS_NUMBA = True
except ImportError:  # pragma: no cover
    HAS_NUMBA = False


if HAS_NUMBA:  # pragma: no cover - compiled path

    @njit(parallel=True, cache=True, fastmath=False)
    def _fill_sse_grid(cum, ages, rates_grid, sse):  # pragma: no cover
        """Fill ``sse[i, g] = sum_j (cum[i, j] - exp(-r_g * ages[j]))**2``."""
        n_sites, k = cum.shape
        n_grid = rates_grid.shape[0]
        for i in prange(n_sites):
            for g in range(n_grid):
                r = rates_grid[g]
                acc = 0.0
                for j in range(k):
                    res = cum[i, j] - np.exp(-r * ages[j])
                    acc += res * res
                sse[i, g] = acc
        return sse

    def _sse_grid(cum: np.ndarray, ages: np.ndarray, grid: np.ndarray) -> np.ndarray:
        sse = np.empty((cum.shape[0], grid.size), dtype=np.float64)
        _fill_sse_grid(cum, ages, grid, sse)
        return sse

else:  # pragma: no cover - pure-numpy path

    def _sse_grid(cum: np.ndarray, ages: np.ndarray, grid: np.ndarray) -> np.ndarray:
        # looping over the grid (rather than a 3-D broadcast) bounds memory.
        sse = np.empty((cum.shape[0], grid.size), dtype=np.float64)
        for gi, gv in enumerate(grid):
            res = cum - np.exp(-gv * ages)[None, :]
            sse[:, gi] = (res * res).sum(axis=1)
        return sse


def normalised_cumulative(profiles: np.ndarray) -> np.ndarray:
    """Row-normalised cumulative profile of a ``(n_sites, k)`` array.

    One-dimensional input (a single site) is promoted to one row, so both
    back-ends see the same 2-D shape.  Rows that sum to zero stay zero rather
    than becoming NaN; non-finite rows are passed through untouched and are
    filtered by the caller.
    """
    profiles = np.atleast_2d(np.asarray(profiles, dtype=np.float64))
    totals = profiles.sum(axis=1, keepdims=True)
    return np.cumsum(profiles, axis=1) / np.where(totals > 0, totals, 1.0)


def batched_grid_fit(profiles: np.ndarray, ages: np.ndarray, rates_grid: np.ndarray):
    """SSE of the one-parameter model on a rate grid, per site.

    Parameters
    ----------
    profiles
        ``(n_sites, k)`` (or ``(k,)`` for a single site) index values per slice.
    ages
        ``(k,)`` older boundary of each slice.
    rates_grid
        ``(n_grid,)`` candidate rates.

    Returns
    -------
    (cumulative, sse_grid)
        ``cumulative`` is the ``(n_sites, k)`` normalised cumulative profile
        that was searched, and ``sse_grid`` the ``(n_sites, n_grid)`` sum of
        squared errors.  The shape is identical whichever back-end is active.
    """
    cum = normalised_cumulative(profiles)
    ages = np.asarray(ages, dtype=np.float64)
    grid = np.asarray(rates_grid, dtype=np.float64)
    if cum.shape[1] != ages.size:
        raise ValueError(f"profile width {cum.shape[1]} does not match {ages.size} slice ages")
    return cum, _sse_grid(cum, ages, grid)


__all__ = ["HAS_NUMBA", "batched_grid_fit", "normalised_cumulative"]
