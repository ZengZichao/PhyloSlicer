"""Robust exponential-accumulation fitting (replaces the single-start nls)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

__all__ = [
    "RateFit",
    "fit_rate",
    "fit_rates_batch",
    "origin_time",
    "adaptive_bounds",
    "DEFAULT_RATE_BOUNDS",
]

#: Fallback search interval for the rate, in units of ``1 / branch-length unit``.
#:
#: ``r`` multiplies slice ages, so the natural scale of a rate is
#: ``1 / root_age``: on a tree normalised to unit age the band ``(1e-6, 50)``
#: spans "essentially no accumulation" to "everything in the youngest 2% of the
#: tree".  On a 300 Myr tree the same numbers mean something else entirely, so
#: pass ``bounds="auto"`` (or an explicit band) for trees whose ages are not
#: normalised - see :func:`adaptive_bounds`.
DEFAULT_RATE_BOUNDS = (1e-6, 50.0)

#: Relative width of the window inside which a solution counts as "at a bound";
#: always combined with the outward-gradient test in :func:`_bound_flags`.
_BOUND_TOL = 1e-4


def _as_bounds(bounds, ages: np.ndarray) -> tuple[float, float]:
    """Resolve ``bounds`` into a concrete ``(lo, hi)`` pair.

    ``"auto"`` rescales :data:`DEFAULT_RATE_BOUNDS` by ``1 / max(ages)``, which
    keeps the band dimensionless (a multiple of the inverse depth scale) and
    reproduces ``(1e-6, 50)`` exactly on unit-age trees.
    """
    if bounds is None or (isinstance(bounds, str) and bounds == "auto"):
        return adaptive_bounds(ages)
    lo, hi = (float(bounds[0]), float(bounds[1]))
    if not (0 < lo < hi):
        raise ValueError(f"bounds must satisfy 0 < lo < hi, got {(lo, hi)}")
    return lo, hi


def adaptive_bounds(ages: np.ndarray) -> tuple[float, float]:
    """``(1e-6, 50) / max(ages)`` - the dimension-preserving default band."""
    ages = np.asarray(ages, dtype=float)
    scale = float(np.nanmax(np.abs(ages))) if ages.size else 0.0
    if not np.isfinite(scale) or scale <= 0:
        return DEFAULT_RATE_BOUNDS
    return (DEFAULT_RATE_BOUNDS[0] / scale, DEFAULT_RATE_BOUNDS[1] / scale)


def _bound_flags(r, lo: float, hi: float, d1_at_bounds=None):
    """Flags for "the reported rate is a boundary optimum, not an interior one".

    The position window is deliberately looser than a solver tolerance (the
    per-row ``least_squares`` stops a few ulp away from the bound), so it is
    always ANDed with the outward-gradient test: being parked at ``lo`` only
    means "clipped" when the objective *rises* there, and likewise at ``hi``
    when it still *falls*.  Without the gradient condition an interior optimum
    that happens to coincide with a bound could not be told apart from a
    solution the bound had swallowed.
    """
    r = np.atleast_1d(np.asarray(r, dtype=float))
    at_lo = np.abs(r - lo) <= _BOUND_TOL * abs(lo)
    at_hi = np.abs(r - hi) <= _BOUND_TOL * abs(hi)
    if d1_at_bounds is not None:
        d1_lo, d1_hi = d1_at_bounds
        at_lo = at_lo & (d1_lo >= 0)
        at_hi = at_hi & (d1_hi <= 0)
    return at_lo, at_hi


def _sse_derivative(cum: np.ndarray, ages: np.ndarray, r) -> np.ndarray:
    """True d SSE / d r = +2 * sum_j ages_j * u_j * (c_j - u_j), u_j = e^(-r a_j)."""
    cum = np.atleast_2d(np.asarray(cum, dtype=float))
    r = np.atleast_1d(np.asarray(r, dtype=float))
    u = np.exp(-r[:, None] * ages[None, :])
    return 2.0 * (ages * (cum - u) * u).sum(-1)


def _bound_reason(at_lo: bool, at_hi: bool, lo: float, hi: float) -> str:
    if at_hi:
        return (
            f"rate clipped at the upper bound {hi:g} of `bounds`; the "
            "unconstrained optimum lies beyond it, so the true accumulation is "
            "faster than reported - widen `bounds` or interpret as a lower bound"
        )
    if at_lo:
        return (
            f"rate clipped at the lower bound {lo:g} of `bounds`; the "
            "profile carries no detectable time structure, so the reported rate "
            "is effectively zero - an upper bound on the true accumulation"
        )
    return ""


@dataclass
class RateFit:
    """Result of fitting ``cumulative(profile) ~ exp(-r * ages)``.

    Attributes
    ----------
    rate : float
        Fitted exponential rate r (NaN when the fit failed; failures carry a
        ``reason`` instead of being silent).  When the row is clipped at a bound
        the bound value itself is returned.
    converged : bool
        ``True`` only when the solver reached a *stationary point of the
        objective in the interior of* ``bounds``.  A boundary optimum - i.e. a
        rate pinned at ``bounds[0]`` or ``bounds[1]`` - is **not** reported as
        convergence; it is flagged by :attr:`at_lower_bound` /
        :attr:`at_upper_bound` and explained in :attr:`reason`.  Clipping and
        convergence are therefore decoupled (the two signals are separate by design).
    r_squared : float
        1 - SSE/SST of the fitted cumulative curve.
    residuals : (k,) ndarray
    n_iter : int
    start_used : float
    bounds : tuple of float
        The ``(lo, hi)`` band the solve was confined to.
    at_lower_bound : bool
        The reported rate equals ``bounds[0]`` (constrained optimum on the
        boundary).
    at_upper_bound : bool
        The reported rate equals ``bounds[1]``.
    reason : str
        Empty when converged; a short explanation otherwise, including which
        bound clipped the estimate.
    """

    rate: float = np.nan
    converged: bool = False
    r_squared: float = np.nan
    residuals: np.ndarray = field(default_factory=lambda: np.empty(0))
    n_iter: int = 0
    start_used: float = np.nan
    bounds: tuple[float, float] = DEFAULT_RATE_BOUNDS
    at_lower_bound: bool = False
    at_upper_bound: bool = False
    reason: str = ""

    def __repr__(self) -> str:  # pragma: no cover - cosmetic
        if self.at_upper_bound:
            status = "clipped at upper bound"
        elif self.at_lower_bound:
            status = "clipped at lower bound"
        elif self.converged:
            status = "ok"
        else:
            status = f"failed ({self.reason})"
        return f"RateFit(rate={self.rate:.6g}, {status})"


def _model(r: float, ages: np.ndarray) -> np.ndarray:
    return np.exp(-r * ages)


def _log_linear_start(cum: np.ndarray, ages: np.ndarray, bounds=DEFAULT_RATE_BOUNDS) -> float:
    """Analytic coarse estimate: slope of log(cum) vs ages (positive cum only)."""
    lo, hi = bounds
    mask = cum > 0
    if mask.sum() < 2:
        return float(np.clip(0.2, lo, hi))
    x = ages[mask]
    y = np.log(cum[mask])
    # guard against -inf from log(very small numbers)
    finite = np.isfinite(y)
    if finite.sum() < 2:
        return float(np.clip(0.2, lo, hi))
    x, y = x[finite], y[finite]
    slope = np.polyfit(x, y, 1)[0]
    r0 = -slope
    if not np.isfinite(r0) or r0 <= 0:
        return float(np.clip(0.2, lo, hi))
    return float(np.clip(r0, lo, hi))


def fit_rate(
    profile: np.ndarray,
    ages: np.ndarray,
    method: str = "lsq",
    multistart: int = 8,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
) -> RateFit:
    """Fit the single-parameter exponential accumulation model.

    Parameters
    ----------
    profile : (k,) array
        Per-slice index values in root -> tip order (the function accumulates
        them internally, matching treesliceR's ``cumsum`` + normalisation).
    ages : (k,) array
        Older boundary of each slice in Ma before present, root -> tip order
        (identical to the reversed ``timeSteps`` vector of treesliceR).
    method : {"lsq", "mle"}
        "lsq" = bounded least squares (Gaussian errors); "mle" = Gaussian
        maximum likelihood (same objective, different optimiser).
    multistart : int
        Number of starting points tried around the analytic coarse estimate;
        the best objective wins.  Unlike ``stats::nls``, a failure on one
        start never aborts the whole analysis.
    bounds : (lo, hi) or "auto"
        Search band for the rate, in ``1 / branch-length unit``.  ``"auto"``
        rescales the default ``(1e-6, 50)`` by ``1 / max(ages)`` so the band
        stays dimensionless on trees whose ages are not normalised.  A solution
        pinned at either end is reported through
        :attr:`RateFit.at_lower_bound` / :attr:`RateFit.at_upper_bound` with an
        explanatory :attr:`~RateFit.reason` and ``converged=False``: clipping
        and convergence are separate signals by design.

    Notes
    -----
    treesliceR hands the problem to ``stats::nls`` *unconstrained*, so this band
    is a deliberate, user-visible difference.  Never read a rate sitting on
    ``hi`` as "the rate is 50": it is only a lower bound on the true
    accumulation rate, and one sitting on ``lo`` only an upper bound.
    """
    profile = np.asarray(profile, dtype=np.float64)
    ages = np.asarray(ages, dtype=np.float64)
    if method not in ("lsq", "mle"):
        # everything else in the package validates its string options; falling
        # through to least squares made method="maximum_likelihood" silently
        # report a Gaussian *least-squares* answer
        raise ValueError(f"method must be 'lsq' or 'mle', got {method!r}")
    if profile.shape != ages.shape:
        raise ValueError("profile and ages must have the same length")
    if profile.size < 2:
        return RateFit(reason="fewer than two slices")
    if not np.isfinite(profile).all() or not np.isfinite(ages).all():
        return RateFit(reason="non-finite input values")
    total = profile.sum()
    if total <= 0:
        return RateFit(reason="profile sums to zero (empty or degenerate assemblage)")

    cum = np.cumsum(profile) / total  # cumulative proportion in (0, 1]
    lo, hi = _as_bounds(bounds, ages)

    r0 = _log_linear_start(cum, ages, (lo, hi))
    grid = np.clip(r0 * np.geomspace(0.1, 10.0, max(1, multistart)), lo, hi)
    grid[0] = r0

    def sse(r: float) -> float:
        res = cum - _model(r, ages)
        return float(res @ res)

    # gradient at the two band edges: decides whether a solution parked on an
    # edge is a genuine boundary optimum (see _bound_flags)
    d1_lo, d1_hi = _sse_derivative(cum, ages, np.array([lo, hi]))

    best = RateFit(reason="no start converged", bounds=(lo, hi))
    best_sse = np.inf
    n_failures = 0
    last_error = ""
    for r_start in grid:
        try:
            if method == "mle":
                from scipy.optimize import minimize

                def nll(r_arr: np.ndarray) -> float:
                    r = float(np.asarray(r_arr).reshape(()))
                    var = max(np.var(cum), 1e-12)
                    res = cum - _model(r, ages)
                    return 0.5 * float(res @ res) / var

                opt = minimize(nll, x0=[r_start], method="Nelder-Mead", options={"maxiter": 1000})
                r_hat = float(np.clip(np.asarray(opt.x).reshape(()), lo, hi))
                n_iter = int(opt.nit)
                ok = bool(opt.success)
            else:
                opt = least_squares(
                    lambda r: cum - _model(np.asarray(r).reshape(()), ages),
                    x0=float(r_start),
                    bounds=(lo, hi),
                    xtol=1e-12,
                    ftol=1e-12,
                    gtol=1e-12,
                    max_nfev=2000,
                )
                r_hat = float(opt.x[0])
                n_iter = int(opt.nfev)
                ok = bool(opt.success)
        except Exception as exc:  # a single bad start must not abort the run
            n_failures += 1
            last_error = str(exc)
            continue
        this_sse = sse(r_hat)
        if this_sse < best_sse:
            at_lo_arr, at_hi_arr = _bound_flags(
                np.asarray([r_hat]), lo, hi, (np.array([d1_lo]), np.array([d1_hi]))
            )
            at_lo, at_hi = bool(at_lo_arr[0]), bool(at_hi_arr[0])
            clipped = at_lo or at_hi
            residuals = cum - _model(r_hat, ages)
            sst = float(((cum - cum.mean()) ** 2).sum())
            if clipped:
                reason = _bound_reason(at_lo, at_hi, lo, hi)
            else:
                reason = "" if ok else "optimiser did not report convergence"
            best = RateFit(
                rate=r_hat,
                # a boundary optimum is a clipped answer, not a convergence
                # success: keep the two signals orthogonal by design
                converged=ok and not clipped,
                r_squared=1.0 - this_sse / sst if sst > 0 else np.nan,
                residuals=residuals,
                n_iter=n_iter,
                start_used=float(r_start),
                bounds=(lo, hi),
                at_lower_bound=at_lo,
                at_upper_bound=at_hi,
                reason=reason,
            )
            best_sse = this_sse
    if not np.isfinite(best.rate) and n_failures:
        best.reason = f"all starts failed ({last_error})"
        return best
    # The constrained optimum may live on an edge of the band, and a bounded
    # solver can stop a little way inside it (for a near-flat objective,
    # noticeably inside).  Rather than trust the reported position, compare the
    # objective at the two bounds explicitly: whichever bound is at least as
    # good as the interior answer *is* the constrained optimum.
    for edge, at_lo_flag, at_hi_flag in ((lo, True, False), (hi, False, True)):
        edge_sse = sse(edge)
        if edge_sse <= best_sse * (1.0 + 1e-12):
            sst = float(((cum - cum.mean()) ** 2).sum())
            best = RateFit(
                rate=float(edge),
                converged=False,
                r_squared=1.0 - edge_sse / sst if sst > 0 else np.nan,
                residuals=cum - _model(float(edge), ages),
                n_iter=0,
                start_used=float(edge),
                bounds=(lo, hi),
                at_lower_bound=at_lo_flag,
                at_upper_bound=at_hi_flag,
                reason=_bound_reason(at_lo_flag, at_hi_flag, lo, hi),
            )
            best_sse = edge_sse
    return best


def origin_time(rate: float, percentile: float = 0.95) -> float:
    """Time (Ma) before the present at which the fitted curve stands at ``1 - percentile``.

    ``pXO = -ln(1 - percentile) / rate``.  Under the model of eqn 3 the fitted curve
    is ``cumulative = exp(-rate * age)``, so this is the moment by which
    ``1 - percentile`` of the extant total had already accumulated and the remaining
    ``percentile`` was still to come.  ``percentile=0.95`` therefore returns
    ``-ln(0.05)/rate`` - the same quantity treesliceR reaches with its ``pDO = 5``
    default, which is likewise read as "the first 5% accumulated".
    """
    if not 0 < percentile < 1:
        raise ValueError("percentile must lie strictly between 0 and 1")
    if not np.isfinite(rate) or rate <= 0:
        return np.nan
    return float(-np.log(1.0 - percentile) / rate)


def fit_rates_batch(
    profiles: np.ndarray,
    ages: np.ndarray,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
    max_iter: int = 100,
    tol: float = 1e-12,
) -> list[RateFit]:
    """Fit every row of an (n_sites, k) profile matrix with one vectorised
    Newton solve per iteration.

    All rows share the same ``ages`` grid, so the 1-D bounded objective is
    solved simultaneously for every site: cost O(iterations x n_sites x k)
    with no per-site optimiser call.  Rows where the vectorised Newton stencil
    fails (non-finite or degenerate profile) are returned unconverged and can be
    refined with :func:`fit_rate`.

    ``bounds`` accepts an explicit ``(lo, hi)`` band or ``"auto"`` (rescale the
    default by ``1 / max(ages)``).  Rows whose constrained optimum sits on a
    bound return that bound as the rate, set :attr:`RateFit.at_lower_bound` /
    :attr:`RateFit.at_upper_bound`, and are *not* certified as converged.
    """
    profiles = np.atleast_2d(np.asarray(profiles, dtype=np.float64))
    ages = np.asarray(ages, dtype=np.float64)
    S, k = profiles.shape
    lo, hi = _as_bounds(bounds, ages)
    totals = profiles.sum(axis=1, keepdims=True)
    finite_rows = np.isfinite(profiles).all(axis=1) & (totals[:, 0] > 0)
    cum = np.cumsum(profiles, axis=1) / np.where(totals > 0, totals, 1.0)

    def sse_grad_hess(rv: np.ndarray):
        # SSE(r) = sum_j (c_j - u_j)^2 with u_j = exp(-r * a_j)
        # SSE'  =  2 sum a_j u_j (c_j - u_j)
        # SSE'' =  2 sum a_j^2 u_j (2 u_j - c_j)
        # NB: the sign of SSE' matters.  An earlier revision returned
        # -2*sum(...) here, i.e. the *negative* gradient; the bisection (which
        # only looks for a sign change) and the |SSE'| <= tol criterion are
        # invariant to that, but the stage-3 Newton update `r - d1/d2` was then
        # an *ascent* step, harmless only because the bracket it is clipped to
        # had already collapsed onto the bisection root.  It is now the true
        # gradient, so the polish descends and the bound tests below can read
        # the direction literally.
        rv = np.atleast_1d(np.asarray(rv, dtype=float))
        if rv.size == 1:
            u = np.exp(-rv[0] * ages)[None, :]
        else:
            u = np.exp(-rv[:, None] * ages[None, :])
        res = cum - u
        sse = (res * res).sum(-1)
        d1 = 2.0 * (ages * res * u).sum(-1)
        d2 = 2.0 * (ages**2 * u * (2.0 * u - cum)).sum(-1)
        return sse, d1, d2

    # ---- stage 1: coarse log-grid search (robust global bracket) ------- #
    # The SSE grid comes from the optional-acceleration layer, whose compiled
    # and pure-numpy back-ends both return the full (S, n_grid) matrix.  The
    # solver below is therefore a single code path whether or not numba is
    # installed (pip install phyloslicer[accel]).
    grid = np.geomspace(lo, hi, 96)
    from .._accel import batched_grid_fit

    _, sse_grid = batched_grid_fit(profiles, ages, grid)
    best = sse_grid.argmin(axis=1)
    bracket = np.column_stack(
        [grid[np.maximum(best - 1, 0)], grid[np.minimum(best + 1, grid.size - 1)]]
    )

    # ---- stage 2: vectorised bisection inside each bracket ------------- #
    # d1 changes sign at an interior optimum; bisection on its sign is
    # unconditionally stable, and 80 halvings exhaust float precision.  A
    # bracket with *no* sign change means the constrained optimum is at a
    # bracket edge, so bisecting there is meaningless (it used to drift to
    # grid[best + 1] and report that neighbour as the answer when the true
    # minimiser lay below bounds[0]).
    a_br = bracket[:, 0].copy()
    b_br = bracket[:, 1].copy()
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        f_a = sse_grad_hess(a_br)[1]
        f_b = sse_grad_hess(b_br)[1]
        has_stationary = np.isfinite(f_a) & np.isfinite(f_b) & (f_a * f_b <= 0)
        for _ in range(max_iter):
            mid = 0.5 * (a_br + b_br)
            f_mid = sse_grad_hess(mid)[1]
            same_as_a = f_a * f_mid <= 0
            b_br = np.where(same_as_a, mid, b_br)
            a_br = np.where(same_as_a, a_br, mid)
            # keep the sign stored for the left endpoint in sync with the
            # endpoint actually kept, otherwise the bracket can be lost
            f_a = np.where(same_as_a, f_a, f_mid)
            if np.all(np.abs(b_br - a_br) <= tol * (1.0 + np.abs(mid))):
                break
        r_bisected = 0.5 * (a_br + b_br)

        # rows without a sign change: keep the bracket endpoint with the lower
        # objective, which is exactly the bound when the grid argmin was an
        # endpoint (best == 0 -> a_br == lo; best == last -> b_br == hi)
        sse_a = sse_grad_hess(bracket[:, 0])[0]
        sse_b = sse_grad_hess(bracket[:, 1])[0]
        r = np.where(
            has_stationary, r_bisected, np.where(sse_a <= sse_b, bracket[:, 0], bracket[:, 1])
        )
        r = np.clip(r, lo, hi)

        # ---- stage 3: a few clipped Newton steps for fast polish -------- #
        # only rows with a genuine stationary point may be polished; snapping a
        # boundary row back onto its bound keeps the reported value exact
        for _ in range(10):
            _, d1, d2 = sse_grad_hess(r)
            step = np.where((d2 > 0) & np.isfinite(d2), d1 / np.where(d2 > 0, d2, 1.0), 0.0)
            r_new = np.where(has_stationary, np.clip(r - step, bracket[:, 0], bracket[:, 1]), r)
            if np.all(np.abs(r_new - r) <= 1e-14 * (1.0 + np.abs(r))):
                break
            r = r_new
        sse, d1, d2 = sse_grad_hess(r)

    # Convergence is certified only for rows that landed on a genuine
    # *interior* stationary point (vanishing gradient).  Hitting a bound is a
    # separate, explicitly reported condition and must never be laundered into
    # `converged=True`.  Comparing the gradient against the
    # post-bisection micro-bracket would be vacuous - every r sits "at the
    # edge" of a bracket that has collapsed onto it - hence the comparison to
    # the objective scale below.
    at_lo, at_hi = _bound_flags(
        r, lo, hi, (sse_grad_hess(np.array([lo]))[1], sse_grad_hess(np.array([hi]))[1])
    )
    scale = np.abs(cum).max(axis=1) + 1e-300
    small_grad = np.abs(d1) <= 1e-6 * scale * (np.abs(ages).max() + 1.0)
    interior_ok = small_grad & ~at_lo & ~at_hi
    converged = finite_rows & np.isfinite(sse) & interior_ok & (totals[:, 0] > 0)

    sst = ((cum - cum.mean(axis=1, keepdims=True)) ** 2).sum(1)
    out = []
    for i in range(S):
        if not finite_rows[i]:
            out.append(RateFit(reason="non-finite or empty profile", bounds=(lo, hi)))
            continue
        residuals = cum[i] - np.exp(-r[i] * ages)
        if at_lo[i] or at_hi[i]:
            reason = _bound_reason(bool(at_lo[i]), bool(at_hi[i]), lo, hi)
        elif converged[i]:
            reason = ""
        else:
            reason = "vectorised Newton did not satisfy criteria"
        out.append(
            RateFit(
                rate=float(r[i]),
                converged=bool(converged[i]),
                r_squared=float(1.0 - sse[i] / sst[i]) if sst[i] > 0 else np.nan,
                residuals=residuals,
                n_iter=max_iter,
                start_used=float(r[i]),
                bounds=(lo, hi),
                at_lower_bound=bool(at_lo[i]),
                at_upper_bound=bool(at_hi[i]),
                reason=reason,
            )
        )
    return out
