"""Cumulative-rate API: cpd_rate / cpe_rate / cpb_rate / rate_profile / sensitivity."""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..core.tree import TreeArray
from ..indices.api import pb_per_slice, pd_per_slice, pe_per_slice
from ..io.matrices import align_tree_matrix, as_matrix
from ..slicing.api import slice_pieces
from .fitting import DEFAULT_RATE_BOUNDS, RateFit, fit_rate, origin_time

__all__ = [
    "cpd_rate",
    "cpe_rate",
    "cpb_rate",
    "rate_profile",
    "sensitivity",
    "RateFit",
    "fit_rate",
    "origin_time",
    "DEFAULT_RATE_BOUNDS",
]


def _validate_common(tree: TreeArray, mat, n_slices: int, percentile: float):
    """Align tree and matrix (must happen *before* slicing) and check args.

    The slicing options (``criterion``, ``allow_non_ultrametric``) are **not**
    taken here: they belong to ``slice_pieces``, and accepting them in this
    helper once made a caller think a flag had been checked when it was dropped.
    """
    if n_slices < 1:
        raise ValueError("n_slices must be a positive integer")
    if not 0 < percentile < 1:
        raise ValueError("percentile must lie strictly between 0 and 1")
    mat = as_matrix(mat)
    return align_tree_matrix(tree, mat)


def _fit_profiles(profiles: np.ndarray, ages: np.ndarray, multistart: int, bounds) -> list[RateFit]:
    """Fit every row of an (n_sites, k) profile matrix.

    The vectorised Newton kernel (:func:`fit_rates_batch`) solves all sites
    simultaneously; rows it cannot certify are re-fit with the robust
    per-row bounded multi-start optimiser so the two paths share the same
    objective and the same optimum.  A row the batch kernel pinned at a bound
    is *not* re-refit: the per-row optimiser reaches the same boundary answer,
    and re-running it there would replace the exact bound value with a nearby
    interior point (the two signals are deliberately decoupled).
    """
    from .fitting import fit_rates_batch

    batch = fit_rates_batch(profiles, ages, bounds=bounds)
    out: list[RateFit] = []
    for i, fit in enumerate(batch):
        if fit.converged or fit.at_lower_bound or fit.at_upper_bound:
            out.append(fit)
        else:  # rare: refine with the robust per-row optimiser
            out.append(
                fit_rate(profiles[i], ages, method="lsq", multistart=multistart, bounds=bounds)
            )
    return out


def _rate_frame(
    fits: list[RateFit],
    totals: np.ndarray,
    percentile: float,
    rate_name: str,
    total_name: str,
    origin_name: str,
    site_ids,
    ultrametric: bool = True,
    root_age: float | None = None,
) -> pd.DataFrame:
    rows = []
    for sid, fit, tot in zip(site_ids, fits, totals):
        pXO = origin_time(fit.rate, percentile=percentile)
        within, reason = _origin_within_tree(pXO, root_age, fit.reason)
        rows.append(
            {
                "site_id": sid,
                rate_name: fit.rate,
                total_name: tot,
                origin_name: pXO,
                "converged": fit.converged,
                "r_squared": fit.r_squared,
                "reason": reason,
                "at_lower_bound": fit.at_lower_bound,
                "at_upper_bound": fit.at_upper_bound,
                "ultrametric": ultrametric,
                "origin_within_tree": within,
            }
        )
    return pd.DataFrame(rows)


def _origin_within_tree(pXO, root_age, reason):
    """Is the fitted origin inside the time the tree actually covers?

    ``pXO = -ln(1 - percentile) / rate`` is an extrapolation of the fitted
    exponential, and a slow accumulation rate pushes it past the root: the
    origin then names a moment the phylogeny does not contain.  The value is
    reported unchanged - clipping it would hide the model failure - but the row
    says so.
    """
    if root_age is None or not np.isfinite(pXO):
        return True, reason
    if pXO <= root_age:
        return True, reason
    note = "origin extrapolated beyond root"
    return False, note if not reason else f"{reason}; {note}"


def cpd_rate(
    tree: TreeArray,
    mat,
    n_slices: int = 100,
    criterion: str = "time",
    percentile: float = 0.95,
    diagnostics: bool = True,
    multistart: int = 8,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
    allow_non_ultrametric: bool = False,
) -> pd.DataFrame:
    """Rate of PD accumulation (mirrors treesliceR's ``CpD``) with fit diagnostics.

    Returns a DataFrame: ``site_id, CpD, PD, pDO, converged, r_squared, reason,
    at_lower_bound, at_upper_bound, ultrametric``.
    ``percentile=0.95`` (accumulated proportion) matches R's default
    ``pDO = 5``; ``pDO = -ln(1 - percentile) / CpD``.

    ``bounds`` is the search band for the rate, in ``1 / branch-length unit``
    (``"auto"`` rescales ``(1e-6, 50)`` by ``1 / root_age``).  A row pinned at
    either end carries ``at_lower_bound`` / ``at_upper_bound``, an explanatory
    ``reason`` and ``converged=False`` - never a silent "converged" rate of 50.
    treesliceR fits unconstrained, so widen the band before concluding that a
    fast rate is real.  Non-ultrametric trees raise unless
    ``allow_non_ultrametric=True``; the verdict is echoed in ``ultrametric``.
    """
    tree, mat = _validate_common(tree, mat, n_slices, percentile)
    slices = slice_pieces(
        tree, n=n_slices, criterion=criterion, allow_non_ultrametric=allow_non_ultrametric
    )
    profiles = pd_per_slice(tree, mat, slices)
    totals = profiles.sum(axis=1)
    ages = slices.ages
    fits = _fit_profiles(profiles, ages, multistart, bounds)
    if not diagnostics:
        for f in fits:
            f.residuals = np.empty(0)
    return _rate_frame(
        fits,
        totals,
        percentile,
        "CpD",
        "PD",
        "pDO",
        mat.index,
        ultrametric=slices.ultrametric,
        root_age=tree.root_age,
    )


def cpe_rate(
    tree: TreeArray,
    mat,
    n_slices: int = 100,
    criterion: str = "time",
    percentile: float = 0.95,
    diagnostics: bool = True,
    multistart: int = 8,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
    allow_non_ultrametric: bool = False,
) -> pd.DataFrame:
    """Rate of phylogenetic-endemism accumulation (mirrors ``CpE``).

    See :func:`cpd_rate` for ``bounds``, the bound flags and the ultrametric
    audit column; the returned columns are otherwise identical with
    ``CpE / PE / pEO`` in place of ``CpD / PD / pDO``.
    """
    tree, mat = _validate_common(tree, mat, n_slices, percentile)
    slices = slice_pieces(
        tree, n=n_slices, criterion=criterion, allow_non_ultrametric=allow_non_ultrametric
    )
    profiles = pe_per_slice(tree, mat, slices)
    totals = profiles.sum(axis=1)
    fits = _fit_profiles(profiles, slices.ages, multistart, bounds)
    if not diagnostics:
        for f in fits:
            f.residuals = np.empty(0)
    return _rate_frame(
        fits,
        totals,
        percentile,
        "CpE",
        "PE",
        "pEO",
        mat.index,
        ultrametric=slices.ultrametric,
        root_age=tree.root_age,
    )


def cpb_rate(
    tree: TreeArray,
    mat,
    adj,
    n_slices: int = 100,
    component: str = "sorensen",
    approach: str = "multisite",
    weighted: bool = False,
    criterion: str = "time",
    percentile: float = 0.95,
    multistart: int = 8,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
    multisite_domain: str = "mixed",
    allow_non_ultrametric: bool = False,
) -> pd.DataFrame:
    """Rate of phylogenetic beta-diversity accumulation (mirrors ``CpB``/``CpB_RW``).

    ``weighted=True`` gives the range-weighted variant (``CpB_RW``;
    Laffan et al. 2018).  Returns
    ``site_id, CpB (or CpB_RW), PB (or PB_RW), pBO, converged, r_squared,
    reason, at_lower_bound, at_upper_bound, ultrametric``.
    Neighbourhoods skipped by treesliceR (fewer than 2 sites, degenerate
    composition) are returned as ``converged=False`` rows with an explicit
    ``reason`` instead of silent ``c(NA, NA, NA)``.

    ``multisite_domain`` chooses the summation domain of the shared term in the
    multisite approach: ``"mixed"`` (the default; also spelled ``"reference"``)
    reproduces treesliceR's arithmetic, whose index grows with the neighbourhood
    size (``m / (m + 2)`` on the reference construction, i.e. 0.5 -> 0.91 as m
    goes 2 -> 20), while ``"paired"`` puts every term in the pair-summed domain
    and removes that drift.  The domain changes the reported ``PB`` total for
    every component, and the fitted ``CpB`` for ``nestedness`` only: for
    ``sorensen`` and ``turnover`` the multisite profile is divided by its own
    row sum, which cancels the domain's scalar denominator exactly, so the
    fitted rate is identical either way.  Keeping one domain per analysis is
    therefore a reporting-consistency rule, not a numerical one for those two
    components.  See ``phyloslicer.indices.MULTISITE_DOMAINS``.
    """
    tree, mat = _validate_common(tree, mat, n_slices, percentile)
    adj_arr = np.asarray(adj)
    if adj_arr.shape[0] != mat.shape[0]:
        raise ValueError("adj must have one row per site of mat")
    slices = slice_pieces(
        tree, n=n_slices, criterion=criterion, allow_non_ultrametric=allow_non_ultrametric
    )
    result = pb_per_slice(
        tree,
        mat,
        adj_arr,
        slices,
        component=component,
        approach=approach,
        weighted=weighted,
        multisite_domain=multisite_domain,
    )
    rate_name = "CpB_RW" if weighted else "CpB"
    total_name = "PB_RW" if weighted else "PB"
    origin_name = "pBO"

    # Collect valid beta profiles and fit them in one vectorised batch
    # (same path as cpd_rate / cpe_rate) for consistency and speed.
    valid_idx = [
        i
        for i in range(adj_arr.shape[0])
        if result["status"][i] == "ok"
        and not np.isnan(result["values"][i]).any()
        and np.isfinite(result["total"][i])
    ]
    if valid_idx:
        profiles = np.vstack([result["values"][i] for i in valid_idx])
        fits = _fit_profiles(profiles, slices.ages, multistart, bounds)
    else:
        fits = []
    fit_map = {valid_idx[j]: fits[j] for j in range(len(valid_idx))}

    rows = []
    for i in range(adj_arr.shape[0]):
        series = result["values"][i]
        total = result["total"][i]
        status = result["status"][i]
        if status != "ok" or np.isnan(series).any() or np.isnan(total):
            rows.append(
                {
                    "site_id": mat.index[i],
                    rate_name: np.nan,
                    total_name: total if np.isfinite(total) else np.nan,
                    origin_name: np.nan,
                    "converged": False,
                    "r_squared": np.nan,
                    "reason": status if status != "ok" else "non-finite beta profile",
                    "at_lower_bound": False,
                    "at_upper_bound": False,
                    "ultrametric": slices.ultrametric,
                    "origin_within_tree": True,
                }
            )
            continue
        fit = fit_map[i]
        pXO = origin_time(fit.rate, percentile=percentile)
        origin_within, origin_reason = _origin_within_tree(pXO, tree.root_age, fit.reason)
        rows.append(
            {
                "site_id": mat.index[i],
                rate_name: fit.rate,
                total_name: total,
                origin_name: pXO,
                "converged": fit.converged,
                "r_squared": fit.r_squared,
                "reason": origin_reason,
                "at_lower_bound": fit.at_lower_bound,
                "at_upper_bound": fit.at_upper_bound,
                "ultrametric": slices.ultrametric,
                "origin_within_tree": origin_within,
            }
        )
    return pd.DataFrame(rows)


def rate_profile(
    tree: TreeArray,
    mat,
    adj=None,
    n_slices: int = 100,
    index: str = "PD",
    criterion: str = "time",
    normalize: bool = True,
    multisite_domain: str = "mixed",
    allow_non_ultrametric: bool = False,
) -> dict[str, np.ndarray]:
    """Per-slice index profile without fitting.

    Parameters
    ----------
    normalize : bool
        ``True`` (default) returns each row divided by its own total, i.e. a
        *relative* profile whose finite rows sum to 1.  ``False`` returns the
        raw per-slice values and reproduces ``treesliceR::r_phylo``, which for
        ``index="PD"`` returns ``sum(edge.length[positions])`` (absolute PD per
        slice, row sums equal to the site's total PD) and for ``"PE"`` the
        range-weighted absolute values.  For the beta indices ``normalize=False``
        returns the curve ``CpB`` feeds into its fit, i.e. ``r_pbd / pbd``; the
        undivided ``r_phylo(index = "PB")`` curve is the ``"series"`` key of
        :func:`phyloslicer.indices.pb_per_slice`, because the R implementation
        relativises PB against the whole-tree index while *naming* it "relative"
        (documented in validation/NOTES.md).
    index : {"PD", "PE", "PB", "PB_RW"}
    multisite_domain : {"mixed", "reference", "paired"}
        Only used by the beta indices; ``"reference"`` is a synonym of the
        default ``"mixed"``.  See
        :data:`phyloslicer.indices.MULTISITE_DOMAINS`.

    Returns ``{"values": (n_sites, k), "ages": (k,), "status": (n_sites,),
    "ultrametric": bool}``, plus ``"normalisation"`` for the beta indices.
    """
    if index not in ("PD", "PE", "PB", "PB_RW"):
        raise ValueError("index must be one of PD, PE, PB, PB_RW")
    if not isinstance(normalize, bool):
        raise TypeError("normalize must be a bool")
    tree, mat = align_tree_matrix(tree, as_matrix(mat))
    slices = slice_pieces(
        tree, n=n_slices, criterion=criterion, allow_non_ultrametric=allow_non_ultrametric
    )
    if index in ("PB", "PB_RW"):
        if adj is None:
            raise ValueError(f"index={index} requires an adjacency matrix")
        if normalize:
            result = pb_per_slice(
                tree,
                mat,
                np.asarray(adj),
                slices,
                component="sorensen",
                approach="multisite",
                weighted=(index == "PB_RW"),
                multisite_domain=multisite_domain,
            )
        else:
            # treesliceR hands `r_pbd` (or the un-normalised turnover curve) to
            # nls; do not re-scale it here
            from ..indices.api import r_compatible_profile

            result = r_compatible_profile(
                tree,
                mat,
                np.asarray(adj),
                slices,
                component="sorensen",
                approach="multisite",
                weighted=(index == "PB_RW"),
                multisite_domain=multisite_domain,
            )
        return {
            "values": result["values"],
            "ages": slices.ages,
            "status": result["status"],
            "normalisation": result["normalisation"],
            "ultrametric": slices.ultrametric,
        }
    prof = pd_per_slice(tree, mat, slices) if index == "PD" else pe_per_slice(tree, mat, slices)
    total = prof.sum(axis=1, keepdims=True)
    status = np.array(["ok"] * prof.shape[0], dtype=object)
    status[(total[:, 0] == 0) | ~np.isfinite(total[:, 0])] = "empty_site"
    if not normalize:
        return {
            "values": prof,
            "ages": slices.ages,
            "status": status,
            "ultrametric": slices.ultrametric,
        }
    with np.errstate(divide="ignore", invalid="ignore"):
        return {
            "values": prof / total,
            "ages": slices.ages,
            "status": status,
            "ultrametric": slices.ultrametric,
        }


def sensitivity(
    tree: TreeArray,
    mat,
    adj=None,
    slice_grid: np.ndarray | list[int] = (10, 25, 50, 100, 200, 400),
    rate: str = "cpd",
    sample_sites: int = 0,
    component: str = "sorensen",
    approach: str = "multisite",
    weighted: bool = False,
    seed: int | None = None,
    bounds: tuple[float, float] | str = DEFAULT_RATE_BOUNDS,
    multisite_domain: str = "mixed",
    allow_non_ultrametric: bool = False,
) -> pd.DataFrame:
    """Sensitivity of a cumulative rate to the number of slices.

    Mirrors treesliceR's ``CpR_sensitivity`` but ``sample_sites=0`` uses every
    site (R requires a positive sample).  Also reports a heuristic suggested
    slice count per site: the smallest grid entry after which the rate stops
    changing monotonically by more than 5%.  The suggestion is a numerical
    stability heuristic, not a convergence diagnostic - ``converged`` and the
    bound flags are reported alongside so a clipped rate cannot be mistaken for
    a stable one.
    """
    if rate not in ("cpd", "cpe", "cpb", "cpb_rw"):
        raise ValueError("rate must be one of cpd, cpe, cpb, cpb_rw")
    if rate in ("cpb", "cpb_rw") and adj is None:
        raise ValueError(f"rate={rate} requires an adjacency matrix")
    if sample_sites < 0:
        raise ValueError("sample_sites must be non-negative (0 = use every site)")
    tree, mat = align_tree_matrix(tree, as_matrix(mat))
    tree.check_ultrametric(allow_non_ultrametric=allow_non_ultrametric)
    rng = np.random.default_rng(seed)
    n_sites = mat.shape[0]
    if sample_sites and sample_sites > n_sites:
        raise ValueError("the number of site samples is bigger than the sites available")
    if sample_sites:
        chosen = np.sort(rng.choice(n_sites, size=sample_sites, replace=False))
    else:
        chosen = np.arange(n_sites)

    records = []
    for n in slice_grid:
        n = int(n)
        slices = slice_pieces(
            tree, n=n, criterion="time", allow_non_ultrametric=allow_non_ultrametric
        )
        if rate == "cpd":
            prof = pd_per_slice(tree, mat, slices)[chosen]
        elif rate == "cpe":
            prof = pe_per_slice(tree, mat, slices)[chosen]
        else:
            result = pb_per_slice(
                tree,
                mat,
                np.asarray(adj),
                slices,
                component=component,
                approach=approach,
                weighted=(rate == "cpb_rw"),
                multisite_domain=multisite_domain,
            )
            prof = result["values"][chosen]
        ages = slices.ages
        for local_i, site_i in enumerate(chosen):
            fit = fit_rate(prof[local_i], ages, bounds=bounds)
            records.append(
                {
                    "site_id": mat.index[site_i],
                    "n_slices": n,
                    "rate": fit.rate,
                    "converged": fit.converged,
                    "at_lower_bound": fit.at_lower_bound,
                    "at_upper_bound": fit.at_upper_bound,
                }
            )
    df = pd.DataFrame.from_records(records)

    # heuristic: smallest slice count where the rate is within 5% of the
    # next (larger) grid estimate
    suggested = []
    for site, grp in df.groupby("site_id", sort=False):
        vals = grp.sort_values("n_slices")["rate"].to_numpy()
        pick = np.nan
        for j in range(len(vals) - 1):
            if np.isfinite(vals[j]) and np.isfinite(vals[j + 1]):
                denom = max(abs(vals[j + 1]), 1e-12)
                if abs(vals[j] - vals[j + 1]) / denom <= 0.05:
                    pick = grp.sort_values("n_slices")["n_slices"].to_numpy()[j]
                    break
        suggested.append(
            {"site_id": site, "suggested_slices": pick, "ultrametric": slices.ultrametric}
        )
    sug = pd.DataFrame(suggested)
    return df.merge(sug, on="site_id", how="left")
