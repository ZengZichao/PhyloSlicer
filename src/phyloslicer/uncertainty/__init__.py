"""Posterior-tree-set uncertainty quantification (novel vs treesliceR).

Trees are consumed through a streaming generator so a posterior pool larger
than memory can be processed; per-tree work reuses the vectorised kernels.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Iterable, Iterator

import numpy as np
import pandas as pd

from ..core.tree import TreeArray
from ..io.matrices import as_matrix
from ..io.newick import (
    _is_existing_dir,
    _is_existing_file,
    _looks_like_newick,
    read_newick,
    read_nexus,
    split_newick_records,
)
from ..rates.api import cpd_rate, cpe_rate

__all__ = [
    "posterior_rate",
    "posterior_origin",
    "bootstrap_sites",
    "iter_trees",
    "ESS_PROXY_CAVEAT",
]


def _trees_from_text(text: str, origin: str) -> Iterator[TreeArray]:
    """Yield every tree in a Newick text, **warning** about records that fail.

    A posterior pool is only meaningful if an unreadable record is reported:
    silently dropping draws shrinks ``n_trees`` with no trace at all, which is
    precisely the failure this module promises not to make.
    """
    try:
        records = split_newick_records(text)
    except ValueError as exc:
        warnings.warn(f"{origin}: {exc}; the whole source was skipped", stacklevel=3)
        return
    for j, record in enumerate(records):
        try:
            yield read_newick(record)
        except Exception as exc:  # a corrupt record must not sink the pool
            warnings.warn(
                f"{origin}: skipping record {j + 1} of {len(records)}: {exc}", stacklevel=3
            )


def iter_trees(source) -> Iterator[TreeArray]:
    """Yield trees from a list, a Newick string, a multi-tree Newick file, or a directory.

    A ``str`` is resolved the same way :func:`phyloslicer.io.read_newick` resolves
    it: an existing file or directory first, then a string that carries Newick
    structure, and only then as a missing path.  Passing a multi-tree Newick
    *string* therefore works instead of reaching :class:`pathlib` and failing
    with the operating system's ``File name too long``.
    """
    if isinstance(source, TreeArray):
        yield source
        return
    if isinstance(source, (str, Path)):
        as_str = str(source)
        p = Path(as_str)
        if _is_existing_dir(as_str):
            files = sorted(
                q
                for q in p.iterdir()
                if q.suffix.lower() in (".tre", ".tree", ".nwk", ".newick", ".nex", ".t")
            )
            if not files:
                raise ValueError(f"no tree files found in directory {p}")
            for f in files:
                if f.suffix.lower() == ".nex":
                    yield from read_nexus(f)
                else:
                    yield from _trees_from_text(f.read_text(), str(f.name))
            return
        if _is_existing_file(as_str):
            if p.suffix.lower() == ".nex":
                yield from read_nexus(p)
            else:
                yield from _trees_from_text(p.read_text(), str(p.name))
            return
        if _looks_like_newick(as_str.strip()):
            yield from _trees_from_text(as_str, "<Newick string>")
            return
        raise FileNotFoundError(
            f"{as_str[:200]!r} is not an existing tree file or directory, nor a "
            "Newick string (no '(', ')', ':' or ';' found); pass a path to a "
            "multi-tree file, a directory of trees, or Newick text"
        )
    if isinstance(source, Iterable):
        for t in source:
            if not isinstance(t, TreeArray):
                raise TypeError("tree collections must contain TreeArray objects")
            yield t
        return
    raise TypeError("source must be a TreeArray, path, directory or iterable of trees")


def _hdi(samples: np.ndarray, cred: float) -> tuple[float, float]:
    """Nearest-neighbour empirical highest-density interval."""
    samples = np.sort(samples[np.isfinite(samples)])
    n = samples.size
    if n == 0:
        return np.nan, np.nan
    k = max(1, int(np.ceil(cred * n)))
    if k >= n:
        return float(samples[0]), float(samples[-1])
    widths = samples[k - 1 :] - samples[: n - k + 1]
    j = int(np.argmin(widths))
    return float(samples[j]), float(samples[j + k - 1])


#: Short text repeated in every `ess_proxy_caveat` output column, so the
#: heuristic cannot be mistaken for a convergence diagnostic once the frame has
#: been written to CSV (DataFrame attrs do not survive a round-trip).
ESS_PROXY_CAVEAT = (
    "heuristic only: assumes on-disk draw order == chain order and stops the "
    "autocorrelation sum at the first rho < 0.05, which biases ESS upward; "
    "not an MCMC convergence diagnostic (use Tracer/coda/arviz)"
)


def _ess_proxy(samples: np.ndarray) -> float:
    """Crude effective-sample-size proxy from lag-1..lag-10 autocorrelation.

    Two limitations are deliberate and are surfaced to the caller through
    :data:`ESS_PROXY_CAVEAT` (written into the result frame as the
    ``ess_proxy_caveat`` column and as ``attrs["ess_proxy"]``):

    1. the initialisation ``rho_sum`` loop stops at the *first* lag with
       ``rho < 0.05``, so later, possibly positive, autocorrelations are ignored
       and the result is an **over**-estimate of ESS (an under-estimate of
       chain correlation);
    2. the samples are assumed to be in chain order simply because that is the
       order they were read from disk.  A thinned, shuffled or
       per-file-partitioned posterior pool breaks that assumption and the
       number becomes meaningless.

    With < 10 draws the raw count is returned, which is the least-bad answer
    for a quantity that cannot be estimated at all from so few points.
    """
    x = np.asarray(samples, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 10:
        return float(n)
    x = x - x.mean()
    denom = float(x @ x)
    if denom <= 0:
        return float(n)
    rho_sum = 0.0
    for lag in range(1, 11):
        if lag >= n:
            break
        rho = float(x[:-lag] @ x[lag:]) / denom
        if rho < 0.05:
            break
        rho_sum += rho
    return float(n / (1.0 + 2.0 * rho_sum))


def _validate_hdi(hdi: float) -> None:
    if not 0 < hdi <= 1:
        raise ValueError("hdi must lie in (0, 1]")


def _audit_row(res):
    """Per-draw audit flags, as ``[ultrametric, origin_within_tree]``.

    ``posterior_rate`` and ``posterior_origin`` used to drop the audit columns
    the single-tree fitters return, so a pooled estimate gave no way to see that
    some draws were non-ultrametric or that every draw put the origin outside
    the tree.  A missing column counts as clean, which keeps third-party
    ``runner`` callables working.
    """

    def flag(name):
        if name not in res.columns:
            return np.ones(len(res), dtype=bool)
        return res[name].to_numpy(dtype=bool)

    return np.column_stack([flag("ultrametric"), flag("origin_within_tree")])


def posterior_rate(
    trees,
    mat,
    rate: str = "cpd",
    n_slices: int = 100,
    hdi: float = 0.95,
    **rate_kwargs,
) -> pd.DataFrame:
    """Per-site rate summary across a posterior tree set.

    Returns ``site_id, rate_mean, rate_median, hdi_low, hdi_high, ess_proxy,
    ess_proxy_caveat, n_trees, n_converged``.  ``hdi_low``/``hdi_high`` are the
    nearest-neighbour highest-density interval covering ``hdi`` of the draws,
    and the applied level is recorded in ``df.attrs["hdi_level"]`` so that
    :func:`phyloslicer.viz.plot_posterior` labels the figure with the level that
    was actually used (a DataFrame read back from CSV loses attrs, so pass
    ``hdi=`` there).  ``ess_proxy`` is a crude heuristic, not a convergence
    diagnostic - see :data:`ESS_PROXY_CAVEAT`.
    """
    if rate not in ("cpd", "cpe"):
        raise ValueError("rate must be 'cpd' or 'cpe'")
    _validate_hdi(hdi)
    mat = as_matrix(mat)
    runner = cpd_rate if rate == "cpd" else cpe_rate
    rate_col = "CpD" if rate == "cpd" else "CpE"

    per_tree: list[np.ndarray] = []
    audit: list[np.ndarray] = []
    n_trees = 0
    for tree in iter_trees(trees):
        n_trees += 1
        try:
            res = runner(tree, mat, n_slices=n_slices, **rate_kwargs)
        except Exception as exc:  # a bad posterior draw must not sink the pool
            warnings.warn(f"skipping tree {n_trees}: {exc}", stacklevel=2)
            continue
        per_tree.append(res[rate_col].to_numpy(dtype=float))
        audit.append(_audit_row(res))
    if not per_tree:
        raise RuntimeError("no posterior tree produced rates")
    if len(per_tree) < n_trees:
        warnings.warn(
            f"{n_trees - len(per_tree)} of {n_trees} posterior trees were skipped; "
            f"every summary below rests on {len(per_tree)} draws",
            stacklevel=2,
        )
    rates = np.vstack(per_tree)  # (n_trees, n_sites)
    audit_flags = np.vstack(audit)  # (n_trees, 2): [ultrametric, origin_within_tree]

    rows = []
    for j, site in enumerate(mat.index):
        col = rates[:, j]
        lo, hi = _hdi(col, hdi)
        rows.append(
            {
                "site_id": site,
                "rate_mean": float(np.nanmean(col)),
                "rate_median": float(np.nanmedian(col)),
                "hdi_low": lo,
                "hdi_high": hi,
                "ess_proxy": _ess_proxy(col),
                "ess_proxy_caveat": ESS_PROXY_CAVEAT,
                "n_trees": int(col.size),
                "n_converged": int(np.isfinite(col).sum()),
                "n_ultrametric_draws": int(audit_flags[:, 0][: col.size].sum()),
                "origin_within_tree": bool(audit_flags[:, 1][: col.size].all()),
            }
        )
    out = pd.DataFrame(rows)
    out.attrs["hdi_level"] = float(hdi)
    out.attrs["ess_proxy"] = ESS_PROXY_CAVEAT
    return out


def posterior_origin(
    trees,
    mat,
    rate: str = "cpd",
    n_slices: int = 100,
    percentile: float = 0.95,
    hdi: float = 0.95,
    **rate_kwargs,
) -> pd.DataFrame:
    """Credible interval of the diversity-origin time pXO across posterior trees.

    ``hdi_low``/``hdi_high`` are computed at the requested ``hdi`` level and the
    level is recorded in ``df.attrs["hdi_level"]``.
    """

    if rate not in ("cpd", "cpe"):
        raise ValueError("rate must be 'cpd' or 'cpe'")
    _validate_hdi(hdi)
    mat = as_matrix(mat)
    runner = cpd_rate if rate == "cpd" else cpe_rate
    origin_col = "pDO" if rate == "cpd" else "pEO"
    per_tree: list[np.ndarray] = []
    audit: list[np.ndarray] = []
    n_trees = 0
    for tree in iter_trees(trees):
        n_trees += 1
        try:
            res = runner(tree, mat, n_slices=n_slices, percentile=percentile, **rate_kwargs)
        except Exception as exc:  # a bad posterior draw must not sink the pool
            warnings.warn(f"skipping tree {n_trees}: {exc}", stacklevel=2)
            continue
        per_tree.append(res[origin_col].to_numpy(dtype=float))
        audit.append(_audit_row(res))
    if not per_tree:
        raise RuntimeError("no posterior tree produced origin times")
    if len(per_tree) < n_trees:
        warnings.warn(
            f"{n_trees - len(per_tree)} of {n_trees} posterior trees were skipped; "
            f"every summary below rests on {len(per_tree)} draws",
            stacklevel=2,
        )
    origins = np.vstack(per_tree)
    audit_flags = np.vstack(audit)
    rows = []
    for j, site in enumerate(mat.index):
        col = origins[:, j]
        lo, hi = _hdi(col, hdi)
        rows.append(
            {
                "site_id": site,
                "origin_median": float(np.nanmedian(col)),
                "hdi_low": lo,
                "hdi_high": hi,
                "n_ultrametric_draws": int(audit_flags[:, 0][: col.size].sum()),
                "origin_within_tree": bool(audit_flags[:, 1][: col.size].all()),
            }
        )
    out = pd.DataFrame(rows)
    out.attrs["hdi_level"] = float(hdi)
    return out


def bootstrap_sites(
    tree: TreeArray,
    mat,
    n_boot: int = 1000,
    rate: str = "cpd",
    n_slices: int = 100,
    seed: int | None = None,
    **rate_kwargs,
) -> pd.DataFrame:
    """Stability of per-site rates under site resampling (with replacement).

    The rate model is refit on the *pooled* assemblage of each bootstrap
    replicate; returns per-replicate pooled rates.
    """
    from ..indices.api import pd_per_slice, pe_per_slice
    from ..rates.fitting import fit_rate
    from ..slicing.api import slice_pieces

    if rate not in ("cpd", "cpe"):
        raise ValueError("rate must be 'cpd' or 'cpe'")
    mat = as_matrix(mat)
    # align first (prune tree tips absent from the matrix, drop matrix
    # species absent from the tree) exactly like cpd_rate/cpe_rate so the
    # bootstrap operates on the same aligned tree as the point estimates
    from ..io.matrices import align_tree_matrix

    tree, mat = align_tree_matrix(tree, mat)
    rng = np.random.default_rng(seed)
    slices = slice_pieces(tree, n=n_slices)
    prof = (pd_per_slice if rate == "cpd" else pe_per_slice)(tree, mat, slices)
    ages = slices.ages
    rates = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, mat.shape[0], size=mat.shape[0])
        pooled = prof[idx].sum(axis=0)
        fit = fit_rate(pooled, ages)
        if not fit.converged:
            warnings.warn(
                f"bootstrap replicate {b} did not converge: {fit.reason}",
                stacklevel=2,
            )
        rates[b] = fit.rate
    return pd.DataFrame({"replicate": np.arange(n_boot), "rate": rates})
