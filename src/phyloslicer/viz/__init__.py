"""Publication-oriented static visualisations (matplotlib backend)."""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["plot_rate_line", "plot_sensitivity", "plot_posterior"]


def plot_rate_line(
    profile, fit=None, ages=None, percentile: float = 0.95, ax=None, direction: str = "from_root"
):
    """Cumulative curve + fitted exponential + origin-time marker.

    ``profile`` is a per-slice vector (root -> tip); ``fit`` is a RateFit;
    ``ages`` the slice boundary ages (Ma).  Mirrors ``CpR_graph`` line mode.

    ``direction`` chooses the x-axis convention.  The default ``"from_root"``
    increases from the root towards the present, which is the convention the
    worked examples in the user guide use and the one a reader gets by calling
    this function with no extra arguments.  ``"from_present"`` reproduces the
    treesliceR view (x is age in Ma before present), which is what a reader
    coming from the R package expects; the two views are mirror images about the
    root age, so mixing them changes the direction the curve runs.
    """
    import matplotlib.pyplot as plt

    if not 0 < percentile < 1:
        # origin_time() rejects this, and drawing it instead parked an axvline
        # at NaN or -inf with a legend reading "150% origin = nan Ma" - which a
        # user who wrote percentile=95 would only notice in the figure
        raise ValueError(f"percentile must lie strictly between 0 and 1, got {percentile!r}")
    profile = np.asarray(profile, dtype=float)
    k = profile.size
    if ages is None:
        import warnings

        warnings.warn(
            "ages not provided; using placeholder spacing — the x-axis "
            "will be incorrect. Pass ages=slices.ages for correct scaling.",
            stacklevel=2,
        )
        ages = np.linspace(k, 1, k)  # placeholder spacing when unknown
    ages = np.asarray(ages, dtype=float)
    if direction not in ("from_root", "from_present"):
        raise ValueError('direction must be "from_root" or "from_present"')
    # x in the chosen convention; the model is always evaluated in age coordinates
    x_obs = ages.max() - ages if direction == "from_root" else ages
    cum = np.cumsum(profile) / max(profile.sum(), 1e-300)

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(x_obs, cum, s=18, color="#1b6b4a", label="observed cumulative", zorder=3)
    if fit is not None and np.isfinite(fit.rate):
        # the fit is drawn only across the observed interval, exactly as in the
        # worked examples: the model's asymptote at the present is beyond the
        # last measured slice boundary, so extending the curve there would plot
        # a region with no data behind it
        grid = np.linspace(ages.min(), ages.max(), 200)
        x_fit = ages.max() - grid if direction == "from_root" else grid
        ax.plot(
            x_fit,
            np.exp(-fit.rate * grid),
            color="#d99b2b",
            lw=2,
            label=f"fit: r = {fit.rate:.3g} (R² = {fit.r_squared:.3f})",
        )
        x0 = -np.log(1 - percentile) / fit.rate
        # the marker is an age in Ma (the label says so), but it has to be drawn
        # in the same coordinates as the points and the curve: mirroring the
        # latter without mirroring this one parked the origin line past the
        # present end of a from-root axis.
        ax.axvline(
            ages.max() - x0 if direction == "from_root" else x0,
            color="0.3",
            ls="--",
            lw=1,
            label=f"{int(percentile * 100)}% origin = {x0:.3g} Ma",
        )
    ax.set_xlabel("time from root (Myr)" if direction == "from_root" else "age (Ma before present)")
    ax.set_ylabel("cumulative proportion")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False, fontsize=8)
    return ax


def plot_sensitivity(sens_df: pd.DataFrame, ax=None):
    """Rate vs number-of-slices sensitivity plot (mirrors ``CpR_sensitivity_plot``)."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    for site, grp in sens_df.groupby("site_id"):
        grp = grp.sort_values("n_slices")
        ax.plot(grp["n_slices"], grp["rate"], marker="o", ms=3, lw=1, alpha=0.8, label=str(site))
    ax.set_xlabel("number of slices")
    ax.set_ylabel("rate")
    ax.set_xscale("log")
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) <= 8:
        ax.legend(frameon=False, fontsize=8, title="site")
    else:
        ax.legend().remove()
    return ax


def plot_posterior(
    posterior_df: pd.DataFrame,
    rate: str = "rate_median",
    kind: str = "ecdf",
    ax=None,
    hdi: float | None = None,
    hdi_label: str | None = None,
):
    """Tree-set rate distribution (ECDF or histogram) with HDI shading.

    ``hdi`` is the credible level the shaded interval was computed at.  When
    omitted it is read from ``posterior_df.attrs["hdi_level"]``, which
    :func:`phyloslicer.posterior_rate` / :func:`~phyloslicer.posterior_origin`
    write.  Passing neither is an error rather than a silent guess: the previous
    version hard-coded a fallback of ``0.95``, so a user who had called
    ``posterior_rate(..., hdi=0.68)`` got a panel annotated "95% HDI" around a
    68% interval - a mislabelled figure.

    The shading needs a single interval, so it is drawn only for one-row
    (single-site) summaries; for a multi-site frame the pooled distribution is
    shown without shading and :func:`plot_posterior` notes it in the axes title.
    """
    import matplotlib.pyplot as plt

    x = posterior_df[rate].to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        raise ValueError(f"no finite values in column {rate!r} to plot")
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 4))
    if kind == "ecdf":
        xs = np.sort(x)
        ys = np.arange(1, xs.size + 1) / xs.size
        ax.step(xs, ys, where="post", color="#1b6b4a")
    elif kind == "hist":
        ax.hist(x, bins=30, color="#1b6b4a", alpha=0.85)
    else:
        raise ValueError("kind must be 'ecdf' or 'hist'")
    if hdi is None:
        hdi = posterior_df.attrs.get("hdi_level")
    if {"hdi_low", "hdi_high"}.issubset(posterior_df.columns) and len(posterior_df) == 1:
        lo, hi = posterior_df["hdi_low"].iloc[0], posterior_df["hdi_high"].iloc[0]
        if hdi is None:
            raise ValueError(
                "cannot label the HDI: posterior_df carries no attrs['hdi_level'] "
                "(it was not produced by posterior_rate/posterior_origin in this "
                "session) - pass hdi=0.95 explicitly"
            )
        if not 0 < hdi <= 1:
            raise ValueError("hdi must lie in (0, 1]")
        label = hdi_label or f"{hdi * 100:g}% HDI"
        ax.axvspan(lo, hi, color="#d99b2b", alpha=0.3, label=label)
        ax.legend(frameon=False, fontsize=8)
    elif hdi is not None and len(posterior_df) > 1:
        ax.set_title(
            f"{rate} over {len(posterior_df)} sites (no single {hdi * 100:g}% HDI)", fontsize=9
        )
    ax.set_xlabel(rate)
    ax.set_ylabel("cumulative proportion" if kind == "ecdf" else "count")
    return ax
