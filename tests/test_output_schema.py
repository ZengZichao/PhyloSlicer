"""The output schema of the frame-returning public APIs is pinned.

The 0.1.0 release audit found that the rate fitters and :func:`posterior_rate`
documented fewer columns than they returned: ``origin_within_tree`` and
``n_ultrametric_draws`` had been added to the frames without the prose being
carried along, and nothing failed.  Prose cannot fail a build, so this module
pins the schema itself and then checks that each docstring still names every
column of the frame it describes - either side drifting alone is a failure.

The workload is the bundled simulation at a small size (24 sites / 60 species /
4 draws, 10 slices): the returned columns do not depend on how many sites or
slices a run uses, so a fast workload pins the same schema as the 208-site
default.
"""

from __future__ import annotations

import pandas as pd

from phyloslicer.datasets import sim_passerines
from phyloslicer.rates.api import cpb_rate, cpd_rate, cpe_rate
from phyloslicer.uncertainty import bootstrap_sites, posterior_origin, posterior_rate

#: columns every cumulative-rate frame carries, whatever index it fitted
_AUDIT_COLUMNS = {
    "site_id",
    "converged",
    "r_squared",
    "reason",
    "at_lower_bound",
    "at_upper_bound",
    "ultrametric",
    "origin_within_tree",
}


def _rate_columns(rate: str, total: str, origin: str) -> set[str]:
    """The eleven columns of a cumulative-rate frame: audit block + three."""
    return _AUDIT_COLUMNS | {rate, total, origin}


def test_rate_and_uncertainty_output_columns_are_pinned():
    """Each frame-returning API returns exactly its documented column set."""
    data = sim_passerines(n_sites=24, n_species=60, n_trees=4, seed=11)
    trees, mat, adj = data["trees"], data["mat"], data["adj"]
    tree = trees[0]

    n_slices = 10
    cases = {
        "cpd_rate": (
            cpd_rate,
            (tree, mat),
            {"n_slices": n_slices},
            _rate_columns("CpD", "PD", "pDO"),
        ),
        "cpe_rate": (
            cpe_rate,
            (tree, mat),
            {"n_slices": n_slices},
            _rate_columns("CpE", "PE", "pEO"),
        ),
        "cpb_rate": (
            cpb_rate,
            (tree, mat, adj),
            {"n_slices": n_slices},
            _rate_columns("CpB", "PB", "pBO"),
        ),
        # the range-weighted variant renames the two value columns and only them
        "cpb_rate(weighted=True)": (
            cpb_rate,
            (tree, mat, adj),
            {"n_slices": n_slices, "weighted": True},
            _rate_columns("CpB_RW", "PB_RW", "pBO"),
        ),
        "posterior_rate": (
            posterior_rate,
            (trees, mat),
            {"rate": "cpd", "n_slices": 10},
            {
                "site_id",
                "rate_mean",
                "rate_median",
                "hdi_low",
                "hdi_high",
                "ess_proxy",
                "ess_proxy_caveat",
                "n_trees",
                "n_converged",
                "n_ultrametric_draws",
                "origin_within_tree",
            },
        ),
        "posterior_origin": (
            posterior_origin,
            (trees, mat),
            {"rate": "cpd", "n_slices": 10},
            {
                "site_id",
                "origin_median",
                "hdi_low",
                "hdi_high",
                "n_ultrametric_draws",
                "origin_within_tree",
            },
        ),
        "bootstrap_sites": (
            bootstrap_sites,
            (tree, mat),
            {"n_boot": 20, "n_slices": 10, "seed": 3},
            {"replicate", "rate"},
        ),
    }

    for name, (func, args, kwargs, expected) in cases.items():
        frame = func(*args, **kwargs)
        assert isinstance(frame, pd.DataFrame), f"{name} no longer returns a DataFrame"
        returned = set(frame.columns)
        assert returned == expected, (
            f"{name} output schema changed: gained {sorted(returned - expected)}, "
            f"lost {sorted(expected - returned)}"
        )
        # and the published prose still lists every one of those columns
        documented = func.__doc__ or ""
        undocumented = sorted(column for column in expected if column not in documented)
        assert not undocumented, f"{name} returns columns its docstring never names: {undocumented}"
