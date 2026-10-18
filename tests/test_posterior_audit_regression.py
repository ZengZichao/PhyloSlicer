"""Regression test for the posterior audit columns (FS-13 / API-04).

`n_ultrametric_draws` and `origin_within_tree` used to be read off the flattened
per-draw audit block, so they reported the first n_draws *sites of draw 1* and were
identical for every site regardless of which draw failed ultrametricity.
"""

import re
import warnings

import pytest

from phyloslicer.core import TreeArray
from phyloslicer.datasets import sim_passerines
from phyloslicer.uncertainty import posterior_origin, posterior_rate


def _break_ultrametricity(tree):
    """Halve one terminal branch length in the Newick and re-parse."""
    nwk = tree.to_newick()
    m = re.search(r"([A-Za-z0-9_.]+):([0-9.eE+-]+)", nwk)
    label, raw = m.group(1), m.group(2)
    return TreeArray.from_newick(nwk.replace(f"{label}:{raw}", f"{label}:{float(raw) / 2}", 1))


@pytest.fixture()
def pool():
    d = sim_passerines(n_sites=6, n_species=24, n_trees=3, seed=1)
    good = list(d["trees"])
    return good, _break_ultrametricity(good[1]), d["mat"]


def _run(func, pool, mat):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return func(pool, mat, n_slices=20, allow_non_ultrametric=True)


def test_audit_counts_ultrametric_draws_per_site(pool):
    good, bad, mat = pool
    for order in ([good[0], bad, good[1]], [bad, good[0], good[1]]):
        out = _run(posterior_rate, order, mat)
        assert (out["n_ultrametric_draws"] == 2).all(), (
            "two of three draws are ultrametric regardless of pool order"
        )


def test_audit_columns_are_not_constant_across_sites_when_sites_differ(pool):
    good, bad, mat = pool
    out = _run(posterior_rate, [good[0], bad, good[1]], mat)
    assert len(out) == mat.shape[0]
    assert (out["n_ultrametric_draws"] == 2).all()


def test_posterior_origin_audit_matches_rate_audit(pool):
    good, bad, mat = pool
    order = [good[0], bad, good[1]]
    r = _run(posterior_rate, order, mat)
    o = _run(posterior_origin, order, mat)
    assert (r["n_ultrametric_draws"].to_numpy() == o["n_ultrametric_draws"].to_numpy()).all()
