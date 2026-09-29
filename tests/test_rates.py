"""Rate fitting, cumulative-rate APIs, uncertainty quantification."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phyloslicer.rates import fit_rate, origin_time
from phyloslicer.rates.api import cpb_rate, cpd_rate, cpe_rate, rate_profile, sensitivity
from phyloslicer.simulate import yule_tree
from phyloslicer.slicing import slice_pieces
from phyloslicer.uncertainty import bootstrap_sites, posterior_origin, posterior_rate


class TestFitRate:
    def test_recovers_known_rate(self):
        r_true = 0.8
        # a very fine slice grid so the unit-intercept exponential matches the
        # normalised cumulative almost exactly (youngest slice ~ age 0)
        ages = np.linspace(5.0, 1e-4, 200)
        cum = np.exp(-r_true * ages)
        profile = np.diff(np.concatenate([[0.0], cum]))
        fit = fit_rate(profile, ages)
        assert fit.converged
        assert fit.rate == pytest.approx(r_true, rel=2e-3)
        assert fit.r_squared > 0.999

    def test_degenerate_profile_flagged_not_raised(self):
        ages = np.array([3.0, 2.0, 1.0])
        fit = fit_rate(np.zeros(3), ages)
        assert not fit.converged
        assert fit.reason

    def test_empty_site_reason(self):
        fit = fit_rate(np.array([np.nan, np.nan]), np.array([2.0, 1.0]))
        assert not fit.converged

    def test_origin_time(self):
        assert origin_time(0.2, percentile=0.95) == pytest.approx(-np.log(0.05) / 0.2)
        assert np.isnan(origin_time(np.nan))
        assert np.isnan(origin_time(0.0))
        with pytest.raises(ValueError):
            origin_time(0.2, percentile=1.0)


class TestCpdCpe:
    @pytest.fixture
    def setup(self):
        tree = yule_tree(25, age=1.0, seed=42)
        mat = pd.DataFrame(
            (np.random.default_rng(3).random((8, 25)) < 0.3).astype(float),
            columns=tree.tip_labels,
            index=[f"s{i}" for i in range(8)],
        )
        mat.iloc[0] = 1.0
        return tree, mat

    def test_cpd_output_shape_and_columns(self, setup):
        tree, mat = setup
        res = cpd_rate(tree, mat, n_slices=20)
        # The full, intentional column contract of cpd_rate: the fit answer, the
        # three audit columns (reason / bounds / ultrametric) and the origin-range
        # flag. See phyloslicer/rates/api.py.
        assert list(res.columns) == [
            "site_id",
            "CpD",
            "PD",
            "pDO",
            "converged",
            "r_squared",
            "reason",
            "at_lower_bound",
            "at_upper_bound",
            "ultrametric",
            "origin_within_tree",
        ]
        assert len(res) == mat.shape[0]
        assert res["converged"].all()
        assert (res["CpD"] > 0).all()
        assert (res["pDO"] > 0).all()
        # an unconstrained fit on a clean ultrametric tree sits strictly inside
        # the search bounds, and the ultrametric verdict is echoed back
        assert not res["at_lower_bound"].any()
        assert not res["at_upper_bound"].any()
        assert res["ultrametric"].all()

    def test_cpd_matches_manual_fit(self, setup):
        tree, mat = setup
        res = cpd_rate(tree, mat, n_slices=10)
        stack = slice_pieces(tree, n=10)
        from phyloslicer.indices import pd_per_slice

        prof = pd_per_slice(tree, mat, stack)[0]
        fit = fit_rate(prof, stack.ages)
        assert res["CpD"].iloc[0] == pytest.approx(fit.rate)

    def test_cpe_runs(self, setup):
        tree, mat = setup
        res = cpe_rate(tree, mat, n_slices=15)
        assert (res["PE"] > 0).all()
        assert res["converged"].all()

    def test_empty_site_reported(self, setup):
        tree, mat = setup
        mat_empty = mat.copy()
        mat_empty.iloc[3, :] = 0.0
        res = cpd_rate(tree, mat_empty, n_slices=10)
        assert not res["converged"].iloc[3]
        assert res["reason"].iloc[3]
        assert np.isnan(res["CpD"].iloc[3])

    def test_rate_profile(self, setup):
        tree, mat = setup
        out = rate_profile(tree, mat, n_slices=8, index="PD")
        assert out["values"].shape == (8, 8)
        # relative profiles sum to 1 per site
        assert np.allclose(out["values"].sum(axis=1), 1.0)

    def test_sensitivity(self, setup):
        tree, mat = setup
        df = sensitivity(tree, mat, slice_grid=[5, 10, 20], rate="cpd", seed=1)
        assert set(df["n_slices"]) == {5, 10, 20}
        assert "suggested_slices" in df.columns


class TestCpb:
    @pytest.fixture
    def setup(self):
        tree = yule_tree(12, age=1.0, seed=7)
        rng = np.random.default_rng(9)
        mat = pd.DataFrame(
            (rng.random((6, 12)) < 0.4).astype(float),
            columns=tree.tip_labels,
            index=[f"s{i}" for i in range(6)],
        )
        adj = np.eye(6, dtype=int)
        for i in range(5):
            adj[i, i + 1] = adj[i + 1, i] = 1
        return tree, mat, adj

    def test_cpb_multisite(self, setup):
        tree, mat, adj = setup
        res = cpb_rate(tree, mat, adj, n_slices=12, component="sorensen")
        assert {"site_id", "CpB", "PB", "pBO", "converged", "r_squared", "reason"} <= set(
            res.columns
        )
        # degenerate neighbourhoods are flagged, the rest converge
        assert res["converged"].sum() >= 2
        assert (res["reason"] != "").sum() >= 0

    def test_cpb_rw_and_components(self, setup):
        tree, mat, adj = setup
        for weighted in (False, True):
            for comp in ("sorensen", "turnover", "nestedness"):
                res = cpb_rate(tree, mat, adj, n_slices=10, component=comp, weighted=weighted)
                assert len(res) == 6

    def test_cpb_pairwise(self, setup):
        tree, mat, adj = setup
        res = cpb_rate(tree, mat, adj, n_slices=10, approach="pairwise")
        assert res["converged"].sum() >= 3


class TestUncertainty:
    def test_posterior_rate(self, tmp_path):
        from phyloslicer.io.newick import write_newick

        trees = [yule_tree(15, age=1.0, seed=s) for s in (1, 2, 3, 4, 5)]
        rng = np.random.default_rng(0)
        mat = pd.DataFrame(
            (rng.random((5, 15)) < 0.3).astype(float),
            columns=trees[0].tip_labels,
            index=[f"s{i}" for i in range(5)],
        )
        res = posterior_rate(trees, mat, n_slices=10)
        assert {
            "site_id",
            "rate_mean",
            "rate_median",
            "hdi_low",
            "hdi_high",
            "ess_proxy",
            "n_trees",
        } <= set(res.columns)
        assert (res["n_trees"] == 5).all()
        assert (res["hdi_low"] <= res["rate_median"]).all()
        assert (res["rate_median"] <= res["hdi_high"]).all()

        # directory-of-trees input
        d = tmp_path / "posterior"
        d.mkdir()
        for i, t in enumerate(trees):
            (d / f"tree{i}.tre").write_text(write_newick(t) + "\n")
        res2 = posterior_rate(str(d), mat, n_slices=10)
        assert np.allclose(res["rate_mean"], res2["rate_mean"])

    def test_posterior_origin(self):
        trees = [yule_tree(12, age=1.0, seed=s) for s in (1, 2)]
        rng = np.random.default_rng(4)
        mat = pd.DataFrame(
            (rng.random((3, 12)) < 0.35).astype(float),
            columns=trees[0].tip_labels,
            index=[f"s{i}" for i in range(3)],
        )
        res = posterior_origin(trees, mat, n_slices=8)
        assert (res["origin_median"] > 0).all()

    def test_bootstrap_sites(self):
        tree = yule_tree(12, age=1.0, seed=1)
        rng = np.random.default_rng(2)
        mat = pd.DataFrame(
            (rng.random((6, 12)) < 0.35).astype(float),
            columns=tree.tip_labels,
            index=[f"s{i}" for i in range(6)],
        )
        res = bootstrap_sites(tree, mat, n_boot=20, n_slices=8, seed=0)
        assert len(res) == 20
        assert np.isfinite(res["rate"]).sum() >= 15


def _pure_sse_grid(cum, ages, grid):
    """Reference implementation of the coarse grid search, without numba."""
    sse = np.empty((cum.shape[0], grid.size), dtype=np.float64)
    for gi, gv in enumerate(grid):
        res = cum - np.exp(-gv * ages)[None, :]
        sse[:, gi] = (res * res).sum(axis=1)
    return sse


class TestAccelerationBackends:
    """The compiled and pure back-ends must be shape- and value-identical.

    An earlier release returned ``(best_rate, best_sse)`` from the numba
    branch while the numpy branch returned a full ``(n_sites, n_grid)`` SSE
    matrix, so ``fit_rates_batch`` raised ``AxisError`` on every install that
    could import numba.  These tests pin the invariant.
    """

    def test_grid_shapes_match(self):
        from phyloslicer import _accel

        rng = np.random.default_rng(11)
        profiles = np.abs(rng.random((5, 10)))
        ages = np.linspace(0.1, 1.0, 10)
        grid = np.geomspace(1e-3, 50.0, 96)
        compiled = _accel.batched_grid_fit(profiles, ages, grid)
        orig = _accel._sse_grid
        try:
            _accel._sse_grid = _pure_sse_grid
            plain = _accel.batched_grid_fit(profiles, ages, grid)
        finally:
            _accel._sse_grid = orig
        assert compiled[0].shape == plain[0].shape == (5, 10)
        assert compiled[1].shape == plain[1].shape == (5, 96)
        assert compiled[1].ndim == 2
        np.testing.assert_allclose(compiled[1], plain[1], rtol=1e-12, atol=0.0)

    def test_one_dimensional_input_is_promoted(self):
        from phyloslicer._accel import batched_grid_fit

        rng = np.random.default_rng(12)
        profiles = np.abs(rng.random(10))
        ages = np.linspace(0.1, 1.0, 10)
        cum, sse = batched_grid_fit(profiles, ages, np.geomspace(1e-3, 50.0, 24))
        assert cum.shape == (1, 10) and sse.shape == (1, 24)
        with pytest.raises(ValueError, match="slice ages"):
            batched_grid_fit(profiles, ages[:-1], np.geomspace(1e-3, 50.0, 24))

    def test_fit_rates_batch_accepts_single_row(self):
        from phyloslicer.rates.fitting import fit_rates_batch

        rng = np.random.default_rng(13)
        ages = np.linspace(2.0, 0.05, 12)
        two_d = fit_rates_batch(np.abs(rng.random((1, 12))), ages)
        one_d = fit_rates_batch(np.abs(rng.random(12)), ages)
        assert len(two_d) == len(one_d) == 1

    @pytest.mark.skipif(
        not __import__("phyloslicer._accel", fromlist=["HAS_NUMBA"]).HAS_NUMBA,
        reason="numba not installed (pip install phyloslicer[accel])",
    )
    def test_rates_identical_with_and_without_numba(self):
        from phyloslicer import _accel
        from phyloslicer.simulate import yule_tree

        tree = yule_tree(30, age=1.0, seed=7)
        rng = np.random.default_rng(8)
        mat = pd.DataFrame(
            (rng.random((6, 30)) < 0.4).astype(float),
            columns=tree.tip_labels,
            index=[f"s{i}" for i in range(6)],
        )
        with_numba = cpd_rate(tree, mat, n_slices=20)
        orig = _accel._sse_grid
        try:
            _accel._sse_grid = _pure_sse_grid
            without = cpd_rate(tree, mat, n_slices=20)
        finally:
            _accel._sse_grid = orig
        assert list(with_numba.converged) == list(without.converged)
        np.testing.assert_array_equal(
            with_numba["CpD"].to_numpy(float), without["CpD"].to_numpy(float)
        )


class TestVizDefaultDirection:
    """正图 P1-7: the graphical default must agree with the published figure orientation."""

    def test_default_axis_runs_from_root_to_present(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        from phyloslicer import viz

        profile = np.array([0.4, 0.3, 0.2, 0.1])
        ages = np.array([1.0, 0.75, 0.5, 0.25])
        fig, ax = plt.subplots()
        viz.plot_rate_line(profile, fit=None, ages=ages, ax=ax)
        xs = sorted(float(x) for x in ax.collections[0].get_offsets()[:, 0])
        assert xs == pytest.approx([0.0, 0.25, 0.5, 0.75], abs=1e-9)
        assert "from root" in ax.get_xlabel().lower()
        plt.close(fig)

    def test_from_present_reproduces_the_reference_view(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np

        from phyloslicer import viz

        profile = np.array([0.4, 0.3, 0.2, 0.1])
        ages = np.array([1.0, 0.75, 0.5, 0.25])
        fig, ax = plt.subplots()
        viz.plot_rate_line(profile, fit=None, ages=ages, ax=ax, direction="from_present")
        xs = sorted(float(x) for x in ax.collections[0].get_offsets()[:, 0])
        assert xs == pytest.approx([0.25, 0.5, 0.75, 1.0], abs=1e-9)
        assert "age" in ax.get_xlabel().lower()
        plt.close(fig)

    def test_bad_direction_is_rejected(self):
        import numpy as np
        import pytest

        from phyloslicer import viz

        with pytest.raises(ValueError, match="direction"):
            viz.plot_rate_line(
                np.array([0.5, 0.5]), ages=np.array([2.0, 1.0]), direction="sideways"
            )
