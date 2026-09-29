"""Regression tests for bugs found and fixed in the 0.1.0 post-review sweep.

Each test pins one concrete defect: the semantic slicing edge case, the
beta-diversity normalisation errors, the dataset misalignment, the
fitting-diagnostics weaknesses, and the I/O robustness gaps.
"""

from __future__ import annotations

import itertools
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phyloslicer.core import TreeArray
from phyloslicer.datasets import sim_passerines
from phyloslicer.indices import pb, pb_per_slice, r_compatible_profile
from phyloslicer.io import (
    align_tree_matrix,
    parse_newick,
    read_newick,
    read_nexus,
    write_newick,
)
from phyloslicer.rates.api import cpd_rate, rate_profile, sensitivity
from phyloslicer.rates.fitting import fit_rate, fit_rates_batch
from phyloslicer.simulate import make_dataset, yule_tree
from phyloslicer.slicing import slice_pieces, slice_rootward, slice_tipward
from phyloslicer.uncertainty import bootstrap_sites, posterior_origin, posterior_rate

REPO_ROOT = Path(__file__).resolve().parents[1]
REFERENCE_DIR = REPO_ROOT / "validation" / "reference"


# ---------------------------------------------------------------------- #
# fixtures / helpers
# ---------------------------------------------------------------------- #
@pytest.fixture
def beta_setup():
    tree = yule_tree(10, age=1.0, seed=5)
    rng = np.random.default_rng(2)
    mat = pd.DataFrame(
        (rng.random((6, 10)) < 0.4).astype(float),
        columns=tree.tip_labels,
        index=[f"s{i}" for i in range(6)],
    )
    adj = np.eye(6, dtype=int)
    for i in range(5):
        adj[i, i + 1] = adj[i + 1, i] = 1
    return tree, mat, adj


def _queen_adj(n_sites: int) -> np.ndarray:
    """The queen-grid adjacency ``generate_reference.R`` used for the scenarios.

    Re-implemented from ``validation/r/generate_reference.R:69-80`` (cell order
    ``x = rep(1:side, times = side)``, ``y = rep(1:side, each = side)``) so the
    beta tests run on the same neighbourhoods the reference CSVs were built on.
    """
    side = int(np.ceil(np.sqrt(n_sites)))
    cells = np.array([(x, y) for y in range(1, side + 1) for x in range(1, side + 1)])
    cells = cells[:n_sites].astype(float)
    adj = np.zeros((n_sites, n_sites), dtype=int)
    for i in range(n_sites):
        near = (np.abs(cells[:, 0] - cells[i, 0]) <= 1) & (np.abs(cells[:, 1] - cells[i, 1]) <= 1)
        adj[i, near] = 1
    return adj


@pytest.fixture(scope="session")
def bd20_random():
    """Committed treesliceR cross-validation scenario: 20 tips, 30 sites.

    The inputs the reference CSVs were produced from, so a profile property
    checked here is checked on the same neighbourhood structure the R side saw
    (30 slices, queen grid) rather than on a toy 6-site stack.
    """
    tree = read_newick(REFERENCE_DIR / "bd20__tree.tre")
    mat = pd.read_csv(REFERENCE_DIR / "bd20__random__mat.csv", index_col=0)
    mat.columns = mat.columns.astype(str)
    mat = (mat > 0).astype(float)
    tree, mat = align_tree_matrix(tree, mat)
    return tree, mat, _queen_adj(mat.shape[0])


class _CorruptTree(TreeArray):
    """Passes iter_trees' isinstance gate but raises on first geometric use,
    simulating a corrupt posterior draw."""

    @property
    def root_age(self):
        raise ValueError("corrupt posterior draw")


def _path_edge_set(tree, tip) -> set[int]:
    parent = tree.parent_of()
    edges: set[int] = set()
    node = int(tip)
    while node in parent:
        for e, (p, c) in enumerate(tree.edges.tolist()):
            if c == node:
                edges.add(e)
                node = p
                break
    return edges


# ---------------------------------------------------------------------- #
# slicing edge cases
# ---------------------------------------------------------------------- #
class TestSliceZero:
    def test_rootward_zero_is_empty(self, balanced4):
        """slice_rootward(t=0) keeps the youngest 0 Ma -> an empty slice."""
        s = slice_rootward(balanced4, 0.0)
        assert s.pd == pytest.approx(0.0)
        assert (np.asarray(s.lengths) == 0).all()

    def test_tipward_zero_is_full_tree(self, balanced4):
        s = slice_tipward(balanced4, 0.0)
        assert s.pd == pytest.approx(balanced4.total_pd())

    def test_rootward_full_age_is_full_tree(self, balanced4):
        T = balanced4.root_age
        s = slice_rootward(balanced4, T)
        assert s.pd == pytest.approx(balanced4.total_pd())

    def test_root_tipward_partition_at_zero(self, random_tree_20):
        t = 0.0
        young = slice_rootward(random_tree_20, t).lengths
        old = slice_tipward(random_tree_20, t).lengths
        assert np.allclose(young + old, random_tree_20.lengths)

    def test_root_age_uses_deepest_tip(self):
        """root_age must be the deepest tip depth, not the last tip's."""
        tree = yule_tree(6, age=1.0, seed=4)
        lengths = tree.lengths.copy()
        tip_edge = int(np.flatnonzero(tree.edges[:, 1] == 0)[0])
        lengths[tip_edge] += 0.5
        stretched = TreeArray(list(tree.tip_labels), tree.edges.copy(), lengths)
        assert stretched.root_age == pytest.approx(
            float(stretched.node_depths()[: stretched.n_tips].max())
        )


# ---------------------------------------------------------------------- #
# beta-diversity normalisation
# ---------------------------------------------------------------------- #
class TestBetaNormalisation:
    def test_rate_profile_pb_not_double_normalised(self, beta_setup):
        """rate_profile('PB') must equal pb_per_slice values (already relative)."""
        tree, mat, adj = beta_setup
        out = rate_profile(tree, mat, adj=adj, n_slices=6, index="PB")
        stack = slice_pieces(tree, n=6)
        direct = pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite")
        assert np.allclose(out["values"], direct["values"], equal_nan=True)
        finite = np.isfinite(out["values"])
        assert finite.any()
        assert np.nanmax(np.abs(out["values"][finite])) <= 1.0 + 1e-12

    def test_all_beta_profiles_partition_the_neighbourhood_beta(self, bd20_random):
        """All six component x approach profiles sum to 1 to 1e-8 (review, A2).

        ``min(b, c)`` is not additive over slices, so the relative profile has to
        be normalised by the sum of its own shares rather than by the whole-tree
        index; :func:`pb_per_slice` does that, which is why the tolerance here is
        a hard ``1e-8`` and not the ``0.05`` this test used to carry.  The old
        version ran on a 6-site / 5-slice toy stack where the deviation happened
        to be small; measured on the committed ``bd20__random`` scenario at 30
        slices the row sums are 1 within 4.5e-16 for every combination.
        """
        tree, mat, adj = bd20_random
        stack = slice_pieces(tree, n=30, criterion="time")
        checked = 0
        for component in ("sorensen", "turnover", "nestedness"):
            for approach in ("multisite", "pairwise"):
                res = pb_per_slice(tree, mat, adj, stack, component, approach)
                finite = np.isfinite(res["total"]) & np.isfinite(res["values"]).all(axis=1)
                assert finite.sum() == mat.shape[0], (component, approach, res["status"][:3])
                sums = res["values"][finite].sum(axis=1)
                assert np.allclose(sums, 1.0, atol=1e-8), (component, approach)
                assert np.abs(sums - 1.0).max() < 1e-12, (component, approach)
                checked += 1
        assert checked == 6

    def test_profile_and_total_share_the_same_quantities(self, bd20_random):
        """``total`` must stay the whole-tree index, identical to :func:`pb`.

        A profile could be made to sum to 1 by rescaling the total as well; the
        assertion that the reported total is *not* rescaled - it reproduces the
        independent no-slicing implementation to 1e-12 - is what stops the
        normalisation from being a cosmetic change.
        """
        tree, mat, adj = bd20_random
        stack = slice_pieces(tree, n=30, criterion="time")
        for component in ("sorensen", "turnover", "nestedness"):
            for approach in ("multisite", "pairwise"):
                for domain in ("mixed", "paired"):
                    res = pb_per_slice(
                        tree,
                        mat,
                        adj,
                        stack,
                        component,
                        approach,
                        multisite_domain=domain,
                    )
                    ref = pb(tree, mat, adj, component, approach, multisite_domain=domain)[
                        "PB"
                    ].to_numpy()
                    ok = np.isfinite(res["total"]) & np.isfinite(ref)
                    assert ok.any(), (component, approach, domain)
                    assert np.allclose(res["total"][ok], ref[ok], atol=1e-12), (
                        component,
                        approach,
                        domain,
                    )

    def test_r_compatible_profile_keeps_the_unnormalised_curve(self, bd20_random):
        """The treesliceR-facing curve still breaks the partition property.

        Pinned so the previous test cannot quietly be satisfied by changing
        *both* normalisations: :func:`r_compatible_profile` must keep summing to
        0.92-1.0 (turnover) and overshooting to 1.54 (nestedness) on this
        scenario, exactly as the reference implementation does.
        """
        tree, mat, adj = bd20_random
        stack = slice_pieces(tree, n=30, criterion="time")
        for component in ("turnover", "nestedness"):
            for approach in ("multisite", "pairwise"):
                rc = r_compatible_profile(tree, mat, adj, stack, component, approach)["values"]
                finite = np.isfinite(rc).all(axis=1)
                dev = np.abs(rc[finite].sum(axis=1) - 1.0)
                assert dev.max() > 1e-3, (component, approach)
                profile = pb_per_slice(tree, mat, adj, stack, component, approach)["values"]
                assert np.abs(profile[finite].sum(axis=1) - 1.0).max() < 1e-12
        for component in ("sorensen",):
            rc = r_compatible_profile(tree, mat, adj, stack, component, "multisite")["values"]
            assert np.abs(rc.sum(axis=1) - 1.0).max() < 1e-12  # b + c is additive

    def test_nestedness_multisite_matches_r_arithmetic(self, beta_setup):
        """Series reproduces treesliceR's literal multisite nestedness math."""
        tree, mat, adj = beta_setup
        stack = slice_pieces(tree, n=5)
        res = pb_per_slice(tree, mat, adj, stack, "nestedness", "multisite")
        C = stack.contribution_matrix()
        k = stack.n_slices
        verified = 0
        for focal in range(mat.shape[0]):
            if not np.isfinite(res["total"][focal]):
                continue
            members = np.flatnonzero(adj[focal] == 1)
            paths = [
                _path_edge_set(tree, t)
                for m in members
                for t in np.flatnonzero(mat.iloc[m].to_numpy() > 0)
            ]
            # regroup: one path set per member
            paths = []
            for m in members:
                present = np.flatnonzero(mat.iloc[m].to_numpy() > 0)
                edges: set[int] = set()
                for t in present:
                    edges |= _path_edge_set(tree, t)
                paths.append(edges)

            def full_pd(es: set) -> float:
                return float(sum(tree.lengths[e] for e in es))

            def slice_pd(es: set, j: int) -> float:
                return float(sum(C[e, j] for e in es))

            pd_single = np.array([full_pd(p) for p in paths])
            pairs = list(itertools.combinations(range(members.size), 2))
            pw_pd = np.array([full_pd(paths[i] | paths[j]) for i, j in pairs])
            ii = np.array([i for i, _ in pairs])
            jj = np.array([j for _, j in pairs])
            # R's recycled (P x 2) matrix arithmetic: row sum = 2*pw - pd_i - pd_j
            bc = float((2 * pw_pd - (pd_single[ii] + pd_single[jj])).sum())
            mn = float((pw_pd - np.maximum(pd_single[ii], pd_single[jj])).sum())
            a_tot = float(pd_single.sum() - full_pd(set().union(*paths)))
            nest = bc / (2 * a_tot + bc) - mn / (a_tot + mn)
            series = []
            for j in range(k):
                piece = np.array([slice_pd(p, j) for p in paths])
                pw_piece = np.array(
                    [slice_pd(paths[i] | paths[jj[n]], j) for n, (i, _) in enumerate(pairs)]
                )
                r_bc = float((2 * pw_piece - (piece[ii] + piece[jj])).sum())
                r_mn = float((pw_piece - np.maximum(piece[ii], piece[jj])).sum())
                series.append(r_bc / (2 * a_tot + bc) - r_mn / (a_tot + mn))
            expected = np.array(series) / nest
            assert np.allclose(res["values"][focal], expected, atol=1e-10)
            verified += 1
        assert verified >= 2

    def test_pb_pairwise_ignores_undefined_pairs(self, beta_setup):
        """pb() pairwise should nanmean over pairs, not propagate one NaN."""
        tree, mat, adj = beta_setup
        mat2 = mat.copy()
        mat2.iloc[1, :] = 0.0  # empty site -> some pair ratios 0/0 = NaN
        res = pb(tree, mat2, adj, component="sorensen", approach="pairwise")
        assert np.isfinite(res["PB"].to_numpy()).any()


# ---------------------------------------------------------------------- #
# fitting diagnostics
# ---------------------------------------------------------------------- #
class TestFittingDiagnostics:
    def test_optimizer_exception_keeps_best_start(self):
        """A failing multistart must not clobber an earlier successful fit."""
        import phyloslicer.rates.fitting as fitting

        ages = np.linspace(2.0, 0.1, 20)
        cum = np.exp(-0.9 * ages)
        profile = np.diff(np.concatenate([[0.0], cum]))
        orig = fitting.least_squares
        calls = {"n": 0}

        def flaky(*a, **k):
            calls["n"] += 1
            if calls["n"] >= 2:
                raise RuntimeError("boom")
            return orig(*a, **k)

        fitting.least_squares = flaky
        try:
            fit = fit_rate(profile, ages, multistart=4)
        finally:
            fitting.least_squares = orig
        assert fit.converged
        assert np.isfinite(fit.rate)

    def test_batch_converged_implies_robust_agreement(self):
        """Whenever the vectorised kernel certifies convergence, its rate
        must agree with the robust per-row optimiser."""
        rng = np.random.default_rng(0)
        ages = np.linspace(2.0, 0.05, 12)
        checked = 0
        for _ in range(30):
            p = rng.random(12)
            p[0] += 0.5
            fb = fit_rates_batch(p[None, :], ages)[0]
            if fb.converged:
                fr = fit_rate(p, ages)
                assert abs(fb.rate - fr.rate) < 1e-3
                checked += 1
        assert checked > 0

    def test_degenerate_profile_reason(self):
        fit = fit_rate(np.zeros(4), np.linspace(2, 0.1, 4))
        assert not fit.converged
        assert fit.reason


# ---------------------------------------------------------------------- #
# datasets
# ---------------------------------------------------------------------- #
class TestDatasets:
    def test_sim_passerines_tree_matrix_align(self):
        data = sim_passerines(n_sites=12, n_species=8, n_trees=2)
        tree, mat = data["trees"][0], data["mat"]
        tree2, mat2 = align_tree_matrix(tree, mat)
        assert list(mat2.columns) == tree2.tip_labels
        assert mat2.shape[0] == 12

    def test_sim_passerines_end_to_end(self):
        data = sim_passerines(n_sites=9, n_species=8, n_trees=2)
        res = cpd_rate(data["trees"][0], data["mat"], n_slices=5)
        assert np.isfinite(res["CpD"]).all()

    def test_sim_passerines_posterior_pool_consistent(self):
        data = sim_passerines(n_sites=9, n_species=8, n_trees=3)
        labels = {tuple(sorted(t.tip_labels)) for t in data["trees"]}
        assert len(labels) == 1


# ---------------------------------------------------------------------- #
# uncertainty robustness
# ---------------------------------------------------------------------- #
class TestUncertaintyRobustness:
    @pytest.fixture
    def corrupt_pool(self, beta_setup):
        tree, mat, _ = beta_setup
        bad = _CorruptTree(
            tip_labels=list(mat.columns[:2]),
            edges=np.array([[2, 0], [2, 1]]),
            lengths=np.array([0.5, 0.5]),
        )
        return [tree, tree, bad], mat

    def test_posterior_origin_skips_bad_tree(self, corrupt_pool):
        trees, mat = corrupt_pool
        with pytest.warns(UserWarning, match="skipping tree 3"):
            res = posterior_origin(trees, mat, n_slices=5)
        assert len(res) == mat.shape[0]
        assert (res["origin_median"] > 0).all()

    def test_posterior_rate_skips_bad_tree(self, corrupt_pool):
        trees, mat = corrupt_pool
        with pytest.warns(UserWarning, match="skipping tree 3"):
            res = posterior_rate(trees, mat, n_slices=5)
        assert (res["n_trees"] == 2).all()

    def test_bootstrap_sites_aligns_tree_and_matrix(self, beta_setup):
        tree, mat, _ = beta_setup
        mat_ghost = mat.copy()
        mat_ghost["ghost1"] = 1.0  # species absent from the tree
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ref = bootstrap_sites(tree, mat, n_boot=5, n_slices=6, seed=0)
            got = bootstrap_sites(tree, mat_ghost, n_boot=5, n_slices=6, seed=0)
        assert np.allclose(ref["rate"], got["rate"], equal_nan=True)

    def test_posterior_rate_hdi_validation(self, beta_setup):
        tree, mat, _ = beta_setup
        with pytest.raises(ValueError, match="hdi"):
            posterior_rate([tree], mat, n_slices=4, hdi=1.5)


# ---------------------------------------------------------------------- #
# API validation & CLI
# ---------------------------------------------------------------------- #
class TestValidationGaps:
    def test_sensitivity_cpb_requires_adj(self, beta_setup):
        tree, mat, _ = beta_setup
        with pytest.raises(ValueError, match="adjacency"):
            sensitivity(tree, mat, slice_grid=[4], rate="cpb")

    def test_sensitivity_negative_sample_sites(self, beta_setup):
        tree, mat, _ = beta_setup
        with pytest.raises(ValueError, match="sample_sites"):
            sensitivity(tree, mat, slice_grid=[4], rate="cpd", sample_sites=-1)

    def test_validate_command_both_na_is_agreement(self, tmp_path, random_tree_20):
        from phyloslicer.cli.main import main

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(write_newick(random_tree_20) + "\n")
        mat, _ = make_dataset(n_sites=4, n_species=random_tree_20.n_tips, seed=1)
        mat.columns = random_tree_20.tip_labels
        mat.iloc[3, :] = 0.0  # this site is NA on both sides
        sites_path = tmp_path / "sites.csv"
        mat.to_csv(sites_path)
        rates_path = tmp_path / "rates.csv"
        assert (
            main(
                [
                    "rates",
                    str(tree_path),
                    str(sites_path),
                    "--n-slices",
                    "5",
                    "--out",
                    str(rates_path),
                ]
            )
            == 0
        )
        df = pd.read_csv(rates_path)
        ref = pd.DataFrame({"CpD_ref": df["CpD"]})
        ref_path = tmp_path / "ref.csv"
        ref.to_csv(ref_path)
        out = tmp_path / "val.csv"
        # This fixture deliberately uses the package's own column as the "reference": what it
        # pins is the both-NA-is-agreement rule, not equivalence with R.  The manifest says so
        # explicitly, which is the point of the provenance gate - a self-supplied reference can no
        # longer be passed off as an external one silently.
        manifest = write_reference_manifest(tmp_path, ref_path)
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "5",
                    "--out",
                    str(out),
                    "--manifest",
                    str(manifest),
                ]
            )
            == 0
        )
        val = pd.read_csv(out)
        # the empty site is NA on both sides -> agreement, not failure
        assert bool(val["equivalent"].iloc[3]) is True

    def test_slice_command_rejects_conflicting_options(self, tmp_path, random_tree_20):
        from phyloslicer.cli.main import main

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(write_newick(random_tree_20) + "\n")
        rc = main(
            [
                "slice",
                str(tree_path),
                "--rootward",
                "1.0",
                "--n-slices",
                "3",
                "--out",
                str(tmp_path / "x.tre"),
            ]
        )
        assert rc == 2


# ---------------------------------------------------------------------- #
# I/O robustness
# ---------------------------------------------------------------------- #
class TestIoRobustness:
    def test_nexus_multiline_tree_statements(self, tmp_path):
        text = (
            "#NEXUS\n\nBEGIN TREES;\n"
            "    TREE t1 =\n        ((A:1,B:1):1,(C:1,D:1):1);\n"
            "    TREE t2 = ((A:1,B:1):2,C:2);\n"
            "END;\n"
        )
        p = tmp_path / "t.nex"
        p.write_text(text)
        trees = read_nexus(p)
        assert len(trees) == 2
        assert trees[0].n_tips == 4
        # ((A:1,B:1):2,C:2): tips sit at depth 2 + 1 = 3
        assert trees[1].root_age == pytest.approx(3.0)

    def test_newick_roundtrip_preserves_pd(self, random_tree_20):
        back = parse_newick(write_newick(random_tree_20))
        assert back.total_pd() == pytest.approx(random_tree_20.total_pd())


def write_reference_manifest(tmp_path, ref_path, rate_column="CpD_ref"):
    """Register an ad-hoc reference CSV through the package's own manifest writer (A13)."""
    import hashlib
    import os

    import pandas as _pd

    from phyloslicer.cli.main import format_reference_manifest

    blob = open(ref_path, "rb").read()
    # the comparator reads references with the first column as the index;
    # record the shape the way it will see it, or A13 reports a false mismatch
    frame = _pd.read_csv(ref_path, index_col=0)
    text = format_reference_manifest(
        {
            "written": "test",
            "generator": "tests (synthetic fixture)",
            "reference_package": "none - this fixture pins the NA rule, not R equivalence",
        },
        {
            os.path.basename(str(ref_path)): {
                "sha256": hashlib.sha256(blob).hexdigest(),
                "rows": str(frame.shape[0]),
                "cols": str(frame.shape[1]),
                "orientation": "site-by-index",
                "columns_pattern": ",".join(map(str, frame.columns[:3])),
                "rate_column": rate_column,
            }
        },
    )
    path = tmp_path / "manifest.txt"
    path.write_text(text)
    return path
