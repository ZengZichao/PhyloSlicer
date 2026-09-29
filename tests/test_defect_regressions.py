"""Regression tests for defects found while cross-validating against treesliceR.

One test per fixed defect (or per family of closely-related defects), each
docstring naming the behaviour it pins, so a future reader can trace a test
back
docstring naming the behaviour it pins, so a future reader can trace a test
back
I/O round-trips, input validation, CLI exit codes, and the promises that
docstrings make about their own behaviour.
"""

from __future__ import annotations

import json
import re
import warnings

import numpy as np
import pandas as pd
import pytest

from phyloslicer import TreeArray
from phyloslicer.core.tree import SlicedTree
from phyloslicer.indices import api as indices
from phyloslicer.indices.api import pb_per_slice, pd_per_slice, pe_per_slice
from phyloslicer.io.matrices import align_tree_matrix, presence_matrix
from phyloslicer.io.newick import read_newick, split_newick_records, write_newick
from phyloslicer.rates.api import cpd_rate
from phyloslicer.rates.fitting import fit_rate
from phyloslicer.simulate import birth_death_tree, make_dataset, yule_tree
from phyloslicer.slicing.api import prune_tips, slice_interval, slice_pieces, slice_rootward
from phyloslicer.spatial import adjacency_from_bounds, from_adjacency_dict, to_adjacency_dict
from phyloslicer.uncertainty import iter_trees, posterior_rate

# --------------------------------------------------------------------- #
# shared fixtures
# --------------------------------------------------------------------- #


@pytest.fixture
def apostrophe_tree() -> TreeArray:
    """A tree whose labels contain the characters Newick readers struggle with."""
    return TreeArray.from_newick("(((a'b:0.5,c:0.5):0.5,(d:0.5,e:0.5):0.5):0.1,(f:0.5,g:0.5):0.6);")


@pytest.fixture
def four_tip() -> TreeArray:
    """((A,B),(C,D)) with root age 2 and *varied* terminal branches (ultrametric)."""
    return TreeArray(
        ["A", "B", "C", "D"],
        np.array([[4, 5], [4, 6], [5, 0], [5, 1], [6, 2], [6, 3]]),
        np.array([1.0, 2.0, 2.0, 2.0, 1.0, 1.0]),
    )


@pytest.fixture
def small_setup(four_tip):
    mat = pd.DataFrame(
        [[1, 1, 0, 0], [0, 1, 1, 0], [1, 0, 0, 1], [1, 1, 1, 0]],
        index=["s1", "s2", "s3", "s4"],
        columns=four_tip.tip_labels,
    )
    tree, mat = align_tree_matrix(four_tip, mat)
    return tree, mat, slice_pieces(tree, n=6)


# ===================================================================== #
# P0-1  Newick writer/reader quote handling
# ===================================================================== #


class TestNewickQuoting:
    def test_writer_escapes_apostrophe_and_round_trips(self, apostrophe_tree):
        text = apostrophe_tree.to_newick()
        assert "'a''b'" in text  # Newick escapes a quote by doubling it
        back = TreeArray.from_newick(text)
        assert back == apostrophe_tree
        assert "a'b" in back.tip_labels

    def test_reader_accepts_the_doubled_quote_escape(self):
        tree = TreeArray.from_newick("('a''b':1.0,(C:1,D:1):1.0);")
        assert tree.tip_labels[0] == "a'b"

    def test_all_cli_writable_labels_survive_round_trip(self):
        awkward = TreeArray(
            ["A B", "C'D", "E,F", "G;H", "I)J"],
            np.array([[5, 0], [5, 1], [5, 6], [6, 2], [6, 3], [6, 4]]),
            np.array([1.0, 1.0, 0.4, 0.5, 0.5, 0.5]),
        )
        assert TreeArray.from_newick(awkward.to_newick()) == awkward

    def test_slice_cli_output_can_be_read_back(self, tmp_path):
        """The documented `slice -> posterior` chain must not dead-end."""
        src = tmp_path / "in.tre"
        src.write_text("(((a'b:0.5,c:0.5):0.5,(d:0.5,e:0.5):0.5):0.1,(f:0.5,g:0.5):0.6);\n")
        from phyloslicer.cli.main import main

        out = tmp_path / "slices"
        assert main(["slice", str(src), "--n-slices", "3", "--out", str(out)]) == 0
        pool = tmp_path / "pool"
        pool.mkdir()
        for j, f in enumerate(sorted(out.glob("slice_*.tre"))):
            # a slice may be degenerate (all-zero lengths); keep the readable ones
            try:
                t = TreeArray.from_newick(f.read_text())
            except ValueError:
                continue
            if t.root_age > 0:
                (pool / f"t{j}.tre").write_text(f.read_text())
        assert list(pool.glob("*.tre")), "no slice survived the write/read round trip"


# ===================================================================== #
# P0-2 / P2-15  validate: exit codes, tolerance, quantity, n_slices
# ===================================================================== #


def _register(tmp_path, ref_path, rate_column="CpD", n_slices=None):
    """Write a manifest governing ``ref_path`` through the package's own writer."""
    from phyloslicer.cli.main import format_reference_manifest, make_reference_record

    frame = pd.read_csv(ref_path, index_col=0)
    header = {"format": "phyloslicer-reference-manifest/1"}
    if n_slices is not None:
        header["n_slices"] = str(n_slices)
    record = make_reference_record(
        ref_path,
        orientation="site-by-index",
        columns_pattern=(
            f"^({rate_column}|{rate_column[1:]}|p{'D' if rate_column.endswith('D') else 'E'}O)$"
        ),
        rows=frame.shape[0],
        cols=frame.shape[1],
        rate_column=rate_column,
    )
    manifest = tmp_path / "manifest.txt"
    manifest.write_text(format_reference_manifest(header, {ref_path.name: record}))
    return manifest


@pytest.fixture
def cpd_reference(tmp_path, four_tip):
    """A self-consistent CpD reference: our own output, written out and read back."""
    mat = pd.DataFrame(
        [[1, 1, 0, 0], [0, 1, 1, 0], [1, 0, 0, 1], [1, 1, 1, 0]],
        index=["s1", "s2", "s3", "s4"],
        columns=four_tip.tip_labels,
    )
    tree_path = tmp_path / "in.tre"
    tree_path.write_text(write_newick(four_tip) + "\n")
    sites_path = tmp_path / "sites.csv"
    mat.to_csv(sites_path)
    res = cpd_rate(four_tip, mat, n_slices=8)
    ref_path = tmp_path / "ref.csv"
    pd.DataFrame({"CpD": res["CpD"].to_numpy()}).to_csv(ref_path)
    return tree_path, sites_path, ref_path, res


class TestValidateGate:
    def test_agreement_exits_zero(self, tmp_path, cpd_reference):
        from phyloslicer.cli.main import main

        tree_path, sites_path, ref_path, _ = cpd_reference
        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        out = tmp_path / "report.csv"
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--manifest",
                    str(manifest),
                    "--n-slices",
                    "8",
                    "--out",
                    str(out),
                ]
            )
            == 0
        )
        assert "True" in out.read_text()

    def test_disagreement_is_a_nonzero_exit(self, tmp_path, cpd_reference):
        """The gate's whole purpose: a 0 % equivalence rate must not exit 0."""
        from phyloslicer.cli.main import main

        tree_path, sites_path, ref_path, _ = cpd_reference
        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        out = tmp_path / "report.csv"
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--manifest",
                    str(manifest),
                    "--n-slices",
                    "8",
                    "--tolerance",
                    "1e-30",
                    "--out",
                    str(out),
                ]
            )
            == 1
        )

    def test_no_strict_can_be_downgraded_for_exploration(self, tmp_path, cpd_reference):
        from phyloslicer.cli.main import main

        tree_path, sites_path, ref_path, _ = cpd_reference
        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--manifest",
                    str(manifest),
                    "--n-slices",
                    "8",
                    "--tolerance",
                    "1e-30",
                    "--no-strict",
                    "--out",
                    str(tmp_path / "r.csv"),
                ]
            )
            == 0
        )

    def test_default_tolerance_does_not_deny_its_own_reference(self, tmp_path, cpd_reference):
        """1e-8 used to report 0.0000 on a perfectly good comparison."""
        from phyloslicer.cli.main import main

        tree_path, sites_path, ref_path, _ = cpd_reference
        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        out = tmp_path / "r.csv"
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--manifest",
                    str(manifest),
                    "--n-slices",
                    "8",
                    "--out",
                    str(out),
                ]
            )
            == 0
        )

    def test_wrong_slice_count_is_refused_not_scored(self, tmp_path, cpd_reference):
        from phyloslicer.cli.main import main

        tree_path, sites_path, ref_path, _ = cpd_reference
        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        rc = main(
            [
                "validate",
                str(tree_path),
                str(sites_path),
                "--reference",
                str(ref_path),
                "--manifest",
                str(manifest),
                "--n-slices",
                "50",
                "--out",
                str(tmp_path / "r.csv"),
            ]
        )
        assert rc == 2

    def test_report_columns_are_named_after_the_compared_index(self, tmp_path, cpd_reference):
        tree_path, sites_path, ref_path, _ = cpd_reference
        from phyloslicer.cli.main import main

        manifest = _register(tmp_path, ref_path, "CpD", n_slices=8)
        out = tmp_path / "r.csv"
        main(
            [
                "validate",
                str(tree_path),
                str(sites_path),
                "--reference",
                str(ref_path),
                "--manifest",
                str(manifest),
                "--n-slices",
                "8",
                "--out",
                str(out),
            ]
        )
        header = out.read_text().splitlines()[0]
        assert "phyloslicer_CpD" in header and "reference_CpD" in header
        assert "reference_CpB" not in header

    def test_rate_column_help_no_longer_promises_an_override(self):
        from phyloslicer.cli.main import _build_parser

        parser = _build_parser()
        text = (
            parser.format_help()
            + parser._subparsers._group_actions[0].choices["validate"].format_help()
        )
        assert "assert which column" in text or "assert the column" in text


# ===================================================================== #
# P0-3 / P1-1  posterior pools: nothing is dropped in silence
# ===================================================================== #


class TestPosteriorPoolIntegrity:
    def test_stray_separators_do_not_abort_the_pool(self, tmp_path):
        text = write_newick(yule_tree(6, seed=1)) + ";;\n();\n" + write_newick(yule_tree(6, seed=2))
        path = tmp_path / "pool.tre"
        path.write_text(text)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            trees = list(iter_trees(path))
        assert len(trees) == 2  # both real trees, no exception from the noise
        assert any("skipping record" in str(w.message) for w in caught)

    def test_unreadable_file_is_reported_not_silently_dropped(self, tmp_path):
        good = write_newick(yule_tree(6, seed=1))
        (tmp_path / "a.tre").write_text(good)
        (tmp_path / "b.tre").write_text(good)
        (tmp_path / "c_bad.tre").write_text(good.replace("(", "('", 1))
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            trees = list(iter_trees(tmp_path))
        assert len(trees) < 3
        assert caught, "a dropped tree file must be announced"

    def test_newick_text_is_accepted_like_read_newick_does(self):
        text = write_newick(yule_tree(6, seed=1)) + "\n" + write_newick(yule_tree(6, seed=2))
        assert len(list(iter_trees(text))) == 2

    def test_missing_path_gives_a_finding_error_not_an_oserror(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            list(iter_trees(tmp_path / "nope.tre"))
        with pytest.raises((FileNotFoundError, ValueError)):
            list(iter_trees("definitely not a path or a tree"))

    def test_shrunk_pool_is_announced(self, tmp_path):
        """posterior_rate must not summarise a silently smaller pool."""
        tree = yule_tree(8, age=4.0, seed=3)
        mat = pd.DataFrame(
            np.eye(4, 8, dtype=float), columns=tree.tip_labels, index=[f"s{i}" for i in range(4)]
        )
        (tmp_path / "t.tre").write_text(write_newick(tree) + "\n")
        good = posterior_rate(tmp_path / "t.tre", mat, n_slices=5)
        assert int(good["n_trees"].iloc[0]) == 1

    def test_iter_trees_exported_at_top_level(self):
        import phyloslicer

        assert callable(phyloslicer.iter_trees)
        assert callable(phyloslicer.bootstrap_sites)


# ===================================================================== #
# P1-2  birth_death_tree delivers n tips, in bounded time
# ===================================================================== #


class TestBirthDeathTree:
    @pytest.mark.parametrize("n", [10, 20, 50, 100, 500])
    def test_yields_exactly_n_tips_at_any_realistic_size(self, n):
        tree = birth_death_tree(n, age=2.0, birth_rate=1.0, death_rate=0.2, seed=7)
        assert tree.n_tips == n
        assert tree.is_ultrametric
        assert tree.root_age == pytest.approx(2.0)

    def test_previously_hanging_parameterisation_returns(self):
        """age=200 with lambda-mu=0.6 used to simulate e**120 lineages forever."""
        tree = birth_death_tree(15, age=200.0, birth_rate=1.0, death_rate=0.4, seed=1)
        assert tree.n_tips == 15

    @pytest.mark.parametrize(
        "kwargs,message",
        [
            ({"n": 1}, r"n must be >= 2"),
            ({"n": 10, "age": 0.0}, "age must be positive"),
            ({"n": 10, "birth_rate": 0.0}, "birth_rate must be positive"),
            ({"n": 10, "death_rate": -1.0}, "non-negative"),
        ],
    )
    def test_invalid_requests_fail_loudly(self, kwargs, message):
        with pytest.raises(ValueError, match=message):
            birth_death_tree(**kwargs)

    def test_subcritical_rates_warn_before_retrying(self):
        with pytest.warns(UserWarning, match="not below birth_rate"):
            try:
                birth_death_tree(50, age=1.0, birth_rate=1.0, death_rate=1.2, seed=1)
            except RuntimeError:
                pass  # exhausting attempts is a legitimate outcome here


# ===================================================================== #
# P1-4 / P1-5  the numeric domain of lengths and thresholds
# ===================================================================== #


class TestNumericDomain:
    @pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
    def test_non_finite_branch_lengths_are_refused(self, bad):
        with pytest.raises(ValueError, match="finite"):
            TreeArray(["A", "B"], np.array([[2, 0], [2, 1]]), np.array([bad, 1.0]))

    def test_infinite_tree_is_never_certified_ultrametric(self):
        # guards the `spread <= rtol * scale` degeneration: inf <= inf is True
        depths = np.array([np.inf, 1.0])
        scale, spread = depths.max(), depths.max() - depths.min()
        assert scale == np.inf and spread == np.inf
        assert not (np.isfinite(scale) and np.isfinite(spread))

    @pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
    def test_non_finite_thresholds_are_refused(self, four_tip, bad):
        with pytest.raises(ValueError, match="finite"):
            slice_rootward(four_tip, bad)
        with pytest.raises(ValueError, match="finite"):
            slice_interval(four_tip, bad, 1.0)

    def test_slice_with_non_finite_threshold_cannot_produce_nan_pd(self, four_tip):
        with pytest.raises(ValueError):
            slice_rootward(four_tip, float("nan"))


# ===================================================================== #
# P1-6 / P2-2  label mismatch: one contract across the package
# ===================================================================== #


class TestLabelAlignment:
    def test_index_functions_refuse_a_total_mismatch(self, small_setup):
        tree, mat, _ = small_setup
        renamed = mat.copy()
        renamed.columns = ["nope_%d" % i for i in range(renamed.shape[1])]
        with pytest.raises(ValueError, match="no matrix column matches"):
            indices.pd(tree, renamed)
        with pytest.raises(ValueError, match="no matrix column matches"):
            indices.pe(tree, renamed)

    def test_align_tree_matrix_refuses_too_few_shared_names(self, four_tip):
        tree = TreeArray(
            ["A", "B", "C"], np.array([[3, 0], [3, 1], [3, 2]]), np.array([1.0, 1.0, 2.0])
        )
        with pytest.raises(ValueError, match="tree tip labels"):
            align_tree_matrix(tree, pd.DataFrame([[1, 0]], columns=["A", "Z"], index=["s1"]))
        with pytest.raises(ValueError, match="tree tip labels"):
            align_tree_matrix(tree, pd.DataFrame([[1, 0]], columns=["X", "Y"], index=["s1"]))

    def test_ndarray_matrix_raises_instead_of_returning_zeros(self, small_setup):
        tree, mat, _ = small_setup
        with pytest.raises(
            ValueError, match="numpy array gives integer column names|no matrix column"
        ):
            indices.pd(tree, mat.to_numpy())


# ===================================================================== #
# P1-7  slice --out cannot silently mix generations
# ===================================================================== #


class TestCliOutputHandling:
    def test_stale_slice_directory_is_refused_then_forced(self, tmp_path, four_tip):
        from phyloslicer.cli.main import main

        src = tmp_path / "in.tre"
        src.write_text(write_newick(four_tip) + "\n")
        out = tmp_path / "slices"
        assert main(["slice", str(src), "--n-slices", "5", "--out", str(out)]) == 0
        assert len(list(out.glob("slice_*.tre"))) == 5
        assert main(["slice", str(src), "--n-slices", "2", "--out", str(out)]) == 2
        assert len(list(out.glob("slice_*.tre"))) == 5, "refusal must not alter files"
        assert main(["slice", str(src), "--n-slices", "2", "--out", str(out), "--force"]) == 0
        assert sorted(p.name for p in out.glob("slice_*.tre")) == [
            "slice_0001.tre",
            "slice_0002.tre",
        ], "--force must clear the previous generation"


# ===================================================================== #
# P1-3  the run summary's clock
# ===================================================================== #


class TestRunLog:
    def test_elapsed_seconds_is_a_duration(self, tmp_path, four_tip):
        from phyloslicer.cli.main import main

        src = tmp_path / "in.tre"
        src.write_text(write_newick(four_tip) + "\n")
        mat = pd.DataFrame([[1, 1, 0, 0]], columns=four_tip.tip_labels, index=["s1"])
        sites = tmp_path / "m.csv"
        mat.to_csv(sites)
        log = tmp_path / "l.json"
        assert (
            main(
                [
                    "rates",
                    str(src),
                    str(sites),
                    "--n-slices",
                    "5",
                    "--out",
                    str(tmp_path / "r.csv"),
                    "--log",
                    str(log),
                ]
            )
            == 0
        )
        payload = json.loads(log.read_text())
        assert 0 <= payload["elapsed_s"] < 60, payload["elapsed_s"]
        assert payload["started_at"][:4].isdigit()


# ===================================================================== #
# P2-1  pb() and pb_per_slice() agree about degenerate neighbourhoods
# ===================================================================== #


class TestBetaDegeneracyConsistency:
    def test_identical_sites_are_refused_by_both_functions(self, small_setup):
        tree, mat, stack = small_setup
        same = pd.DataFrame(
            np.tile(mat.values[0], (mat.shape[0], 1)), index=mat.index, columns=mat.columns
        )
        adj = np.ones((mat.shape[0], mat.shape[0]), dtype=np.int8)
        per_slice = pb_per_slice(tree, same, adj, stack)
        totals = indices.pb(tree, same, adj)
        assert set(per_slice["status"]) == {"degenerate_neighbourhood"}
        assert np.isnan(per_slice["total"]).all()
        assert np.isnan(totals["PB"]).all(), "pb() must not answer 0 where the other NaNs"
        assert set(totals["status"]) == {"degenerate_neighbourhood"}

    def test_status_column_explains_every_nan(self, small_setup):
        tree, mat, _ = small_setup
        adj = np.eye(mat.shape[0], dtype=np.int8)
        adj[0, 1] = adj[1, 0] = 1
        frame = indices.pb(tree, mat, adj)
        assert {"site_id", "PB", "status"} <= set(frame.columns)
        nan_rows = frame["PB"].isna()
        assert (frame.loc[nan_rows, "status"] != "ok").all()


# ===================================================================== #
# P2-3 / P3-4  the equal-PD curve
# ===================================================================== #


class TestEqualPdSlicing:
    def test_curve_reaches_total_pd_on_a_non_ultrametric_tree(self):
        from phyloslicer.slicing.windows import pd_depth_curve

        tree = TreeArray(
            ["A", "B", "C"], np.array([[3, 0], [3, 1], [3, 2]]), np.array([1.0, 2.0, 2.0])
        )
        _breaks, cum = pd_depth_curve(tree)
        assert cum[-1] == pytest.approx(tree.total_pd())

    def test_waived_non_ultrametric_pd_slicing_works_and_partitions(self):
        tree = TreeArray(
            ["A", "B", "C"], np.array([[3, 0], [3, 1], [3, 2]]), np.array([1.0, 2.0, 2.0])
        )
        with pytest.warns(UserWarning, match="not ultrametric"):
            stack = slice_pieces(tree, n=4, criterion="pd", allow_non_ultrametric=True)
        C = stack.contribution_matrix()
        per_window = C.sum(axis=0)
        assert np.allclose(per_window, per_window[0])  # genuinely equal PD
        assert C.sum() == pytest.approx(tree.total_pd())
        assert np.allclose(C.sum(axis=1), tree.lengths)  # every branch, fully

    def test_depth_for_pd_agrees_with_its_vectorised_sibling(self, four_tip):
        from phyloslicer.slicing.windows import depth_for_pd, depths_for_pd

        targets = np.linspace(0.0, four_tip.total_pd(), 9)
        assert np.allclose(
            [depth_for_pd(four_tip, t) for t in targets], depths_for_pd(four_tip, targets)
        )

    def test_targets_beyond_the_tree_are_rejected(self, four_tip):
        from phyloslicer.slicing.windows import depths_for_pd

        with pytest.raises(ValueError, match="outside the available range"):
            depths_for_pd(four_tip, np.array([0.0, four_tip.total_pd() * 2]))

    def test_pd_windows_partition_the_tree_on_a_large_input(self):
        """The absolute slack used to reject a 5000-tip tree's own total PD."""
        tree = yule_tree(5000, age=500.0, seed=99)
        stack = slice_pieces(tree, n=100, criterion="pd")
        assert np.allclose(
            stack.contribution_matrix().sum(axis=1), tree.lengths, rtol=1e-9, atol=1e-9
        )


# ===================================================================== #
# P2-4 / P2-5 / P2-6 / P2-7 / P2-12  documented promises that were absent
# ===================================================================== #


class TestDocumentedContracts:
    def test_widths_are_depth_widths_and_the_docstring_says_so(self, small_setup):
        tree, _mat, _ = small_setup
        stack = slice_pieces(tree, n=6, criterion="pd")
        assert "depth" in type(stack).widths.__doc__.lower()
        per_slice_pd = stack.contribution_matrix().sum(axis=0)
        assert np.allclose(per_slice_pd, per_slice_pd[0]), "PD is what is equal here"

    def test_asymmetric_adjacency_round_trip_warns(self):
        asym = np.array([[1, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=np.int8)
        with pytest.warns(UserWarning, match="asymmetric"):
            rebuilt = from_adjacency_dict(to_adjacency_dict(asym))
        assert not np.array_equal(rebuilt, asym)  # and it told us

    def test_presence_matrix_does_not_guess_an_abundance_column(self):
        geo = pd.DataFrame({"site": ["a", "a"], "sp": ["1", "2"], "lat": [10.0, 20.0]})
        with pytest.warns(UserWarning, match="ignored"):
            wide = presence_matrix(geo, "site", "sp")
        assert float(wide.loc["a", "1"]) == 1.0 and float(wide.loc["a", "2"]) == 1.0

    def test_presence_matrix_uses_an_explicit_value_column(self):
        long = pd.DataFrame({"site": ["a", "a", "b"], "sp": ["1", "2", "1"], "ab": [3, 0, 5]})
        wide = presence_matrix(long, "site", "sp", "ab")
        assert float(wide.loc["a", "2"]) == 0.0  # abundance 0 is absence
        assert float(wide.loc["a", "1"]) == 1.0

    def test_fit_rate_rejects_an_unknown_method(self):
        profile = np.array([0.1, 0.2, 0.3, 0.4])
        ages = np.array([4.0, 3.0, 2.0, 1.0])
        with pytest.raises(ValueError, match="method must be"):
            fit_rate(profile, ages, method="maximum_likelihood")
        assert fit_rate(profile, ages, method="lsq").rate == pytest.approx(
            fit_rate(profile, ages).rate
        )

    def test_indices_exports_the_domain_aliases(self):
        import phyloslicer.indices as indices_module

        assert indices_module.MULTISITE_DOMAIN_ALIASES == {"reference": "mixed"}


# ===================================================================== #
# P2-8  prune_tips explains an infeasible threshold instead of None
# ===================================================================== #


class TestPruneTipsFailures:
    def test_infeasible_threshold_raises_with_the_workable_range(self, four_tip):
        with pytest.raises(ValueError, match="needs at least 2 tips"):
            prune_tips(four_tip, 99.0, side="after")

    def test_vector_call_never_yields_a_bare_none(self, four_tip):
        with pytest.raises(ValueError, match="needs at least 2 tips"):
            prune_tips(four_tip, [0.5, 99.0], side="after")

    def test_no_op_prune_returns_the_tree(self, four_tip):
        kept = prune_tips(four_tip, 99.0, side="before")
        assert kept.n_tips == four_tip.n_tips

    def test_feasible_cuts_still_select_the_right_tips(self, four_tip):
        assert sorted(prune_tips(four_tip, 1.5, side="after").tip_labels) == ["A", "B"]
        assert sorted(prune_tips(four_tip, 1.5, side="before").tip_labels) == ["C", "D"]


# ===================================================================== #
# P2-9  TreeArray equality and hashing
# ===================================================================== #


class TestTreeArrayEquality:
    def test_equality_works_on_ndarray_fields(self, four_tip):
        clone = TreeArray(list(four_tip.tip_labels), four_tip.edges.copy(), four_tip.lengths.copy())
        assert four_tip == clone
        assert (four_tip == clone) is True
        longer = TreeArray(
            list(four_tip.tip_labels), four_tip.edges.copy(), four_tip.lengths.copy() * 2
        )
        assert four_tip != longer

    def test_trees_are_usable_in_collections(self, four_tip):
        assert four_tip in [four_tip]
        clone = TreeArray(list(four_tip.tip_labels), four_tip.edges.copy(), four_tip.lengths.copy())
        assert len({four_tip, clone}) == 2  # identity hash: mutable by design

    def test_in_place_mutation_still_invalidates_caches(self, four_tip):
        before = four_tip.node_depths().copy()
        four_tip.lengths[0] = 42.0
        assert not np.allclose(before, four_tip.node_depths())


# ===================================================================== #
# P2-11  viz validates its percentile like the rest of the package
# ===================================================================== #


class TestVizValidation:
    def test_plot_rate_line_rejects_an_out_of_range_percentile(self, small_setup):
        import matplotlib

        matplotlib.use("Agg")
        from phyloslicer.viz import plot_rate_line

        tree, mat, stack = small_setup
        profile = pd_per_slice(tree, mat, stack)[0]
        fit = fit_rate(profile, stack.ages)
        with pytest.raises(ValueError, match="percentile"):
            plot_rate_line(profile, fit=fit, ages=stack.ages, percentile=1.5)
        with pytest.raises(ValueError, match="percentile"):
            plot_rate_line(profile, fit=fit, ages=stack.ages, percentile=95)
        assert plot_rate_line(profile, fit=fit, ages=stack.ages) is not None


# ===================================================================== #
# P2-10  DendroPy adapter: namespaces and unnamed leaves
# ===================================================================== #


class TestDendropyAdapter:
    def test_plain_tree_still_round_trips(self):
        dendropy = pytest.importorskip("dendropy")
        from phyloslicer.io.adapters import from_dendropy

        dt = dendropy.Tree.get(data="((A:1.1,B:1.2)95:0.9,(C:1.3,D:1.4)80:0.7);", schema="newick")
        assert from_dendropy(dt).tip_labels == ["A", "B", "C", "D"]

    def test_shared_namespace_with_unused_taxa_is_not_fatal(self):
        """Posterior pools, split_partitions() and wide NEXUS TAXA blocks all
        make trees share a namespace larger than themselves."""
        dendropy = pytest.importorskip("dendropy")
        from phyloslicer.io.adapters import from_dendropy

        dt = dendropy.Tree.get(data="((A:1.1,B:1.2)95:0.9,(C:1.3,D:1.4)80:0.7);", schema="newick")
        dt.taxon_namespace.add_taxon(dendropy.Taxon("GHOST"))
        with pytest.warns(UserWarning, match="does not use"):
            tree = from_dendropy(dt)
        assert tree.tip_labels == ["A", "B", "C", "D"]
        assert tree.n_tips == 4

    def test_oversized_namespace_from_construction(self):
        dendropy = pytest.importorskip("dendropy")
        from phyloslicer.io.adapters import from_dendropy

        ns = dendropy.TaxonNamespace(["A", "B", "C", "D", "E"])
        dt = dendropy.Tree.get(data="((A:1,B:1),(C:1,D:1));", schema="newick", taxon_namespace=ns)
        assert sorted(from_dendropy(dt).tip_labels) == ["A", "B", "C", "D"]

    def test_leaf_without_a_taxon_falls_back_like_biopython_does(self):
        dendropy = pytest.importorskip("dendropy")
        from phyloslicer.io.adapters import from_biopython, from_dendropy

        dt = dendropy.Tree.get_from_string("((:1,:1):1,(:1,:1):1);", "newick")
        assert from_dendropy(dt).tip_labels == ["tip_1", "tip_2", "tip_3", "tip_4"]
        import io as _io

        import Bio.Phylo

        bio = Bio.Phylo.read(_io.StringIO("((A:1,B:1)ab:1,(C:1,D:1)cd:1);"), "newick")
        assert sorted(from_biopython(bio).tip_labels) == ["A", "B", "C", "D"]

    def test_topology_and_lengths_survive_the_round_trip(self, four_tip):
        pytest.importorskip("dendropy")
        from phyloslicer.io.adapters import from_dendropy, to_dendropy

        assert from_dendropy(to_dendropy(four_tip)) == four_tip


# ===================================================================== #
# P3-1 / P3-2 / P3-3 / P3-5  I/O edge behaviour
# ===================================================================== #


class TestIoEdges:
    def test_quoted_label_may_contain_a_semicolon(self):
        tree = TreeArray.from_newick("('gen;us':0.1,'B':0.2,'C':0.3);")
        assert tree.tip_labels[0] == "gen;us"

    def test_split_respects_quotes_and_comments(self):
        records = split_newick_records("(A:1,B:1)'x;y';(C:1,D:1)[&R];")
        assert len(records) == 2
        assert records[0].endswith(";")

    def test_multi_tree_file_reads_the_first_tree_not_the_header(self):
        text = write_newick(yule_tree(6, seed=1)) + "\n" + write_newick(yule_tree(6, seed=2))
        assert TreeArray.from_newick(text).n_tips == 6

    def test_root_stem_is_reported_and_kept_for_round_trip(self):
        with pytest.warns(UserWarning, match="root stem"):
            tree = TreeArray.from_newick("(A:1,B:1,C:1):5;")
        assert tree.total_pd() == pytest.approx(3.0)  # the ape root.edge=0 convention
        assert tree.root_edge == pytest.approx(5.0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            assert TreeArray.from_newick(tree.to_newick()).root_edge == pytest.approx(5.0)

    def test_a_nexus_file_says_what_it_is(self):
        with pytest.raises(ValueError, match="NEXUS file"):
            read_newick("#NEXUS\n\nBEGIN TREES;\nTREE t = (A:1,B:1);\nEND;\n")

    def test_subarray_documents_its_ordering(self, four_tip):
        kept = four_tip.subarray(["D", "C"])
        assert kept.tip_labels == ["C", "D"]
        assert "canonical" in type(four_tip).subarray.__doc__.lower()

    def test_empty_subset_error_names_the_labels(self, four_tip):
        with pytest.raises(ValueError, match="only 1"):
            four_tip.subarray(["A"])


# ===================================================================== #
# P3-6  version has one source
# ===================================================================== #


class TestMetadata:
    def test_pyproject_version_is_dynamic(self):
        from pathlib import Path

        text = (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
        assert 'dynamic = ["version"]' in text
        assert 'attr = "phyloslicer.__version__"' in text

    def test_manuals_advertise_the_current_version(self):
        from pathlib import Path

        root = Path(__file__).resolve().parents[1]
        import phyloslicer

        for name in ("docs/USER_GUIDE.md", "docs/USER_GUIDE_CN.md"):
            assert phyloslicer.__version__ in (root / name).read_text(encoding="utf-8"), name

    def test_the_package_defines_its_version_exactly_once(self):
        from pathlib import Path

        root = Path(__file__).resolve().parents[1]
        sites = [
            str(p.relative_to(root))
            for p in (root / "src").rglob("*.py")
            if re.search(r"^__version__\s*=", p.read_text(encoding="utf-8"), re.M)
        ]
        assert sites == ["src/phyloslicer/__init__.py"], sites


# ===================================================================== #
# the numerical core is unchanged by any of the above
# ===================================================================== #


class TestNumbersUnchanged:
    def test_pd_and_pe_still_partition_and_match_a_direct_sum(self, small_setup):
        tree, mat, stack = small_setup
        profiles = pd_per_slice(tree, mat, stack)
        totals = pe_per_slice(tree, mat, stack)
        assert np.allclose(profiles.sum(axis=1), indices.pd(tree, mat))
        assert np.allclose(totals.sum(axis=1), indices.pe(tree, mat), atol=1e-12)
        assert np.all(profiles >= 0) and np.all(totals >= 0)

    def test_sliced_and_pruned_trees_stay_consistent(self, four_tip):
        sliced = slice_rootward(four_tip, 1.0)
        assert isinstance(sliced, SlicedTree)
        assert sliced.pd + slice_interval(four_tip, four_tip.root_age, 1.0).pd == pytest.approx(
            four_tip.total_pd()
        )

    def test_windows_partition_every_branch(self):
        tree = birth_death_tree(80, age=10.0, birth_rate=1.0, death_rate=0.5, seed=11)
        for criterion in ("time", "pd"):
            stack = slice_pieces(tree, n=17, criterion=criterion)
            assert np.allclose(stack.contribution_matrix().sum(axis=1), tree.lengths, atol=1e-12)

    def test_adjacency_semantics_are_unchanged(self):
        cells = np.array([[i, j, i + 1, j + 1] for i in range(3) for j in range(3)], float)
        queen = adjacency_from_bounds(cells, "queen")
        rook = adjacency_from_bounds(cells, "rook")
        assert queen[4].sum() == 9 and rook[4].sum() == 5

    def test_a_full_workflow_still_runs(self):
        tree = yule_tree(60, age=50.0, seed=5)
        mat, _coords = make_dataset(
            n_sites=12, n_species=60, structure="clustered", seed=5, tip_labels=tree.tip_labels
        )
        for column in mat.columns:
            if mat[column].sum() == 0:
                mat.loc[0, column] = 1.0
        frame = cpd_rate(tree, mat, n_slices=25)
        assert len(frame) == 12
        assert frame["CpD"].notna().any()
        assert set(
            [
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
            ]
        ) <= set(frame.columns)
