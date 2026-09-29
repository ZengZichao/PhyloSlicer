"""I/O, spatial, viz and CLI smoke tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phyloslicer.core import TreeArray
from phyloslicer.io import (
    align_tree_matrix,
    parse_newick,
    presence_matrix,
    read_newick,
    read_nexus,
    write_newick,
    write_nexus,
)
from phyloslicer.simulate import yule_tree
from phyloslicer.spatial import (
    adjacency_from_bounds,
    adjacency_from_coords,
    adjacency_from_grid,
    from_adjacency_dict,
    rate_map,
    resolve_adjacency_tol,
    to_adjacency_dict,
)


class TestNewick:
    def test_roundtrip(self, random_tree_20):
        text = write_newick(random_tree_20)
        back = parse_newick(text)
        # tip order follows newick traversal and may differ from the original
        assert sorted(back.tip_labels) == sorted(random_tree_20.tip_labels)
        assert np.allclose(sorted(back.lengths), sorted(random_tree_20.lengths))
        assert back.root_age == pytest.approx(random_tree_20.root_age)
        assert back.total_pd() == pytest.approx(random_tree_20.total_pd())

    def test_parse_with_support_values_and_internal_labels(self):
        tree = parse_newick("((A:1.0,B:1.0)95:1.0,(C:1.0,D:1.0)100:1.0);")
        assert tree.n_tips == 4
        assert tree.root_age == pytest.approx(2.0)
        assert tree.tip_labels == ["A", "B", "C", "D"]

    def test_quoted_names(self):
        tree = parse_newick("(('sp one':1,'sp two':1):1,'sp three':1);")
        assert "sp one" in tree.tip_labels

    def test_file_roundtrip(self, tmp_path, random_tree_20):
        p = tmp_path / "t.tre"
        p.write_text(write_newick(random_tree_20) + "\n")
        back = read_newick(p)
        assert back.n_tips == random_tree_20.n_tips

    def test_nexus_multi_tree(self, tmp_path):
        trees = [yule_tree(6, seed=s) for s in (1, 2)]
        text = write_nexus(trees)
        p = tmp_path / "t.nex"
        p.write_text(text)
        back = read_nexus(p)
        assert len(back) == 2
        assert back[0].n_tips == 6

    def test_errors(self):
        with pytest.raises(ValueError):
            parse_newick("(A:1,B:1)")  # missing ;


class TestMatrixIO:
    def test_presence_matrix_long_format(self):
        long = pd.DataFrame(
            {
                "site": ["a", "a", "b"],
                "species": ["S1", "S2", "S1"],
                "count": [1, 2, 0],
            }
        )
        wide = presence_matrix(long, site_col="site", species_col="species", value_col="count")
        assert wide.shape == (2, 2)
        assert wide.loc["a", "S1"] == 1.0
        assert wide.loc["b", "S1"] == 0.0
        assert wide.loc["a", "S2"] == 1.0

    def test_align_reports_drops(self, random_tree_20):
        mat = pd.DataFrame(
            np.eye(3, random_tree_20.n_tips),
            columns=random_tree_20.tip_labels,
            index=["s1", "s2", "s3"],
        )
        mat["ghost_species"] = 1.0  # absent from tree
        tree2, mat2 = align_tree_matrix(random_tree_20, mat)
        assert "ghost_species" not in mat2.columns
        assert list(mat2.columns) == tree2.tip_labels


class TestAdapters:
    def test_dendropy_roundtrip(self, random_tree_20):
        _dendropy = pytest.importorskip("dendropy")
        dnd = random_tree_20.to_dendropy()
        back = TreeArray.from_dendropy(dnd)
        assert back.tip_labels == random_tree_20.tip_labels
        assert back.root_age == pytest.approx(random_tree_20.root_age)
        assert back.total_pd() == pytest.approx(random_tree_20.total_pd())

    def test_biopython_roundtrip(self, random_tree_20):
        _bio = pytest.importorskip("Bio")
        bp = random_tree_20.to_biopython()
        back = TreeArray.from_biopython(bp)
        assert set(back.tip_labels) == set(random_tree_20.tip_labels)
        assert back.total_pd() == pytest.approx(random_tree_20.total_pd())


class TestSpatial:
    def test_bounds_rook_vs_queen(self):
        # 2x2 grid of unit squares
        bounds = np.array(
            [
                [0, 0, 1, 1],
                [1, 0, 2, 1],
                [0, 1, 1, 2],
                [1, 1, 2, 2],
            ],
            dtype=float,
        )
        queen = adjacency_from_bounds(bounds, method="queen")
        rook = adjacency_from_bounds(bounds, method="rook")
        # diagonal cells touch only at a corner -> queen only
        assert queen[0, 3] == 1
        assert rook[0, 3] == 0
        # side neighbours are rook-adjacent
        assert rook[0, 1] == 1 and rook[0, 2] == 1
        assert (np.diag(queen) == 1).all()

    def test_knn(self):
        x = np.array([0.0, 1, 2, 10])
        y = np.zeros(4)
        adj = adjacency_from_coords(x, y, method="knn", k=1)
        assert adj[0, 1] == 1 and adj[1, 0] == 1
        assert adj[0, 3] == 0
        assert (np.diag(adj) == 1).all()

    def test_radius(self):
        x = np.array([0.0, 1.0, 5.0])
        adj = adjacency_from_coords(x, np.zeros(3), method="radius", r=2.0)
        assert adj[0, 1] == 1 and adj[0, 2] == 0

    def test_adjacency_dict_roundtrip(self):
        adj = np.array([[1, 1, 0], [1, 1, 1], [0, 1, 1]])
        d = to_adjacency_dict(adj)
        back = from_adjacency_dict(d)
        assert (back == adj).all()

    def test_rate_map_renders(self, tmp_path):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        bounds = pd.DataFrame(
            {
                "minx": [0.0, 1.0],
                "miny": [0.0, 0.0],
                "maxx": [1.0, 2.0],
                "maxy": [1.0, 1.0],
            }
        )
        rates = pd.DataFrame({"CpD": [0.4, np.nan]})
        ax = rate_map(rates, bounds, rate="CpD")
        p = tmp_path / "map.png"
        ax.figure.savefig(p)
        plt.close(ax.figure)
        assert p.exists()


class TestCli:
    def test_slice_command(self, tmp_path, random_tree_20):
        from phyloslicer.cli.main import main
        from phyloslicer.slicing.windows import contribution_matrix

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(write_newick(random_tree_20) + "\n")
        out_path = tmp_path / "out.tre"
        log = tmp_path / "log.json"
        rc = main(
            [
                "slice",
                str(tree_path),
                "--rootward",
                "1.0",
                "--out",
                str(out_path),
                "--log",
                str(log),
            ]
        )
        assert rc == 0
        assert out_path.exists()
        sliced = read_newick(out_path)
        expected = contribution_matrix(
            random_tree_20,
            np.array([[random_tree_20.root_age - 1.0, random_tree_20.root_age]]),
        )[:, 0].sum()
        assert sliced.total_pd() == pytest.approx(expected)
        assert log.exists()

    def test_rates_and_validate_commands(self, tmp_path, random_tree_20):
        from phyloslicer.cli.main import main
        from phyloslicer.simulate import make_dataset

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(write_newick(random_tree_20) + "\n")
        mat, _ = make_dataset(n_sites=6, n_species=random_tree_20.n_tips, seed=1)
        mat.columns = random_tree_20.tip_labels
        sites_path = tmp_path / "sites.csv"
        mat.to_csv(sites_path)
        rates_path = tmp_path / "rates.csv"
        rc = main(
            ["rates", str(tree_path), str(sites_path), "--n-slices", "10", "--out", str(rates_path)]
        )
        assert rc == 0
        df = pd.read_csv(rates_path)
        assert {"CpD", "PD", "pDO", "converged"} <= set(df.columns)

        ref = pd.DataFrame({"CpD_ref": df["CpD"] * 1.0})
        ref_path = tmp_path / "ref.csv"
        ref.to_csv(ref_path)
        out = tmp_path / "val.csv"
        manifest = write_reference_manifest(tmp_path, ref_path, rate_column="CpD_ref")
        rc = main(
            [
                "validate",
                str(tree_path),
                str(sites_path),
                "--reference",
                str(ref_path),
                "--n-slices",
                "10",
                "--out",
                str(out),
                "--tolerance",
                "1e-8",
                "--manifest",
                str(manifest),
            ]
        )
        assert rc == 0
        # A13: the very same reference must be refused when it carries no manifest entry
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "10",
                    "--out",
                    str(tmp_path / "val2.csv"),
                    "--manifest",
                    str(tmp_path / "absent.txt"),
                ]
            )
            == 2
        )
        val = pd.read_csv(out)
        assert (val["equivalent"] == True).all()  # noqa: E712

    def test_version(self, capsys):
        from phyloslicer.cli.main import main

        assert main(["--version"]) == 0
        assert "phyloslicer" in capsys.readouterr().out


def write_reference_manifest(tmp_path, ref_path, rate_column="CpD_ref"):
    """Register an ad-hoc reference CSV through the package's own manifest writer (A13).

    `cli validate` refuses any reference that is not recorded in a provenance manifest, so a
    test that wants the success path has to declare provenance exactly as a real run does.
    """
    import hashlib
    import os

    import pandas as _pd

    from phyloslicer.cli.main import format_reference_manifest

    blob = open(ref_path, "rb").read()
    # the comparator reads references with the first column as the index;
    # record the shape the way it will see it, or A13 reports a false mismatch
    frame = _pd.read_csv(ref_path, index_col=0)
    text = format_reference_manifest(
        {"written": "test", "generator": "tests (synthetic fixture)", "reference_package": "none"},
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


class TestAdjacencyFromGrid:
    """``adjacency_from_grid`` was completely untested before this sweep."""

    @staticmethod
    def _frame(scale: float = 1.0) -> pd.DataFrame:
        """2x2 unit squares (0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)) scaled by ``scale``."""
        grid = pd.DataFrame(
            {
                "minx": [0.0, 1.0, 0.0, 1.0],
                "miny": [0.0, 0.0, 1.0, 1.0],
                "maxx": [1.0, 2.0, 1.0, 2.0],
                "maxy": [1.0, 1.0, 2.0, 2.0],
            }
        )
        return grid * scale

    def test_matches_the_bounds_reference_and_separates_rook_from_queen(self):
        grid = self._frame()
        queen = adjacency_from_grid(grid, method="queen")
        rook = adjacency_from_grid(grid, method="rook")
        assert np.array_equal(queen, adjacency_from_bounds(grid.to_numpy(float), method="queen"))
        assert np.array_equal(rook, adjacency_from_bounds(grid.to_numpy(float), method="rook"))
        assert queen[0, 3] == 1 and rook[0, 3] == 0  # corner contact is queen-only
        assert (np.diag(queen) == 1).all()
        # the diagonal is part of the contract (a focal neighbourhood holds the
        # focal site, treesliceR convention), so the counts include it: 4 self +
        # 12 queen links, 4 self + 8 rook links.
        assert queen.sum() == 16 and rook.sum() == 12
        assert (queen.sum(axis=1) == 4).all() and (rook.sum(axis=1) == 3).all()

    def test_contact_criterion_is_unit_invariant(self):
        """Degrees, metres and kilometres must give the same adjacency.

        The contact tolerance is scaled by the median cell edge, so a fixed
        absolute epsilon (the bug this guards) would glue cells together at
        1e-9 km while still separating them at 1e-9 degrees.
        """
        base_queen = adjacency_from_grid(self._frame(1.0), method="queen")
        base_rook = adjacency_from_grid(self._frame(1.0), method="rook")
        for scale in (1e-3, 1.0, 1e3, 1e6):
            assert np.array_equal(
                adjacency_from_grid(self._frame(scale), method="queen"), base_queen
            ), scale
            assert np.array_equal(
                adjacency_from_grid(self._frame(scale), method="rook"), base_rook
            ), scale

    def test_explicit_tol_is_absolute_by_design(self):
        """A5: the default tol scales with the cell; an explicit ``tol=`` never does.

        "Relative by default, absolute on request" cannot be shown on a unit grid
        alone (1e-9 x 1 == 1e-9), so the metre-scale copy of the same grid is
        judged against it: a gap that is a hair off contact in metres is a real
        gap in absolute terms, and the two calls must disagree there.
        """
        unit = self._frame()
        metric = self._frame(1e6)
        # the default: rel_tol x median cell edge
        assert resolve_adjacency_tol(unit.to_numpy(float)) == pytest.approx(1e-9)
        assert resolve_adjacency_tol(metric.to_numpy(float)) == pytest.approx(1e-3)
        # an explicit tol is returned verbatim, at any coordinate scale
        assert resolve_adjacency_tol(metric.to_numpy(float), tol=1e-9) == 1e-9
        assert resolve_adjacency_tol(unit.to_numpy(float), tol=1e-9) == 1e-9

        def gapped(scale: float, gap: float) -> pd.DataFrame:
            grid = self._frame(scale).copy()
            grid.loc[1, "minx"] += gap  # cell 1 pulled ``gap`` away from cell 0
            return grid

        # 1e-4 m on a 1e6 m cell: inside the relative default, outside tol=1e-9
        assert adjacency_from_grid(gapped(1e6, 1e-4), method="queen")[0, 1] == 1
        assert adjacency_from_grid(gapped(1e6, 1e-4), method="queen", tol=1e-9)[0, 1] == 0
        # the same *absolute* tol=1e-9 on the unit grid calls a 1e-12 gap (a
        # relative hair there) touching: raw coordinates, whatever the cell
        # measure, so the escape hatch stays absolute by design
        assert adjacency_from_grid(gapped(1.0, 1e-12), method="queen", tol=1e-9)[0, 1] == 1
        assert adjacency_from_grid(gapped(1.0, 1e-6), method="queen")[0, 1] == 0
        assert adjacency_from_grid(gapped(1.0, 1e-6), method="queen", tol=1e-5)[0, 1] == 1
        # a gap beyond every tolerance stays a gap on both grids
        assert adjacency_from_grid(gapped(1e6, 1.0), method="queen")[0, 1] == 0
        assert adjacency_from_grid(gapped(1.0, 1e-3), method="queen", tol=1e-9)[0, 1] == 0

    def test_input_validation(self):
        grid = self._frame()
        with pytest.raises(ValueError, match="minx"):
            adjacency_from_grid(pd.DataFrame({"a": [1.0, 2.0]}))
        with pytest.raises(ValueError, match="method"):
            adjacency_from_grid(grid, method="knn")

    def test_geopandas_polygon_semantics(self):
        """L-shaped cells: bounding boxes overlap, geometry may not."""
        gpd = pytest.importorskip("geopandas")
        from shapely.geometry import Polygon

        a = Polygon([(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)])
        b = Polygon([(2, 1), (3, 1), (3, 3), (2, 3)])  # touches a at the corner (2,1)
        gdf = gpd.GeoDataFrame(geometry=[a, b])
        queen = adjacency_from_grid(gdf, method="queen")
        rook = adjacency_from_grid(gdf, method="rook")
        assert queen[0, 1] == 1  # corner contact counts for queen
        assert rook[0, 1] == 0  # a point is not a shared boundary segment
        # ... while the bounds-only fallback would call them overlapping anyway
        bounds = pd.DataFrame(
            np.asarray([list(g.bounds) for g in (a, b)], dtype=float),
            columns=["minx", "miny", "maxx", "maxy"],
        )
        assert adjacency_from_bounds(bounds.to_numpy(float), method="rook")[0, 1] == 1


class TestViz:
    @pytest.fixture(autouse=True)
    def _agg(self):
        import matplotlib

        matplotlib.use("Agg")

    @staticmethod
    def _legend_texts(ax):
        """The strings the legend actually draws.

        ``Text.get_label()`` is the *artist* label (always ``''`` for a legend
        entry); only ``Text.get_text()`` returns the entry's string, so the
        naive accessor made every legend assertion below pass or fail on
        nothing at all.
        """
        legend = ax.get_legend()
        return [] if legend is None else [t.get_text() for t in legend.get_texts()]

    def test_plot_rate_line_uses_the_given_ages_and_labels_the_fit(self):
        from phyloslicer.rates.fitting import fit_rate
        from phyloslicer.viz import plot_rate_line

        ages = np.linspace(2.0, 0.2, 10)  # slice upper boundaries, root -> tip
        # a *positive* exponentially accumulating profile: the previous fixture
        # diffed the curve the other way round, so it summed to < 0, fit_rate
        # refused it (rate = nan) and plot_rate_line drew neither the fitted
        # curve nor the origin marker - the two legend assertions below then
        # checked an empty legend.
        profile = np.diff(np.concatenate([[0.0], np.exp(-1.3 * ages)]))
        assert (profile > 0).all()
        fit = fit_rate(profile, ages)
        assert fit.converged and np.isfinite(fit.rate) and fit.rate > 0
        ax = plot_rate_line(profile, fit=fit, ages=ages)
        assert ax.collections, "the observed cumulative points are missing"
        xy = ax.collections[0].get_offsets()
        # P1-7: the shipped default runs root -> present, so the x positions are
        # the caller's ages measured from the root, never a placeholder index.
        # Element-wise (not sorted): point i must pair with ages[i].
        assert np.allclose(xy[:, 0], ages.max() - ages)  # x axis really is age
        assert "from root" in ax.get_xlabel().lower()
        labels = self._legend_texts(ax)
        assert any(label.startswith("fit: r =") for label in labels), labels
        assert any(label.startswith("95% origin") for label in labels), labels
        assert ax.get_ylim()[1] == pytest.approx(1.02)
        # the origin marker must sit where the drawn curve crosses the
        # percentile *in the same coordinates as the points*
        x0 = -np.log(1 - 0.95) / fit.rate
        origin = next(line for line in ax.lines if line.get_label().startswith("95% origin"))
        assert origin.get_xdata()[0] == pytest.approx(ages.max() - x0)
        assert f"{x0:.3g} Ma" in origin.get_label()

        mirrored = plot_rate_line(profile, fit=fit, ages=ages, direction="from_present")
        assert np.allclose(mirrored.collections[0].get_offsets()[:, 0], ages)
        origin_back = next(
            line for line in mirrored.lines if line.get_label().startswith("95% origin")
        )
        assert origin_back.get_xdata()[0] == pytest.approx(x0)  # raw age again
        assert "age" in mirrored.get_xlabel().lower()

    def test_plot_rate_line_warns_about_placeholder_ages(self):
        from phyloslicer.viz import plot_rate_line

        with pytest.warns(UserWarning, match="placeholder spacing"):
            plot_rate_line(np.full(6, 1 / 6))

    def test_plot_sensitivity_is_a_log_scale_rate_plot(self):
        from phyloslicer.viz import plot_sensitivity

        df = pd.DataFrame(
            {
                "site_id": ["a"] * 4 + ["b"] * 4,
                "n_slices": [10, 25, 50, 100] * 2,
                "rate": [1.2, 1.1, 1.05, 1.02, 0.8, 0.79, 0.78, 0.77],
            }
        )
        ax = plot_sensitivity(df)
        assert ax.get_xscale() == "log"
        assert ax.get_xlabel() == "number of slices"
        assert len(ax.lines) == 2
        assert np.allclose(ax.lines[0].get_xdata(), [10, 25, 50, 100])

    def test_plot_posterior_labels_the_actual_hdi_level(self, span_x_extent):
        """The shading must be named for the level that produced it (A9-era bug)."""
        from phyloslicer.viz import plot_posterior

        frame = pd.DataFrame(
            {
                "rate_median": [1.2],
                "hdi_low": [1.0],
                "hdi_high": [1.4],
            }
        )
        frame.attrs["hdi_level"] = 0.68
        ax = plot_posterior(frame)
        labels = self._legend_texts(ax)
        assert "68% HDI" in labels, labels
        assert not any("95" in label for label in labels), labels
        assert span_x_extent(ax, "68% HDI") == (pytest.approx(1.0), pytest.approx(1.4))

    def test_plot_posterior_explicit_hdi_overrides_attrs(self):
        from phyloslicer.viz import plot_posterior

        frame = pd.DataFrame({"rate_median": [1.2], "hdi_low": [1.0], "hdi_high": [1.4]})
        frame.attrs["hdi_level"] = 0.95
        ax = plot_posterior(frame, hdi=0.9)
        labels = self._legend_texts(ax)
        assert "90% HDI" in labels, labels
        assert not any("95%" in label for label in labels), labels

    def test_plot_posterior_hdi_label_override_is_verbatim(self):
        """``hdi_label=`` lets the caller word the legend; the level still gates it."""
        from phyloslicer.viz import plot_posterior

        frame = pd.DataFrame({"rate_median": [1.2], "hdi_low": [1.0], "hdi_high": [1.4]})
        ax = plot_posterior(frame, hdi=0.95, hdi_label="95% credible")
        assert self._legend_texts(ax) == ["95% credible"]
        with pytest.raises(ValueError, match="hdi must lie"):
            plot_posterior(frame, hdi=1.5)

    def test_plot_posterior_refuses_to_guess_the_hdi_level(self):
        from phyloslicer.viz import plot_posterior

        frame = pd.DataFrame({"rate_median": [1.2], "hdi_low": [1.0], "hdi_high": [1.4]})
        with pytest.raises(ValueError, match="cannot label the HDI"):
            plot_posterior(frame)

    def test_plot_posterior_validates_inputs(self):
        from phyloslicer.viz import plot_posterior

        rng = np.random.default_rng(0)
        frame = pd.DataFrame({"rate_median": rng.random(12)})
        with pytest.raises(ValueError, match="kind must be"):
            plot_posterior(frame, kind="density")
        with pytest.raises(ValueError, match="no finite values"):
            plot_posterior(pd.DataFrame({"rate_median": [np.nan, np.nan]}))

    def test_posterior_multi_site_frame_says_so(self):
        from phyloslicer.viz import plot_posterior

        frame = pd.DataFrame(
            {
                "rate_median": [1.0, 2.0],
                "hdi_low": [0.5, 1.5],
                "hdi_high": [1.5, 2.5],
            }
        )
        ax = plot_posterior(frame, hdi=0.9)
        assert "no single 90% HDI" in ax.get_title()


class TestCliPosterior:
    def test_posterior_command_reports_hdi_and_tree_counts(self, tmp_path):
        from phyloslicer.cli.main import main

        trees_dir = tmp_path / "trees"
        trees_dir.mkdir()
        for seed in (1, 2, 3):
            tree = yule_tree(6, age=1.0, seed=seed)
            (trees_dir / f"t{seed}.tre").write_text(write_newick(tree) + "\n")
        labels = sorted(yule_tree(6, age=1.0, seed=1).tip_labels)
        rng = np.random.default_rng(3)
        mat = pd.DataFrame(
            (rng.random((4, len(labels))) < 0.7).astype(float),
            columns=labels,
            index=["s1", "s2", "s3", "s4"],
        )
        for col in mat.columns:  # every species present at least once
            if mat[col].sum() == 0:
                mat.loc["s1", col] = 1.0
        sites_path = tmp_path / "sites.csv"
        mat.to_csv(sites_path)
        out = tmp_path / "post.csv"
        assert (
            main(
                [
                    "posterior",
                    str(trees_dir),
                    str(sites_path),
                    "--n-slices",
                    "5",
                    "--hdi",
                    "0.68",
                    "--out",
                    str(out),
                ]
            )
            == 0
        )
        res = pd.read_csv(out)
        assert len(res) == 4
        assert {"rate_mean", "rate_median", "hdi_low", "hdi_high", "n_trees", "n_converged"} <= set(
            res.columns
        )
        assert (res["n_trees"] == 3).all()
        assert (res["hdi_low"] <= res["rate_median"]).all()
        assert (res["rate_median"] <= res["hdi_high"]).all()
        assert np.isfinite(res["rate_median"]).all()
        # an explicit hdi level really changes the interval (it is not a fixed 95%)
        res95 = tmp_path / "post95.csv"
        main(
            [
                "posterior",
                str(trees_dir),
                str(sites_path),
                "--n-slices",
                "5",
                "--hdi",
                "0.99",
                "--out",
                str(res95),
            ]
        )
        wide = pd.read_csv(res95)
        assert ((wide["hdi_high"] - wide["hdi_low"]) >= (res["hdi_high"] - res["hdi_low"])).all()

    def test_posterior_command_rejects_bad_hdi(self, tmp_path):
        from phyloslicer.cli.main import main

        tree = yule_tree(6, age=1.0, seed=1)
        (tmp_path / "trees.tre").write_text(write_newick(tree) + "\n")
        mat = pd.DataFrame(np.ones((2, 6)), columns=tree.tip_labels, index=["a", "b"])
        sites = tmp_path / "sites.csv"
        mat.to_csv(sites)
        assert (
            main(
                [
                    "posterior",
                    str(tmp_path / "trees.tre"),
                    str(sites),
                    "--hdi",
                    "1.5",
                    "--n-slices",
                    "3",
                    "--out",
                    str(tmp_path / "x.csv"),
                ]
            )
            == 1
        )  # CLI boundary: reported as an error, not a traceback
