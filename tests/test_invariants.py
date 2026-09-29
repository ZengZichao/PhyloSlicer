"""Property tests for the invariants the reference implementation never asserted.

Every test below is written against the signatures the modules actually expose
(checked in ``slicing/{api,windows}.py``, ``indices/api.py``, ``spatial``,
``viz``, ``uncertainty``, ``cli/main.py``, ``io/matrices.py``):

* ``contribution_matrix(tree, windows)`` takes an **(k, 2) depth-window array**,
  not a :class:`SliceStack` - the stack's own ``contribution_matrix()`` method
  forwards ``stack.depth_windows`` to it (pinned by
  :func:`test_contribution_matrix_accepts_a_window_array_not_a_stack`);
* ``pd_per_slice`` / ``pe_per_slice`` return plain ``(n_sites, k)`` ndarrays,
  ``rate_profile`` returns a ``dict`` of arrays, ``pb`` returns a ``DataFrame``
  with ``site_id`` / ``PB`` columns, and ``adjacency_from_grid`` takes a
  **frame** (bounds columns, or a GeoDataFrame), not ``rows=``/``cols=``;
* ``pb`` / ``pb_per_slice`` take a dense 0/1 adjacency whose *rows* are the
  focal neighbourhoods; the diagonal marks a site as a member of its own
  neighbourhood, which is the treesliceR convention tested below.

The randomised cases use seeded ``numpy`` generators over Yule and birth-death
trees rather than ``hypothesis``: the dev extra was declared and never used,
and randomised property checking of exactly that kind needs no
dependency the test job cannot import.

Each assertion is meant to *bite* - a bound that is approached but never
crossed, window boundaries that must line up, a legend that must carry the
caller's HDI level, a pruned edge set that must equal an independently
reconstructed one.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

from phyloslicer import posterior_rate
from phyloslicer.indices import pb as pb_total
from phyloslicer.indices import pb_per_slice, pd_per_slice, pe_per_slice
from phyloslicer.indices import pd as pd_total
from phyloslicer.indices import pe as pe_total
from phyloslicer.io import align_tree_matrix, write_newick
from phyloslicer.rates import rate_profile
from phyloslicer.simulate import birth_death_tree, yule_tree
from phyloslicer.slicing import SliceStack, slice_pieces
from phyloslicer.slicing.windows import contribution_matrix
from phyloslicer.spatial import adjacency_from_bounds, adjacency_from_grid

# ------------------------------------------------------------------ fixtures
# (name, tree builder, seed) - explicit seeds: ``hash(str)`` is salted per
# process and would make these "randomised" cases irreproducible.
TREES = [
    ("yule-13", lambda seed: yule_tree(13, age=4.0, seed=seed), 101),
    ("yule-24", lambda seed: yule_tree(24, age=7.5, seed=seed), 202),
    ("bd-18", lambda seed: birth_death_tree(18, age=3.0, seed=seed), 303),
]


def _community(tree, seed: int, density: float = 0.45) -> pd.DataFrame:
    """Random presence/absence frame with no empty site and no empty species.

    ``align_tree_matrix`` drops all-absent species and refuses an all-absent
    tree/matrix pair, so the randomisation is clipped to stay inside the
    contract the index kernels are documented for.
    """
    rng = np.random.default_rng(seed)
    values = (rng.random((6, tree.n_tips)) < density).astype(float)
    values[0, :] = 1.0  # one ubiquitous assemblage -> every tip occurs somewhere
    values[:, 0] = 1.0
    return pd.DataFrame(
        values, columns=tree.tip_labels, index=[f"site{i}" for i in range(values.shape[0])]
    )


def _aligned(builder, seed: int):
    tree = builder(seed)
    assert tree.is_ultrametric, "the randomiser must produce ultrametric trees"
    assert tree.root_age > 0.0
    return align_tree_matrix(tree, _community(tree, seed + 1))


def _neighbourhood(n_sites: int, seed: int) -> np.ndarray:
    """A random symmetric 0/1 adjacency with a full diagonal."""
    rng = np.random.default_rng(seed)
    adj = (rng.random((n_sites, n_sites)) < 0.5).astype(np.int8)
    adj = np.maximum(adj, adj.T)
    np.fill_diagonal(adj, 1)
    for i in range(n_sites - 1):  # every focal row really has >= 2 members
        adj[i, i + 1] = adj[i + 1, i] = 1
    return adj


# ------------------------------------------- 1. per-slice PD / PE range bounds
@pytest.mark.parametrize("name,builder,seed", TREES)
@pytest.mark.parametrize("n_slices", [4, 9])
def test_per_slice_pd_is_bounded_by_zero_and_by_the_tree(name, builder, seed, n_slices):
    """§5: no test ever asserted ``0 <= PD_slice <= whole-tree PD``."""
    tree, mat = _aligned(builder, seed)
    slices = slice_pieces(tree, n=n_slices)
    prof = pd_per_slice(tree, mat, slices)
    whole_tree_pd = tree.total_pd()
    per_site_pd = pd_total(tree, mat)  # (n_sites,) Faith's PD of each assemblage

    assert prof.shape == (mat.shape[0], n_slices)
    assert np.all(prof >= 0.0), f"negative per-slice PD: {prof.min()}"
    # a slice of a site cannot store more PD than the site, nor than the tree
    assert np.all(prof <= per_site_pd[:, None] + 1e-12)
    assert np.all(prof <= whole_tree_pd + 1e-12)
    # the bound is approached but never crossed: slices rebuild the site total
    np.testing.assert_allclose(prof.sum(axis=1), per_site_pd, rtol=1e-12, atol=1e-12)
    assert 0.0 < prof.max() < whole_tree_pd, "a slice holding the whole tree is not a slice"


@pytest.mark.parametrize("name,builder,seed", TREES)
@pytest.mark.parametrize("n_slices", [4, 9])
def test_per_slice_pe_is_bounded_by_zero_pd_and_the_tree(name, builder, seed, n_slices):
    tree, mat = _aligned(builder, seed)
    slices = slice_pieces(tree, n=n_slices)
    pe_prof = pe_per_slice(tree, mat, slices)
    pd_prof = pd_per_slice(tree, mat, slices)
    per_site_pe = pe_total(tree, mat)

    assert np.all(pe_prof >= 0.0), f"negative per-slice PE: {pe_prof.min()}"
    # range weighting divides each covered edge by the number of sites sharing
    # it, so per-slice endemism never exceeds the same slice's PD ...
    assert np.all(pe_prof <= pd_prof + 1e-12)
    # ... nor the site's total endemism, nor the whole tree
    assert np.all(pe_prof <= per_site_pe[:, None] + 1e-12)
    assert np.all(pe_prof <= tree.total_pd() + 1e-12)
    np.testing.assert_allclose(pe_prof.sum(axis=1), per_site_pe, rtol=1e-12, atol=1e-12)
    assert 0.0 < pe_prof.max() < tree.total_pd()
    assert np.all(per_site_pe <= pd_total(tree, mat) + 1e-12)


def test_rate_profile_returns_absolute_arrays_that_match_the_kernels():
    """``rate_profile`` is a dict of arrays, and ``normalize=False`` really is
    the absolute (treesliceR ``r_phylo``) view rather than a re-scaled one."""
    tree, mat = _aligned(TREES[0][1], TREES[0][2])
    slices = slice_pieces(tree, n=6)
    out = rate_profile(tree, mat, n_slices=6, index="PD", normalize=False)
    assert isinstance(out, dict)
    assert set(out) == {"values", "ages", "status", "ultrametric"}
    np.testing.assert_allclose(out["values"], pd_per_slice(tree, mat, slices))
    np.testing.assert_allclose(out["ages"], slices.ages)
    relative = rate_profile(tree, mat, n_slices=6, index="PD", normalize=True)["values"]
    np.testing.assert_allclose(relative.sum(axis=1), np.ones(len(relative)))
    assert not np.allclose(relative, out["values"])


# ------------------------------------------------ 2. beta index inside [0, 1]
@pytest.mark.parametrize("component", ["sorensen", "turnover", "nestedness"])
@pytest.mark.parametrize("approach", ["multisite", "pairwise"])
def test_beta_index_stays_inside_the_unit_interval(component, approach):
    """§5: nothing pinned ``0 <= PB <= 1`` for the six component x approach cells."""
    tree, mat = _aligned(TREES[1][1], TREES[1][2])
    adj = _neighbourhood(mat.shape[0], 901)
    domains = ["mixed", "paired"] if approach == "multisite" else ["mixed"]
    for domain in domains:
        frame = pb_total(
            tree, mat, adj, component=component, approach=approach, multisite_domain=domain
        )
        assert isinstance(frame, pd.DataFrame)
        # `status` carries why a neighbourhood was refused, so that pb()'s NaN
        # pattern is the same one pb_per_slice() reports for its `total`
        assert list(frame.columns) == ["site_id", "PB", "status"]
        assert set(frame["status"]) == {"ok"}, frame
        values = frame["PB"].to_numpy(dtype=float)
        finite = values[np.isfinite(values)]
        assert len(finite) == mat.shape[0], frame  # not a vacuous "all NaN" pass
        assert np.all(finite >= 0.0), (component, approach, domain, finite.min())
        assert np.all(finite <= 1.0), (component, approach, domain, finite.max())
        # a constant profile would satisfy the bounds trivially
        assert finite.max() - finite.min() > 1e-6, (component, approach, domain, finite)
        assert 0.0 < finite.min() and finite.max() < 1.0, (component, approach, finite)


def test_per_slice_beta_shares_and_the_whole_tree_index_agree():
    tree, mat = _aligned(TREES[1][1], TREES[1][2] + 10)
    adj = _neighbourhood(mat.shape[0], 911)
    slices = slice_pieces(tree, n=7)
    result = pb_per_slice(tree, mat, adj, slices, component="sorensen")
    values = result["values"]
    usable = np.isfinite(values).all(axis=1)
    assert usable.sum() == mat.shape[0]
    np.testing.assert_allclose(values[usable].sum(axis=1), 1.0, rtol=1e-12)
    np.testing.assert_allclose(
        result["total"][usable],
        pb_total(tree, mat, adj)["PB"].to_numpy(dtype=float)[usable],
        rtol=1e-12,
    )


# ------------------------------------------ 3. windows tile [0, T], C rebuilds L
@pytest.mark.parametrize("criterion", ["time", "pd"])
@pytest.mark.parametrize("n_slices", [1, 2, 5, 11])
@pytest.mark.parametrize("name,builder,seed", TREES)
def test_windows_partition_the_root_age_interval(name, builder, seed, criterion, n_slices):
    tree = builder(seed)
    stack = slice_pieces(tree, n=n_slices, criterion=criterion)
    windows = np.asarray(stack.depth_windows, dtype=float)

    assert windows.shape == (n_slices, 2)
    assert windows[0, 0] == pytest.approx(0.0, abs=1e-12)  # starts at the root
    assert windows[-1, 1] == pytest.approx(tree.root_age)  # ends at the tips
    assert np.all(windows[:, 1] > windows[:, 0]), "empty or reversed window"
    # consecutive windows share their boundary exactly: no gap, no overlap
    np.testing.assert_allclose(windows[:-1, 1], windows[1:, 0], rtol=1e-12, atol=1e-12)
    assert np.all(np.diff(windows[:, 0]) > 0)
    # the Ma-before-present view is the same tiling, flipped and re-anchored
    np.testing.assert_allclose(stack.windows[:, 0] - stack.windows[:, 1], stack.widths)
    assert stack.windows[0, 0] == pytest.approx(tree.root_age)
    assert stack.windows[-1, 1] == pytest.approx(0.0, abs=1e-12)
    if criterion == "time":
        np.testing.assert_allclose(stack.widths, tree.root_age / n_slices, rtol=1e-12)


@pytest.mark.parametrize("criterion", ["time", "pd"])
@pytest.mark.parametrize("name,builder,seed", TREES)
def test_contribution_matrix_rows_reconstruct_every_edge_length(name, builder, seed, criterion):
    """``Sum_j C[e, j] == L_e`` for every edge - the partition identity."""
    tree = builder(seed)
    stack = slice_pieces(tree, n=6, criterion=criterion)
    C = contribution_matrix(tree, stack.depth_windows)

    assert C.shape == (tree.n_edges, stack.n_slices)
    assert np.all(C >= 0.0)
    np.testing.assert_allclose(C.sum(axis=1), tree.lengths, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(C.sum(), tree.total_pd(), rtol=1e-12)
    np.testing.assert_allclose(stack.contribution(), C.sum(axis=0))
    np.testing.assert_allclose(stack.contribution_matrix(), C)
    if criterion == "pd":  # equal-PD windows store equal amounts of branch length
        np.testing.assert_allclose(C.sum(axis=0), tree.total_pd() / 6, rtol=1e-9)
    else:  # equal-time windows: equal widths, but unequal PD
        np.testing.assert_allclose(stack.widths, tree.root_age / 6, rtol=1e-12)
        assert C.sum(axis=0).max() > 1.5 * C.sum(axis=0).min()


def test_contribution_matrix_accepts_a_window_array_not_a_stack():
    """Pin the signature: an (k, 2) array of depth windows, never the stack."""
    tree, _mat = _aligned(TREES[0][1], TREES[0][2] + 20)
    stack = slice_pieces(tree, n=3)
    direct = contribution_matrix(tree, stack.depth_windows)
    np.testing.assert_allclose(direct, stack.contribution_matrix())
    # a stack is not a window array: the kernel must refuse it rather than
    # silently reinterpreting it; this test pins which of the two it takes
    with pytest.raises((TypeError, ValueError), match=r"float|\(k, 2\)"):
        contribution_matrix(tree, stack)
    with pytest.raises(ValueError, match="hi must be"):
        contribution_matrix(tree, np.array([[1.0, 0.2]]))


# --------------------------------------------- 4. spatial.adjacency_from_grid
# a 2 x 3 grid of unit squares, row-major: 0=(r0,c0) 1=(r0,c1) 2=(r0,c2)
#                                          3=(r1,c0) 4=(r1,c1) 5=(r1,c2)
QUEEN_2X3 = [
    [1, 1, 0, 1, 1, 0],
    [1, 1, 1, 1, 1, 1],
    [0, 1, 1, 0, 1, 1],
    [1, 1, 0, 1, 1, 0],
    [1, 1, 1, 1, 1, 1],
    [0, 1, 1, 0, 1, 1],
]
ROOK_2X3 = [
    [1, 1, 0, 1, 0, 0],
    [1, 1, 1, 0, 1, 0],
    [0, 1, 1, 0, 0, 1],
    [1, 0, 0, 1, 1, 0],
    [0, 1, 0, 1, 1, 1],
    [0, 0, 1, 0, 1, 1],
]


def _grid_frame(scale: float = 1.0) -> pd.DataFrame:
    rows, cols = 2, 3
    cell_r = np.repeat(np.arange(rows), cols).astype(float)
    cell_c = np.tile(np.arange(cols), rows).astype(float)
    return pd.DataFrame(
        {
            "minx": cell_c * scale,
            "miny": cell_r * scale,
            "maxx": (cell_c + 1) * scale,
            "maxy": (cell_r + 1) * scale,
        }
    )


def test_adjacency_from_grid_separates_rook_from_queen():
    """§5: ``spatial.adjacency_from_grid`` had no test at all."""
    grid = _grid_frame()
    queen = adjacency_from_grid(grid, method="queen")
    rook = adjacency_from_grid(grid, method="rook")

    assert queen.shape == rook.shape == (6, 6)
    np.testing.assert_array_equal(queen, np.array(QUEEN_2X3, dtype=queen.dtype))
    np.testing.assert_array_equal(rook, np.array(ROOK_2X3, dtype=rook.dtype))
    # corner contacts are queen-only: cells 1=(r0,c1) and 3=(r1,c0) meet at the
    # corner of their 2 x 2 block, cells 2=(r0,c2) and 3=(r1,c0) do not touch
    assert queen[1, 3] == 1 and rook[1, 3] == 0
    assert queen[0, 4] == 1 and rook[0, 4] == 0
    assert queen[2, 3] == 0 and rook[2, 3] == 0
    # queen strictly refines rook, and both are symmetric
    assert np.all((rook == 1) <= (queen == 1))
    np.testing.assert_array_equal(queen, queen.T)
    np.testing.assert_array_equal(rook, rook.T)
    # 11 corner-or-edge contacts plus 6 self-links / 7 shared edges plus 6
    assert queen.sum() == 2 * 11 + 6
    assert rook.sum() == 2 * 7 + 6


def test_adjacency_from_grid_treats_every_site_as_its_own_neighbour():
    """The reflexive diagonal is the treesliceR convention the indices assume."""
    for method in ("rook", "queen"):
        grid = _grid_frame(0.5)
        adj = adjacency_from_grid(grid, method=method)
        assert np.all(np.diag(adj) == 1), f"{method}: a focal row holds its own site"
        assert (adj.sum(axis=1) >= 2).all()
        # ... and the grid path is the bounds path for rectangular cells
        bounds = grid[["minx", "miny", "maxx", "maxy"]].to_numpy(dtype=float)
        np.testing.assert_array_equal(
            adj, adjacency_from_bounds(bounds, method=method), err_msg=method
        )


def test_adjacency_from_grid_rejects_an_unknown_method_and_a_bad_frame():
    with pytest.raises(ValueError, match="method must be"):
        adjacency_from_grid(_grid_frame(), method="knight")
    with pytest.raises(ValueError, match="minx"):
        adjacency_from_grid(pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}))


# ------------------------------------------------------------------- 5. viz
def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _legend_texts(ax):
    legend = ax.get_legend()
    return [] if legend is None else [t.get_text() for t in legend.get_texts()]


def _posterior_frame(hdi: float):
    trees = [yule_tree(12, age=2.0, seed=s) for s in (1, 2, 3, 4, 5)]
    mat = pd.DataFrame(
        np.eye(3, 12) + np.eye(3, 12, k=1),
        columns=trees[0].tip_labels,
        index=["s0", "s1", "s2"],
    )
    post = posterior_rate(trees, mat, n_slices=8, hdi=hdi)
    assert post.attrs["hdi_level"] == pytest.approx(hdi)
    return post


def test_plot_rate_line_scales_the_x_axis_by_the_ages_it_was_given():
    plt = _pyplot()
    from phyloslicer.rates.fitting import fit_rate
    from phyloslicer.viz import plot_rate_line

    try:
        ages = np.linspace(2.0, 0.2, 10)  # root -> tip, like SliceStack.ages
        profile = np.diff(np.concatenate([[0.0], np.exp(-1.3 * ages)]))
        fit = fit_rate(profile, ages)
        ax = plot_rate_line(profile, fit=fit, ages=ages)
        xy = ax.collections[0].get_offsets()
        # the default direction is "from_root": x counts up from the root
        np.testing.assert_allclose(np.sort(xy[:, 0]), np.sort(ages.max() - ages))
        np.testing.assert_allclose(xy[:, 1], np.cumsum(profile) / profile.sum())
        assert xy[0, 0] == pytest.approx(0.0)
        assert xy[-1, 1] == pytest.approx(1.0)
        labels = _legend_texts(ax)
        assert any(text.startswith("fit: r =") for text in labels), labels
        assert any(text.startswith("95% origin") for text in labels), labels
        assert ax.get_ylim()[1] == pytest.approx(1.02)
        ax2 = plot_rate_line(profile, fit=fit, ages=ages, direction="from_present")
        np.testing.assert_allclose(np.sort(ax2.collections[0].get_offsets()[:, 0]), np.sort(ages))
        assert ax2.get_xlabel() == "age (Ma before present)"
    finally:
        plt.close("all")


def test_plot_sensitivity_draws_one_rate_trace_per_site():
    plt = _pyplot()
    from phyloslicer.viz import plot_sensitivity

    try:
        df = pd.DataFrame(
            {
                "site_id": ["a"] * 3 + ["b"] * 3,
                "n_slices": [4, 8, 16] * 2,
                "rate": [1.1, 1.05, 1.0, 0.8, 0.78, 0.77],
            }
        )
        ax = plot_sensitivity(df)
        assert ax.get_xscale() == "log"
        assert ax.get_xlabel() == "number of slices"
        assert len(ax.lines) == 2
        np.testing.assert_allclose(ax.lines[0].get_xdata(), [4, 8, 16])
        np.testing.assert_allclose(ax.lines[0].get_ydata(), [1.1, 1.05, 1.0])
        np.testing.assert_allclose(ax.lines[1].get_ydata(), [0.8, 0.78, 0.77])
    finally:
        plt.close("all")


def test_plot_posterior_labels_the_hdi_level_it_was_given(span_x_extent):
    """viz used to annotate every panel "95% HDI", whatever the caller asked."""
    plt = _pyplot()
    from phyloslicer.viz import plot_posterior

    try:
        single = _posterior_frame(0.68).iloc[[0]].reset_index(drop=True)
        ax = plot_posterior(single)
        labels = _legend_texts(ax)
        assert any("68% HDI" in text for text in labels), labels
        assert not any("95" in text for text in labels), labels
        # the shaded span really is the interval the frame carries
        assert span_x_extent(ax, "68% HDI") == (
            pytest.approx(single["hdi_low"].iloc[0]),
            pytest.approx(single["hdi_high"].iloc[0]),
        )
        # an explicit level beats the recorded one, in either direction
        frame = _posterior_frame(0.95)
        for hdi, wanted in ((0.9, "90% HDI"), (0.5, "50% HDI"), (0.68, "68% HDI")):
            one = frame.iloc[[0]].reset_index(drop=True)
            ax_hdi = plot_posterior(one, hdi=hdi)
            texts = _legend_texts(ax_hdi)
            assert any(wanted in text for text in texts), (hdi, texts)
            assert not any("95%" in text for text in texts), (hdi, texts)
            plt.close("all")
        # a frame with no recorded level must fail, not guess
        bare = pd.DataFrame({"rate_median": [1.2], "hdi_low": [1.0], "hdi_high": [1.4]})
        with pytest.raises(ValueError, match="cannot label the HDI"):
            plot_posterior(bare)
        # a multi-site frame has no single interval to shade, and must say so
        ax_multi = plot_posterior(frame, hdi=0.9)
        assert "no single 90% HDI" in ax_multi.get_title()
        assert ax_multi.get_legend() is None
    finally:
        plt.close("all")


def test_plot_posterior_validates_its_inputs():
    plt = _pyplot()
    from phyloslicer.viz import plot_posterior

    try:
        with pytest.raises(ValueError, match="no finite values"):
            plot_posterior(pd.DataFrame({"rate_median": [np.nan, np.nan]}))
        frame = _posterior_frame(0.95)
        with pytest.raises(ValueError, match="kind must be"):
            plot_posterior(frame, kind="density")
        ax = plot_posterior(frame, kind="hist", hdi=0.95, rate="rate_mean")
        assert len(ax.patches) >= 2  # histogram bars, not an HDI span
    finally:
        plt.close("all")


# ------------------------------------------------------------ 6. CLI posterior
def test_cli_posterior_subcommand_matches_the_library(tmp_path):
    """§5: the ``posterior`` subcommand had no test at all."""
    from phyloslicer import cpd_rate
    from phyloslicer.cli.main import main

    trees = [yule_tree(12, age=3.0, seed=s) for s in range(1, 9)]
    pool = tmp_path / "pool.tre"
    pool.write_text("".join(write_newick(t) + "\n" for t in trees), encoding="utf-8")
    mat = pd.DataFrame(
        np.eye(3, 12) + np.eye(3, 12, k=1),
        columns=trees[0].tip_labels,
        index=["alpha", "beta", "gamma"],
    )
    sites = tmp_path / "sites.csv"
    mat.to_csv(sites)

    def run(out_name: str, extra: list[str]) -> pd.DataFrame:
        out = tmp_path / out_name
        rc = main(
            ["posterior", str(pool), str(sites), "--out", str(out), "--n-slices", "6", *extra]
        )
        assert rc == 0, f"posterior CLI exited with {rc}"
        return pd.read_csv(out)

    narrow = run("post-068.csv", ["--hdi", "0.68"])
    assert list(narrow.columns) == [
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
    ]
    assert list(narrow["site_id"]) == ["alpha", "beta", "gamma"]
    assert (narrow["n_trees"] == len(trees)).all()
    assert (narrow["n_converged"] == len(trees)).all()
    assert np.isfinite(narrow[["rate_mean", "rate_median", "hdi_low", "hdi_high"]].to_numpy()).all()
    assert (narrow["hdi_low"] <= narrow["hdi_high"]).all()
    assert narrow["ess_proxy_caveat"].str.startswith("heuristic only").all()
    # the CLI must run the same kernel as the library over the same pool
    expected = np.nanmean(
        np.vstack([cpd_rate(t, mat, n_slices=6)["CpD"].to_numpy() for t in trees]), axis=0
    )
    np.testing.assert_allclose(narrow["rate_mean"].to_numpy(), expected, rtol=1e-10)

    default = run("post-default.csv", [])
    wide = run("post-095.csv", ["--hdi", "0.95"])
    # --hdi is wired through: a higher level cannot give a narrower interval
    assert (
        narrow["hdi_high"].to_numpy() - narrow["hdi_low"].to_numpy()
        < wide["hdi_high"].to_numpy() - wide["hdi_low"].to_numpy()
    ).all()
    np.testing.assert_allclose(default["hdi_low"].to_numpy(), wide["hdi_low"].to_numpy())
    np.testing.assert_allclose(default["hdi_high"].to_numpy(), wide["hdi_high"].to_numpy())


def test_cli_posterior_reports_a_failure_code_for_an_unusable_pool(tmp_path, capsys):
    from phyloslicer.cli.main import main

    tree = yule_tree(6, age=1.0, seed=2)
    sites = tmp_path / "sites.csv"
    pd.DataFrame([[1.0] * 6], columns=tree.tip_labels, index=["s0"]).to_csv(sites)
    empty = tmp_path / "pool.tre"
    empty.write_text("", encoding="utf-8")
    rc = main(["posterior", str(empty), str(sites), "--out", str(tmp_path / "o.csv")])
    assert rc == 1
    assert "error" in capsys.readouterr().err.lower()


# ----------------------------------------------- 7. align_tree_matrix pruning
def _expected_pruned_edge_index(tree, keep_labels):
    """Independent brute-force reference for the pruned topology.

    Returns ``(kept_edge_indices, mrca)``: every edge on some root -> kept-tip
    path, minus the root -> MRCA stem (``ape::keep.tip(root.edge = 0)``).
    """
    positions = {label: i for i, label in enumerate(tree.tip_labels)}
    parent = tree.parent_of()
    edge_of_child = {int(child): i for i, child in enumerate(tree.edges[:, 1])}
    depths = tree.node_depths()
    spanning: set[int] = set()
    common: set[int] | None = None
    for label in keep_labels:
        node = positions[label]
        path_nodes: list[int] = []
        while node in parent:
            spanning.add(edge_of_child[node])
            path_nodes.append(node)
            node = parent[node]
        path_nodes.append(tree.root)
        common = set(path_nodes) if common is None else (common & set(path_nodes))
    mrca = max(common, key=lambda node: depths[node])
    stem: set[int] = set()
    cursor = mrca
    while cursor in parent:
        stem.add(edge_of_child[cursor])
        cursor = parent[cursor]
    return sorted(spanning - stem), mrca


def test_align_tree_matrix_pruning_keeps_a_self_consistent_tree():
    """§5: the auto-pruning branch never ran, because every committed fixture
    matrix contained every tip.  This checks the pruned topology, node depths,
    ``tip_intervals`` structure and stem removal - not just a tip count."""
    tree = yule_tree(20, age=6.0, seed=61)
    # scattered *inside* one root child, so the kept tips have an MRCA below the
    # source root and the pruning really has to discard a stem
    keep = [tree.tip_labels[i] for i in (4, 7, 9, 12, 15, 18)]
    mat = pd.DataFrame(
        [[1.0, 1.0, 0.0, 1.0, 0.0, 1.0], [0.0, 1.0, 1.0, 1.0, 1.0, 0.0]],
        columns=keep,
        index=["a", "b"],
    )
    with pytest.warns(UserWarning, match="pruning 14 tree tips"):
        cut, aligned = align_tree_matrix(tree, mat)

    assert set(cut.tip_labels) == set(keep)
    assert list(cut.tip_labels) == sorted(keep, key=lambda lab: tree.tip_labels.index(lab))
    assert list(aligned.columns) == list(cut.tip_labels)  # reindexed to the pruned tips
    assert list(aligned.index) == ["a", "b"]
    assert cut.n_tips == len(keep)
    assert cut.is_ultrametric

    depths = cut.node_depths()
    assert np.all(depths[: cut.n_tips] == pytest.approx(cut.root_age, rel=1e-12))
    assert depths[cut.root] == pytest.approx(0.0)

    # every tip edge carries exactly its own unit interval, and the whole set is
    # laminar (nested or disjoint) - the property site membership is built on
    intervals = cut.tip_intervals()
    assert intervals.shape == (cut.n_edges, 2)
    assert (intervals[:, 0] < intervals[:, 1]).all()
    tip_edges = cut.edges[:, 1] < cut.n_tips
    order = np.argsort(cut.edges[tip_edges, 1])
    np.testing.assert_array_equal(
        intervals[tip_edges][order],
        np.column_stack([np.arange(cut.n_tips), np.arange(1, cut.n_tips + 1)]),
    )
    assert intervals.max() == cut.n_tips and intervals.min() == 0
    for a, b in itertools.combinations(map(tuple, intervals), 2):
        nested = (a[0] <= b[0] and b[1] <= a[1]) or (b[0] <= a[0] and a[1] <= b[1])
        disjoint = a[1] <= b[0] or b[1] <= a[0]
        assert nested or disjoint, (a, b)

    # the pruned edge set equals the independent reconstruction, and the
    # root -> MRCA stem is gone rather than folded into the new root
    edge_of_child = {int(c): i for i, c in enumerate(tree.edges[:, 1])}
    expected_edges, mrca = _expected_pruned_edge_index(tree, keep)
    assert mrca != tree.root, "this case must exercise root -> MRCA stem removal"
    np.testing.assert_allclose(np.sort(cut.lengths), np.sort(tree.lengths[expected_edges]))
    assert cut.total_pd() == pytest.approx(tree.lengths[expected_edges].sum(), rel=1e-12)
    missing = set(range(tree.n_edges)) - set(expected_edges)
    assert edge_of_child[mrca] in missing  # the stem edge itself is discarded
    assert tree.total_pd() - cut.total_pd() == pytest.approx(
        tree.lengths[sorted(missing)].sum(), rel=1e-12
    )
    assert cut.root_age == pytest.approx(tree.root_age - tree.node_depths()[mrca], rel=1e-12)
    assert cut.root_age < tree.root_age

    # windows, contributions and PD stay consistent on the pruned tree
    stack = slice_pieces(cut, n=5)
    assert isinstance(stack, SliceStack)
    C = contribution_matrix(cut, stack.depth_windows)
    np.testing.assert_allclose(C.sum(axis=1), cut.lengths, rtol=1e-12, atol=1e-12)
    assert stack.depth_windows[0, 0] == pytest.approx(0.0, abs=1e-12)
    assert stack.depth_windows[-1, 1] == pytest.approx(cut.root_age)
    prof = pd_per_slice(cut, aligned, stack)
    np.testing.assert_allclose(prof.sum(axis=1), pd_total(cut, aligned), rtol=1e-12)


def test_align_tree_matrix_prunes_a_whole_missing_nested_clade():
    """A nested clade absent from the matrix must collapse and re-root deeper."""
    tree = yule_tree(20, age=6.0, seed=70)
    edge_of_child = {int(c): i for i, c in enumerate(tree.edges[:, 1])}
    intervals = tree.tip_intervals()
    # drop *all* tips of the smaller root child: the survivors form one clade, so
    # the pruned tree must be re-rooted on that clade, not on the source root
    root_children = [int(c) for c in tree.edges[tree.edges[:, 0] == tree.root][:, 1]]
    spans = {c: tuple(int(v) for v in intervals[edge_of_child[c]]) for c in root_children}
    dropped_child = min(root_children, key=lambda c: spans[c][1] - spans[c][0])
    kept_child = next(c for c in root_children if c != dropped_child)
    lo, hi = spans[dropped_child]
    dropped = list(tree.tip_labels[lo:hi])
    assert 2 <= len(dropped) <= tree.n_tips - 2, "need a proper nested clade"
    keep = [t for t in tree.tip_labels if t not in dropped]

    mat = pd.DataFrame(np.ones((2, len(keep))), columns=keep, index=["x", "y"])
    with pytest.warns(UserWarning, match=f"pruning {len(dropped)} tree tips"):
        cut, aligned = align_tree_matrix(tree, mat)

    assert set(cut.tip_labels) == set(keep)
    assert not set(cut.tip_labels) & set(dropped)
    assert list(aligned.columns) == list(cut.tip_labels)
    assert cut.is_ultrametric
    # the survivors' MRCA is the untouched root child; its stem is removed
    expected_edges, mrca = _expected_pruned_edge_index(tree, keep)
    assert mrca == kept_child
    assert mrca != tree.root, "the pruned tree must be re-rooted below the source root"
    stem_edge = edge_of_child[mrca]
    assert stem_edge not in expected_edges  # the root -> MRCA stem is discarded
    assert cut.total_pd() == pytest.approx(tree.lengths[expected_edges].sum(), rel=1e-12)
    # what disappeared is the stem *plus* every edge of the dropped clade
    gone = sorted(set(range(tree.n_edges)) - set(expected_edges))
    assert tree.total_pd() - cut.total_pd() == pytest.approx(tree.lengths[gone].sum(), rel=1e-12)
    assert set(gone) > {stem_edge, edge_of_child[dropped_child]}
    assert cut.root_age == pytest.approx(tree.root_age - tree.node_depths()[mrca])
    assert cut.root_age < tree.root_age
    # nothing of the dropped clade survives in the interval bookkeeping
    assert cut.tip_intervals().max() == cut.n_tips
    assert cut.n_edges == len(expected_edges)
    stack = slice_pieces(cut, n=4, criterion="pd")
    C = contribution_matrix(cut, stack.depth_windows)
    np.testing.assert_allclose(C.sum(axis=1), cut.lengths, rtol=1e-12, atol=1e-12)
    assert np.all(C >= 0)
    prof = pd_per_slice(cut, aligned, stack)
    assert np.all(prof >= 0) and np.all(prof <= cut.total_pd())


# ------------------------------------- 8. prune_tips: quantiles as real numbers
# terminal branch lengths of conftest's ``caterpillar5``: A=1, B=1, C=2, D=3, E=4
TERMINAL_CATERPILLAR = [1.0, 1.0, 2.0, 3.0, 4.0]
#: marker for a threshold that cannot produce a sliceable tree
RAISES = object()


@pytest.mark.parametrize(
    "probability, threshold, after, before",
    [
        # side="after" keeps terminal branches >= the quantile (R method=1),
        # "before" keeps the complement (method=2).  RAISES marks a selection
        # that leaves fewer than two tips: prune_tips explains the infeasible
        # threshold instead of returning R's bare NULL.
        (0.25, 1.0, ["A", "B", "C", "D", "E"], RAISES),
        (0.50, 2.0, ["C", "D", "E"], ["A", "B"]),
        (0.75, 3.0, ["D", "E"], ["A", "B", "C"]),
        (0.90, 3.6, RAISES, ["A", "B", "C", "D"]),
        (1.00, 4.0, RAISES, ["A", "B", "C", "D"]),
    ],
)
def test_prune_tips_quantiles_resolve_to_the_declared_lengths(
    probability, threshold, after, before, caterpillar5
):
    """§5 named the quantile test as merely asserting ``kept.n_tips >= 1``.

    The probabilities are pinned to the type-7 quantiles of the terminal branch
    lengths - the interpolation ``stats::quantile`` uses by default, i.e.
    1.0 / 2.0 / 3.0 / 3.6 / 4.0 at p = 0.25 / 0.5 / 0.75 / 0.9 / 1.0 - three
    times over: as numbers, as the tip sets they select, and as the identical
    result of an explicit-length call.
    """
    from phyloslicer import prune_tips

    assert float(np.quantile(TERMINAL_CATERPILLAR, probability)) == pytest.approx(threshold)
    lengths = dict(zip(caterpillar5.tip_labels, TERMINAL_CATERPILLAR, strict=True))
    assert lengths == {"A": 1.0, "B": 1.0, "C": 2.0, "D": 3.0, "E": 4.0}
    for expected, side in ((after, "after"), (before, "before")):
        if expected is RAISES:
            with pytest.raises(ValueError, match="at least 2 tips"):
                prune_tips(caterpillar5, probability, quantiles=True, side=side)
            with pytest.raises(ValueError, match="at least 2 tips"):
                prune_tips(caterpillar5, threshold, side=side)
            continue
        by_q = prune_tips(caterpillar5, probability, quantiles=True, side=side)
        by_value = prune_tips(caterpillar5, threshold, side=side)
        assert sorted(by_q.tip_labels) == expected
        assert sorted(by_value.tip_labels) == expected
        if side == "after":
            assert all(lengths[t] >= threshold for t in by_q.tip_labels)
        else:
            assert all(lengths[t] < threshold for t in by_q.tip_labels)
    if after is not RAISES and before is not RAISES:
        assert set(after) | set(before) == set(caterpillar5.tip_labels)
        assert not set(after) & set(before)


def test_prune_tips_rejects_probabilities_outside_the_unit_interval(caterpillar5):
    from phyloslicer import prune_tips

    with pytest.raises(ValueError, match="quantile values"):
        prune_tips(caterpillar5, 1.5, quantiles=True)
    with pytest.raises(ValueError, match="quantile values"):
        prune_tips(caterpillar5, -0.1, quantiles=True)
    with pytest.raises(ValueError, match="side must be"):
        prune_tips(caterpillar5, 0.5, quantiles=True, side="around")
