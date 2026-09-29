"""Indices: PD / PE / PB against brute-force reference implementations."""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
import pytest

from phyloslicer.core import TreeArray
from phyloslicer.indices import (
    MULTISITE_DOMAINS,
    edge_range_sizes,
    pb,
    pb_per_slice,
    pd_per_slice,
    pe,
    pe_per_slice,
    r_compatible_profile,
    site_edge_membership,
)
from phyloslicer.indices import (
    pd as faith_pd,
)
from phyloslicer.rates.api import cpb_rate
from phyloslicer.simulate import yule_tree
from phyloslicer.slicing import slice_pieces


# ---------------------------------------------------------------------- #
# brute-force references (slow, obviously-correct versions)
# ---------------------------------------------------------------------- #
def brute_force_path_edges(tree, tip):
    parent = tree.parent_of()
    edges = set()
    node = tip
    while node in parent:
        for e, (p, c) in enumerate(tree.edges.tolist()):
            if c == node:
                edges.add(e)
                node = p
                break
    return edges


def brute_force_pd_per_slice(tree, mat, stack, weighted=False):
    n_sites = mat.shape[0]
    k = stack.n_slices
    C = stack.contribution_matrix()
    out = np.zeros((n_sites, k))
    for s in range(n_sites):
        present = np.flatnonzero(mat.iloc[s].to_numpy() > 0)
        path = set()
        for t in present:
            path |= brute_force_path_edges(tree, int(t))
        for j in range(k):
            vals = [C[e, j] for e in path]
            if weighted:
                ranges = edge_range_sizes(tree, mat)
                vals = [v / ranges[e] for e, v in zip(path, vals)]
            out[s, j] = sum(vals)
    return out


def brute_force_neighbourhood_beta(tree, mat, adj, stack, component, weighted=False):
    """Literal translation of the treesliceR arithmetic (combn + sums)."""
    n_sites = mat.shape[0]
    k = stack.n_slices
    C = stack.contribution_matrix()
    if weighted:
        ranges = edge_range_sizes(tree, mat)
        lengths = tree.lengths / ranges
    else:
        lengths = tree.lengths
    out = np.full((n_sites, k), np.nan)
    totals = np.full(n_sites, np.nan)

    def pd_of(mask_edges):
        return float(sum(lengths[e] for e in mask_edges))

    def pd_slice_of(mask_edges, j):
        return float(sum(C[e, j] for e in mask_edges))

    for focal in range(n_sites):
        members = np.flatnonzero(adj[focal] == 1)
        if members.size < 2:
            continue
        paths = [
            set().union(
                *[
                    brute_force_path_edges(tree, int(t))
                    for t in np.flatnonzero(mat.iloc[m].to_numpy() > 0)
                ]
                or set()
            )
            for m in members
        ]
        pd_single = np.array([pd_of(p) for p in paths])
        pairs = list(itertools.combinations(range(members.size), 2))
        pw_paths = [paths[i] | paths[j] for i, j in pairs]
        pw_pd = np.array([pd_of(p) for p in pw_paths])
        ii = np.array([i for i, _ in pairs])
        jj = np.array([j for _, j in pairs])
        _shared = pd_single[ii] + pd_single[jj] - pw_pd  # documented identity

        if component == "sorensen":
            # sum over the R-recycled (P x 2) matrix => sum_pairs (2*pw - pd_i - pd_j)
            bc_tot = float((2 * pw_pd - (pd_single[ii] + pd_single[jj])).sum())
            a_tot = float(pd_single.sum() - pd_of(set().union(*paths)))
            den = 2 * a_tot + bc_tot
            if den == 0:
                continue
            totals[focal] = bc_tot / den
            for j in range(k):
                pd_piece = np.array([pd_slice_of(p, j) for p in paths])
                pw_piece = np.array([pd_slice_of(p, j) for p in pw_paths])
                r_bc = float((2 * pw_piece - (pd_piece[ii] + pd_piece[jj])).sum())
                out[focal, j] = (r_bc / den) / totals[focal]
        else:
            raise NotImplementedError
    return out, totals


class TestMembership:
    def test_membership_matches_bruteforce(self, random_tree_20, random_mat_20):
        A = site_edge_membership(random_tree_20, random_mat_20).toarray()
        for s in range(random_mat_20.shape[0]):
            present = np.flatnonzero(random_mat_20.iloc[s].to_numpy() > 0)
            path = set()
            for t in present:
                path |= brute_force_path_edges(random_tree_20, int(t))
            brute = np.zeros(random_tree_20.n_edges)
            brute[list(path)] = 1.0
            assert np.allclose(A[s], brute)

    def test_range_sizes(self, balanced4):
        mat = pd.DataFrame(
            [[1, 0, 1, 0], [1, 1, 0, 0]], columns=["A", "B", "C", "D"], index=["s1", "s2"]
        )
        ranges = edge_range_sizes(balanced4, mat)
        # edges: [root->AB, root->CD, AB->A, AB->B, CD->C, CD->D]
        expected = np.array([2, 1, 2, 1, 1, 0], dtype=float)
        assert np.array_equal(ranges, expected)


class TestPd:
    def test_total_pd(self, random_tree_20, random_mat_20):
        fast = faith_pd(random_tree_20, random_mat_20)
        for s in range(random_mat_20.shape[0]):
            present = np.flatnonzero(random_mat_20.iloc[s].to_numpy() > 0)
            path = set()
            for t in present:
                path |= brute_force_path_edges(random_tree_20, int(t))
            brute = random_tree_20.lengths[list(path)].sum()
            assert fast[s] == pytest.approx(brute)

    def test_pd_per_slice_matches_bruteforce(self, random_tree_20, random_mat_20):
        stack = slice_pieces(random_tree_20, n=7)
        fast = pd_per_slice(random_tree_20, random_mat_20, stack)
        brute = brute_force_pd_per_slice(random_tree_20, random_mat_20, stack)
        assert np.allclose(fast, brute, atol=1e-10)

    def test_pe_per_slice_matches_weighted_bruteforce(self, random_tree_20, random_mat_20):
        stack = slice_pieces(random_tree_20, n=6)
        fast = pe_per_slice(random_tree_20, random_mat_20, stack)
        brute = brute_force_pd_per_slice(random_tree_20, random_mat_20, stack, weighted=True)
        assert np.allclose(fast, brute, atol=1e-10)

    def test_profiles_sum_to_totals(self, random_tree_20, random_mat_20):
        stack = slice_pieces(random_tree_20, n=9)
        prof = pd_per_slice(random_tree_20, random_mat_20, stack)
        assert np.allclose(prof.sum(axis=1), faith_pd(random_tree_20, random_mat_20))


class TestBeta:
    @pytest.fixture
    def small_setup(self):
        tree = yule_tree(10, age=1.0, seed=5)
        rng = np.random.default_rng(2)
        mat = pd.DataFrame(
            (rng.random((6, 10)) < 0.4).astype(float),
            columns=tree.tip_labels,
            index=[f"s{i}" for i in range(6)],
        )
        mat.iloc[0] = 1.0
        adj = np.eye(6, dtype=int)
        adj[0, 1] = adj[1, 0] = 1
        adj[1, 2] = adj[2, 1] = 1
        adj[2, 3] = adj[3, 2] = 1
        adj[3, 4] = adj[4, 3] = 1
        adj[4, 5] = adj[5, 4] = 1
        return tree, mat, adj

    def test_multisite_sorensen_matches_bruteforce(self, small_setup):
        tree, mat, adj = small_setup
        stack = slice_pieces(tree, n=5)
        res = pb_per_slice(tree, mat, adj, stack, component="sorensen", approach="multisite")
        brute, totals = brute_force_neighbourhood_beta(tree, mat, adj, stack, component="sorensen")
        mask = np.isfinite(brute)
        assert mask.any()
        assert np.allclose(res["values"][mask], brute[mask], atol=1e-10)
        mask_focal = np.isfinite(totals)
        assert np.allclose(res["total"][mask_focal], totals[mask_focal])

    def test_totals_consistent_between_pb_and_pb_per_slice(self, small_setup):
        tree, mat, adj = small_setup
        stack = slice_pieces(tree, n=4)
        res = pb_per_slice(tree, mat, adj, stack, component="sorensen", approach="multisite")
        totals = pb(tree, mat, adj, component="sorensen", approach="multisite")
        mask = np.isfinite(res["total"])
        fast = totals.set_index("site_id")["PB"].to_numpy()[mask]
        assert np.allclose(res["total"][mask], fast)

    def test_components_complementarity(self, small_setup):
        """sorensen = turnover + nestedness for the multisite decomposition."""
        tree, mat, adj = small_setup
        stack = slice_pieces(tree, n=4)
        sor = pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite")["total"]
        tur = pb_per_slice(tree, mat, adj, stack, "turnover", "multisite")["total"]
        nest = pb_per_slice(tree, mat, adj, stack, "nestedness", "multisite")["total"]
        mask = np.isfinite(sor) & np.isfinite(tur) & np.isfinite(nest)
        assert np.allclose(sor[mask], tur[mask] + nest[mask], atol=1e-10)

    def test_weighted_runs(self, small_setup):
        tree, mat, adj = small_setup
        stack = slice_pieces(tree, n=4)
        res = pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite", weighted=True)
        assert np.isfinite(res["total"]).any()

    def test_degenerate_neighbourhood_flagged(self, small_setup):
        tree, mat, adj = small_setup
        mat2 = mat.copy()
        mat2.iloc[0, :] = 0.0
        mat2.iloc[1, :] = 0.0  # the whole focal+neighbourhood is empty
        adj2 = np.eye(mat.shape[0], dtype=int)
        adj2[0, 1] = adj2[1, 0] = 1
        stack = slice_pieces(tree, n=4)
        res = pb_per_slice(tree, mat2, adj2, stack, "sorensen", "multisite")
        assert res["status"][0] == "degenerate_neighbourhood"
        assert np.isnan(res["total"][0])


# ---------------------------------------------------------------------- #
# multisite beta: summation domain of the shared term
# ---------------------------------------------------------------------- #
def broom(m: int) -> TreeArray:
    """One shared skeleton edge of length 1 plus one private terminal edge per tip.

    ``root -(1)-> A -(1)-> t1 .. tm``: site *i* holds tip *i* only, so every pair
    has ``a = 1``, ``b = c = 1`` (Sørensen 0.5) no matter how many sites the
    neighbourhood holds.  This is the construction the review used to show that
    summing ``a`` over sites while summing ``b + c`` over pairs makes the
    multisite index equal ``m / (m + 2)`` instead of staying at 0.5.
    """
    edges = np.array([[m, m + 1]] + [[m + 1, i] for i in range(m)], dtype=int)
    return TreeArray([f"t{i + 1}" for i in range(m)], edges, np.ones(m + 1))


def broom_setup(m: int):
    tree = broom(m)
    mat = pd.DataFrame(np.eye(m), columns=tree.tip_labels, index=[f"p{i}" for i in range(m)])
    adj = np.ones((m, m), dtype=int)  # one neighbourhood of exactly m sites
    return tree, mat, adj


def chain_adjacency(n: int) -> np.ndarray:
    """Sites 0-1-2-...-n-1 as a queen-style chain, each with itself and two mates."""
    adj = np.eye(n, dtype=int)
    for i in range(n - 1):
        adj[i, i + 1] = adj[i + 1, i] = 1
    return adj


class TestMultisiteSummationDomain:
    @pytest.mark.parametrize("m", [2, 3, 5, 10, 20])
    def test_default_follows_the_documented_m_over_m_plus_2_drift(self, m):
        """Reference arithmetic: 0.5 at m=2, 0.909 at m=20, i.e. m/(m+2)."""
        tree, mat, adj = broom_setup(m)
        for component in ("sorensen", "turnover"):
            got = pb(tree, mat, adj, component, "multisite")["PB"].iloc[0]
            assert got == pytest.approx(m / (m + 2), abs=1e-12), (component, m)
        # the profile path (what cpb_rate fits) reports the same total
        stack = slice_pieces(tree, n=4)
        prof = pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite")
        assert prof["total"][0] == pytest.approx(m / (m + 2), abs=1e-12)
        assert prof["values"][0].sum() == pytest.approx(1.0, abs=1e-12)

    def test_default_is_the_reference_arithmetic_and_is_unchanged(self):
        """The option must not move the default: mixed == reference == no kwarg."""
        tree, mat, adj = broom_setup(6)
        stack = slice_pieces(tree, n=5)
        for domain in (None, "mixed", "reference"):
            kwargs = {} if domain is None else {"multisite_domain": domain}
            got = pb(tree, mat, adj, "sorensen", "multisite", **kwargs)["PB"].iloc[0]
            prof = pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite", **kwargs)
            assert got == pytest.approx(6 / 8, abs=1e-12)
            assert prof["total"][0] == pytest.approx(6 / 8, abs=1e-12)

    @pytest.mark.parametrize("m", [3, 5, 20])
    @pytest.mark.parametrize("component", ["sorensen", "turnover"])
    def test_paired_domain_agrees_with_pairwise_aggregation(self, m, component):
        """With every term summed over pairs, m > 2 multisite == pairwise == 0.5."""
        tree, mat, adj = broom_setup(m)
        ms = pb(tree, mat, adj, component, "multisite", multisite_domain="paired")["PB"].iloc[0]
        pw = pb(tree, mat, adj, component, "pairwise")["PB"].iloc[0]
        mixed = pb(tree, mat, adj, component, "multisite")["PB"].iloc[0]
        assert ms == pytest.approx(0.5, abs=1e-12)
        assert ms == pytest.approx(pw, abs=1e-12)
        assert abs(mixed - 0.5) > 0.05  # the drift the option exists to remove

    def test_paired_domain_leaves_two_site_neighbourhoods_untouched(self):
        """m = 2 has a single pair, so both domains must coincide exactly."""
        tree, mat, adj = broom_setup(2)
        stack = slice_pieces(tree, n=5)
        for component in ("sorensen", "turnover", "nestedness"):
            for approach in ("multisite", "pairwise"):
                mixed = pb(tree, mat, adj, component, approach)["PB"].iloc[0]
                paired = pb(tree, mat, adj, component, approach, multisite_domain="paired")[
                    "PB"
                ].iloc[0]
                assert paired == pytest.approx(mixed, abs=1e-15)
                prof_m = pb_per_slice(tree, mat, adj, stack, component, approach)
                prof_p = pb_per_slice(
                    tree, mat, adj, stack, component, approach, multisite_domain="paired"
                )
                # nestedness is identically 0 here, so both sides are NaN rows
                assert np.allclose(prof_m["values"], prof_p["values"], atol=1e-14, equal_nan=True)
                assert np.array_equal(prof_m["status"], prof_p["status"])

    def test_r_compatible_profile_keeps_the_reference_summation(self):
        """The R-facing curve is unaffected by the paired domain by design."""
        tree, mat, adj = broom_setup(7)
        stack = slice_pieces(tree, n=5)
        mixed = r_compatible_profile(tree, mat, adj, stack, "sorensen", "multisite")
        assert mixed["values"][0].sum() == pytest.approx(1.0, abs=1e-12)
        paired = r_compatible_profile(
            tree, mat, adj, stack, "sorensen", "multisite", multisite_domain="paired"
        )
        assert paired["values"][0].sum() == pytest.approx(1.0, abs=1e-12)
        assert paired["total"][0] == pytest.approx(0.5, abs=1e-12)
        assert mixed["total"][0] == pytest.approx(7 / 9, abs=1e-12)

    def test_cpb_rate_threads_the_domain(self):
        """cpb_rate(multisite_domain=...) must reach the same index as pb()."""
        tree, mat, adj = broom_setup(8)
        default = cpb_rate(tree, mat, adj, n_slices=4, component="sorensen")
        paired = cpb_rate(
            tree, mat, adj, n_slices=4, component="sorensen", multisite_domain="paired"
        )
        reference = cpb_rate(
            tree, mat, adj, n_slices=4, component="sorensen", multisite_domain="reference"
        )
        assert default["PB"].iloc[0] == pytest.approx(8 / 10, abs=1e-12)
        assert reference["PB"].iloc[0] == pytest.approx(default["PB"].iloc[0], abs=1e-15)
        assert paired["PB"].iloc[0] == pytest.approx(0.5, abs=1e-12)
        assert np.isfinite(paired["CpB"]).all()

    def test_unknown_domain_is_rejected_everywhere(self):
        tree, mat, adj = broom_setup(4)
        stack = slice_pieces(tree, n=3)
        with pytest.raises(ValueError, match="multisite_domain"):
            pb(tree, mat, adj, "sorensen", "multisite", multisite_domain="baselga")
        with pytest.raises(ValueError, match="multisite_domain"):
            pb_per_slice(tree, mat, adj, stack, "sorensen", "multisite", multisite_domain="")
        with pytest.raises(ValueError, match="multisite_domain"):
            cpb_rate(tree, mat, adj, n_slices=3, multisite_domain="site")
        # the canonical spellings stay exactly the two documented domains
        assert MULTISITE_DOMAINS == ("mixed", "paired")


class TestIndexInvariants:
    """Bounds the review listed as never asserted (section 5, test gaps)."""

    def test_pd_per_slice_bounded_by_site_and_tree_totals(self, random_tree_20, random_mat_20):
        stack = slice_pieces(random_tree_20, n=11)
        prof = pd_per_slice(random_tree_20, random_mat_20, stack)
        site_total = faith_pd(random_tree_20, random_mat_20)
        assert (prof >= -1e-15).all()
        assert prof.sum(axis=1) == pytest.approx(site_total, abs=1e-12)
        assert (prof <= site_total[:, None] + 1e-12).all()
        assert (prof <= random_tree_20.total_pd() + 1e-12).all()
        assert prof.sum() == pytest.approx(site_total.sum(), abs=1e-12)

    def test_pe_per_slice_bounded_by_site_totals(self, random_tree_20, random_mat_20):
        stack = slice_pieces(random_tree_20, n=11)
        prof = pe_per_slice(random_tree_20, random_mat_20, stack)
        site_total = pe(random_tree_20, random_mat_20)
        assert (prof >= -1e-15).all()
        assert prof.sum(axis=1) == pytest.approx(site_total, abs=1e-12)
        assert (prof <= site_total[:, None] + 1e-12).all()

    @pytest.mark.parametrize("component", ["sorensen", "turnover", "nestedness"])
    @pytest.mark.parametrize("approach", ["multisite", "pairwise"])
    @pytest.mark.parametrize("domain", ["mixed", "paired"])
    def test_pb_lies_in_unit_interval(
        self, random_tree_20, random_mat_20, component, approach, domain
    ):
        adj = chain_adjacency(random_mat_20.shape[0])
        stack = slice_pieces(random_tree_20, n=6)
        totals = pb(
            random_tree_20, random_mat_20, adj, component, approach, multisite_domain=domain
        )["PB"].to_numpy()
        finite = np.isfinite(totals)
        assert finite.any()
        assert (totals[finite] >= -1e-12).all() and (totals[finite] <= 1 + 1e-12).all()
        prof = pb_per_slice(
            random_tree_20, random_mat_20, adj, stack, component, approach, multisite_domain=domain
        )
        ok = np.isfinite(prof["total"])
        assert ok.any()
        assert (prof["total"][ok] >= -1e-12).all() and (prof["total"][ok] <= 1 + 1e-12).all()
