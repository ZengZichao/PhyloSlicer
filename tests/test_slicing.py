"""Slicing API and interval-intersection kernels, against analytic answers."""

from __future__ import annotations

import numpy as np
import pytest

from phyloslicer.compat import phylo_pieces, squeeze_int, squeeze_root
from phyloslicer.slicing import (
    depth_for_pd,
    pd_depth_curve,
    prune_tips,
    slice_interval,
    slice_pieces,
    slice_rootward,
    slice_tipward,
)


class TestAnalyticSlices:
    def test_rootward_full_window(self, balanced4):
        s = slice_rootward(balanced4, 2.0)  # whole tree
        assert s.pd == pytest.approx(6.0)

    def test_rootward_young(self, balanced4):
        s = slice_rootward(balanced4, 1.0)  # keep [1, 2] depth = 4 tip branches
        assert s.pd == pytest.approx(4.0)
        assert (s.lengths[[0, 1]] == 0).all()  # root edges removed
        assert (s.lengths[[2, 3, 4, 5]] == 1).all()

    def test_tipward_old(self, balanced4):
        s = slice_tipward(balanced4, 1.0)  # keep [0, 1] = root edges
        assert s.pd == pytest.approx(2.0)

    def test_complement_property(self, random_tree_20):
        """slice_rootward(t) and slice_tipward(t) partition the tree."""
        tree = random_tree_20
        t = 3.7
        young = slice_rootward(tree, t).lengths
        old = slice_tipward(tree, t).lengths
        assert np.allclose(young + old, tree.lengths)

    def test_interval(self, caterpillar5):
        s = slice_interval(caterpillar5, start=3.0, stop=1.0)  # depth [1, 3]
        assert s.pd == pytest.approx(7.0)
        with pytest.raises(ValueError):
            slice_interval(caterpillar5, start=1.0, stop=3.0)

    def test_interval_invert(self, caterpillar5):
        older, younger = slice_interval(caterpillar5, 3.0, 1.0, invert=True)
        assert older.pd == pytest.approx(2.0)  # depth [0, 1]
        assert younger.pd == pytest.approx(5.0)  # depth [3, 4]
        interior = slice_interval(caterpillar5, 3.0, 1.0)
        assert older.pd + younger.pd + interior.pd == pytest.approx(caterpillar5.total_pd())

    def test_pd_criterion_rootward(self, balanced4):
        # remove root-side PD of 2 -> keep 4
        s = slice_rootward(balanced4, 2.0, criterion="pd")
        assert s.pd == pytest.approx(4.0)

    def test_pd_criterion_tipward(self, balanced4):
        # keep root-side PD of 2
        s = slice_tipward(balanced4, 2.0, criterion="pd")
        assert s.pd == pytest.approx(2.0)

    def test_threshold_validation(self, balanced4):
        with pytest.raises(ValueError):
            slice_rootward(balanced4, 2.5)  # >= tree age
        with pytest.raises(ValueError):
            slice_rootward(balanced4, 7.0, criterion="pd")  # > total PD
        with pytest.raises(ValueError):
            slice_rootward(balanced4, -1.0)

    def test_interval_bounds_validation(self, balanced4):
        T = balanced4.root_age
        # time criterion: interval must lie inside [0, T]
        with pytest.raises(ValueError, match="outside"):
            slice_interval(balanced4, start=T + 0.5, stop=0.5)
        with pytest.raises(ValueError, match="outside"):
            slice_interval(balanced4, start=1.5, stop=-0.5)
        # pd criterion: root-side amounts must lie inside [0, total PD]
        total = balanced4.total_pd()
        with pytest.raises(ValueError, match="outside"):
            slice_interval(balanced4, start=-0.5, stop=1.0, criterion="pd")
        with pytest.raises(ValueError, match="outside"):
            slice_interval(balanced4, start=1.0, stop=total + 1.0, criterion="pd")


class TestPdCurve:
    def test_monotone_and_total(self, random_tree_20):
        breaks, cum = pd_depth_curve(random_tree_20)
        assert (np.diff(cum) >= 0).all()
        assert cum[-1] == pytest.approx(random_tree_20.total_pd())

    def test_inverse_roundtrip(self, random_tree_20):
        total = random_tree_20.total_pd()
        for q in (0.1, 0.35, 0.5, 0.9):
            x = depth_for_pd(random_tree_20, q * total)
            breaks, cum = pd_depth_curve(random_tree_20)
            # re-evaluate cumulative PD at x by linear interpolation
            v = np.interp(x, breaks, cum)
            assert v == pytest.approx(q * total, rel=1e-10)


class TestSlicePieces:
    def test_equal_time_pieces_partition_pd(self, random_tree_20):
        stack = slice_pieces(random_tree_20, n=10)
        contrib = stack.contribution()
        assert contrib.sum() == pytest.approx(random_tree_20.total_pd())
        # each slice has equal width
        assert np.allclose(stack.widths, random_tree_20.root_age / 10)

    def test_windows_order_and_ages(self, random_tree_20):
        stack = slice_pieces(random_tree_20, n=4)
        T = random_tree_20.root_age
        expected = np.array([[T, 3 * T / 4], [3 * T / 4, T / 2], [T / 2, T / 4], [T / 4, 0]])
        assert np.allclose(stack.windows, expected)
        assert stack.ages[0] == pytest.approx(T)

    def test_width_rounding_like_r(self, random_tree_20):
        """The rounding rule is decisive, so assert which piece count it selects.

        ``assert n_slices in {floor, ceil}`` (the previous assertion) is true for
        either branch of the rule and for a wrong implementation too; the counts
        and widths below are what treesliceR's ``phylo_pieces`` actually picks.
        """
        T = random_tree_20.root_age
        assert T == pytest.approx(10.0)  # yule_tree(age=10)

        # 7.3 nominal pieces.  The two achievable widths are T/7 = 1.428571 and
        # T/8 = 1.25; the request T/7.3 = 1.369863 is nearer T/7, so the floor
        # branch wins and every slice is 10/7 Ma wide.
        stack = slice_pieces(random_tree_20, width=T / 7.3)
        assert stack.n_slices == 7
        assert np.allclose(stack.widths, T / 7, rtol=1e-12)
        # windows are reported root -> tip, so they span exactly [0, T]
        assert stack.windows.min() == 0.0
        assert stack.windows.max() == pytest.approx(T)

        # 7.8 nominal pieces: T/7.8 = 1.282051 is nearer T/8 -> ceiling branch.
        stack = slice_pieces(random_tree_20, width=T / 7.8)
        assert stack.n_slices == 8
        assert np.allclose(stack.widths, T / 8, rtol=1e-12)

        # exactly between the two candidates (7.5) the rule takes the smaller
        # width, i.e. more slices; a width equal to the root age is one slice.
        assert slice_pieces(random_tree_20, width=T / 7.5).n_slices == 8
        assert slice_pieces(random_tree_20, width=T / 6.0).n_slices == 6
        single = slice_pieces(random_tree_20, width=T)
        assert single.n_slices == 1
        assert single.widths[0] == pytest.approx(T)
        # every slice of a time stack has the same width and they tile [0, T]
        stack = slice_pieces(random_tree_20, width=T / 7.3)
        assert stack.widths.sum() == pytest.approx(T)
        assert np.allclose(np.abs(np.diff(stack.windows[:, 0])), T / 7)

    def test_pd_criterion_pieces(self, random_tree_20):
        stack = slice_pieces(random_tree_20, n=5, criterion="pd")
        # each slice must store equal PD
        contrib = stack.contribution()
        assert np.allclose(contrib, contrib[0], rtol=1e-9)

    def test_contribution_subset(self, balanced4):
        stack = slice_pieces(balanced4, n=2)
        c_all = stack.contribution()
        c_ab = stack.contribution(["A", "B"])
        # window 1 (root side) contains the root edges shared by A and B
        assert c_all[0] == pytest.approx(2.0)
        assert c_ab[0] == pytest.approx(1.0)
        assert c_ab[1] == pytest.approx(2.0)

    def test_neither_n_nor_width_raises(self, balanced4):
        with pytest.raises(ValueError):
            slice_pieces(balanced4)
        with pytest.raises(ValueError):
            slice_pieces(balanced4, n=3, width=0.5)

    def test_to_sliced_trees(self, random_tree_20):
        stack = slice_pieces(random_tree_20, n=3)
        pieces = stack.to_sliced_trees()
        assert len(pieces) == 3
        total = sum(p.pd for p in pieces)
        assert total == pytest.approx(random_tree_20.total_pd())


class TestPruneTips:
    def test_side_after(self, random_tree_20):
        tip_edge = random_tree_20.edges[:, 1] < random_tree_20.n_tips
        terminal = np.zeros(random_tree_20.n_tips)
        terminal[random_tree_20.edges[tip_edge, 1]] = random_tree_20.lengths[tip_edge]
        thr = float(np.median(terminal))
        kept = prune_tips(random_tree_20, thr, side="after")
        assert kept.n_tips == int((terminal >= thr).sum())

    def test_side_before(self, random_tree_20):
        tip_edge = random_tree_20.edges[:, 1] < random_tree_20.n_tips
        terminal = np.zeros(random_tree_20.n_tips)
        terminal[random_tree_20.edges[tip_edge, 1]] = random_tree_20.lengths[tip_edge]
        thr = float(np.median(terminal))
        kept = prune_tips(random_tree_20, thr, side="before")
        assert kept.n_tips == int((terminal < thr).sum())

    def test_quantiles(self, caterpillar5):
        """A quantile threshold selects named tips, not merely "some" tips.

        ``caterpillar5`` has terminal branch lengths A=1, B=1, C=2, D=3, E=4, so
        the type-7 quantiles - the same linear interpolation ``stats::quantile``
        uses by default in R - are 1.0, 2.0, 3.0 and 3.6 at p = 0.25, 0.5, 0.75
        and 0.9.
        """
        terminal = {"A": 1.0, "B": 1.0, "C": 2.0, "D": 3.0, "E": 4.0}
        observed = {
            caterpillar5.tip_labels[int(c)]: float(length)
            for (_, c), length in zip(caterpillar5.edges.tolist(), caterpillar5.lengths.tolist())
            if c < caterpillar5.n_tips
        }
        assert observed == pytest.approx(terminal)

        # side="after" keeps terminal branches >= the quantile (R method=1)
        assert sorted(prune_tips(caterpillar5, 0.25, quantiles=True).tip_labels) == [
            "A",
            "B",
            "C",
            "D",
            "E",
        ]
        assert sorted(prune_tips(caterpillar5, 0.5, quantiles=True).tip_labels) == ["C", "D", "E"]
        assert sorted(prune_tips(caterpillar5, 0.75, quantiles=True).tip_labels) == ["D", "E"]
        # side="before" is the complement (R method=2)
        assert sorted(prune_tips(caterpillar5, 0.75, quantiles=True, side="before").tip_labels) == [
            "A",
            "B",
            "C",
        ]
        # a selection of fewer than two tips raises, naming the threshold and
        # the usable range, rather than returning a bare None (treesliceR's NULL)
        with pytest.raises(ValueError, match="needs at least 2 tips"):
            prune_tips(caterpillar5, 0.25, quantiles=True, side="before")
        with pytest.raises(ValueError, match="needs at least 2 tips"):
            prune_tips(caterpillar5, 0.9, quantiles=True)
        # a vector of probabilities gives one result per entry, in order
        many = prune_tips(caterpillar5, [0.25, 0.5, 0.75], quantiles=True)
        assert [t.n_tips for t in many] == [5, 3, 2]
        # the survivors stay ultrametric, at the source root age
        after = prune_tips(caterpillar5, 0.5, quantiles=True)
        assert after.is_ultrametric
        assert after.root_age == pytest.approx(caterpillar5.root_age)


class TestCompatLayer:
    def test_squeeze_root_matches_slice_rootward(self, random_tree_20):
        a = squeeze_root(random_tree_20, 2.0)
        b = slice_rootward(random_tree_20, 2.0)
        assert np.allclose(a.lengths, b.lengths)

    def test_squeeze_int(self, random_tree_20):
        a = squeeze_int(random_tree_20, 3.0, 1.0)
        b = slice_interval(random_tree_20, 3.0, 1.0)
        assert np.allclose(a.lengths, b.lengths)

    def test_phylo_pieces_r_semantics(self, random_tree_20):
        pieces = phylo_pieces(random_tree_20, n=4, criterion="my", method=1)
        assert len(pieces) == 4
        out = phylo_pieces(random_tree_20, n=4, timeSteps=True, returnTree=True)
        assert len(out) == 3  # pieces, timesteps, original tree

    def test_phylo_pieces_pd(self, random_tree_20):
        pieces = phylo_pieces(random_tree_20, n=5, criterion="pd", method=1)
        vals = [p.pd for p in pieces]
        assert np.allclose(vals, vals[0], rtol=1e-9)
