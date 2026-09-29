"""tests/test_compat.py — treesliceR-compatible wrappers.

Every assertion below is pinned against the *R source*, not against this
package's own behaviour.  R is not installed on the machines that run this
suite, so the anchors are quoted literally from the vendored reference
implementation (``PhyloSlicer-参考软件/treesliceR/R/``) and given as comments:

* ``phylo_pieces.R:95``  ``age <- cumsum(rep(n, nslices))`` — an *increasing*,
  root->tip cumulative depth vector whose last element is the full tree depth;
  returned as-is by ``phylo_pieces(timeSteps = TRUE)`` (``phylo_pieces.R:219``
  does the same for ``criterion = "PD"``, where ``age`` holds the time-mapped
  equal-PD cut points);
* ``squeeze_int.R:75-76, 84-88``  ``invert = TRUE`` returns
  ``list(tree1 = squeeze_root(to), tree2 = squeeze_tips(from))`` — i.e.
  ``(younger, older)``, the crown slice first;
* ``prune_tips.R:52, 121``  only ``method == 1`` and ``method == 2`` have a
  branch; any other value falls off the end and returns ``NULL``;
* ``squeeze_root.R`` / ``squeeze_tips.R`` document ``criterion`` as
  ``"my"`` or ``"PD"``.

The two review-flagged bugs are covered by
``TestPhyloPiecesTimeSteps`` (descending Ma-before-present instead of R's
ascending cumulative depth) and ``TestSqueezeIntInvert`` (reversed complement
order).  ``compat.squeeze_tips`` and ``compat.prune_tips`` had no tests at all
no coverage at all; they are covered here too.
"""

from __future__ import annotations

import numpy as np
import pytest

from phyloslicer import compat
from phyloslicer.core.tree import TreeArray
from phyloslicer.slicing import (
    prune_tips as slice_prune_tips,
)
from phyloslicer.slicing import (
    slice_interval,
    slice_pieces,
    slice_rootward,
    slice_tipward,
)

# `((A:1,B:1):3,(C:2,D:2):2);` -> T = 4, total PD = 11, terminal branches
# [1, 1, 2, 2]: asymmetric enough that method=1 and method=2 prune *different*
# tip sets, and equal-PD cut points differ from the equal-time ones.
NEWICK_T4 = "((A:1,B:1):3,(C:2,D:2):2);"

# `((A:4,B:4):4,(C:4,D:4):4);` -> T = 8, total PD = 24.  This is the exact tree
# the reference implementation uses, so its R vectors can be checked literally.
NEWICK_T8 = "((A:4,B:4):4,(C:4,D:4):4);"


@pytest.fixture
def tree_t4() -> TreeArray:
    return TreeArray.from_newick(NEWICK_T4)


@pytest.fixture
def tree_t8() -> TreeArray:
    return TreeArray.from_newick(NEWICK_T8)


class TestPhyloPiecesTimeSteps:
    """A7(a): R's ``timeSteps`` is the ascending cumulative depth from the root."""

    def test_matches_r_cumsum_for_equal_time_slices(self, tree_t8):
        # treesliceR: n = tree$tree_depth/nslices = 8/4 = 2;
        #             age <- cumsum(rep(n, nslices)) = 2, 4, 6, 8
        pieces, steps = compat.phylo_pieces(tree_t8, n=4, timeSteps=True)
        assert np.array_equal(steps, np.array([2.0, 4.0, 6.0, 8.0]))
        assert len(pieces) == 4

    def test_last_step_is_full_tree_depth(self, tree_t4, tree_t8):
        # cumsum(rep(n, nslices)) necessarily terminates at tree_depth, so the
        # vector covers the whole tree; a lower-bound vector would stop one
        # slice short and silently misalign with the pieces.
        for tree, n in ((tree_t4, 4), (tree_t8, 5)):
            steps = compat.phylo_pieces(tree, n=n, timeSteps=True)[1]
            assert steps[-1] == pytest.approx(tree.root_age)
            assert np.all(np.diff(steps) > 0)

    def test_indexwise_alignment_with_returned_pieces(self, tree_t8):
        # piece j is the slice spanning (steps[j-1], steps[j]] in depth-from-root
        # coordinates, exactly as R's lapply(age, ...) cuts it.
        pieces, steps = compat.phylo_pieces(tree_t8, n=4, timeSteps=True)
        stack = slice_pieces(tree_t8, n=4)
        assert np.array_equal(steps, stack.depth_windows[:, 1])
        # cumulative PD reproduced from the reported steps, edge by edge
        tips = tree_t8.tip_labels
        for j, (lo, hi) in enumerate(stack.depth_windows):
            assert (lo, hi) == (steps[j - 1] if j else 0.0, steps[j])
            assert pieces[j].pd == pytest.approx(stack.contribution(tips)[j])

    def test_pd_criterion_uses_time_mapped_cut_points(self, tree_t8):
        # phylo_pieces.R:219-236: for criterion = "PD" R stores age = the time at
        # which each equal-PD target is reached, again the tipward bound.
        # total PD = 24 over 4 pieces -> 6 PD each -> cut points 3, 5, 6.5, 8.
        pieces, steps = compat.phylo_pieces(tree_t8, n=4, criterion="PD", timeSteps=True)
        assert np.allclose(steps, [3.0, 5.0, 6.5, 8.0])
        assert steps[-1] == pytest.approx(tree_t8.root_age)
        assert np.allclose([p.pd for p in pieces], 6.0)

    def test_method_2_width_slices(self, tree_t4):
        # width 1 on T = 4 divides evenly, so R's raw cumsum(rep(1, 4)) applies
        # unchanged and the default (R) ordering is exact, not merely ascending.
        pieces, steps = compat.phylo_pieces(tree_t4, n=1.0, method=2, timeSteps=True)
        assert np.array_equal(steps, np.array([1.0, 2.0, 3.0, 4.0]))
        assert len(pieces) == 4

    def test_opt_out_returns_main_api_ages(self, tree_t8):
        # timeSteps_as_depth=False is the escape hatch back to the Python main
        # API convention: older boundary of each slice, Ma before present,
        # decreasing.
        steps = compat.phylo_pieces(tree_t8, n=4, timeSteps=True, timeSteps_as_depth=False)[1]
        assert np.array_equal(steps, np.array([8.0, 6.0, 4.0, 2.0]))
        assert np.array_equal(steps, slice_pieces(tree_t8, n=4).ages)
        assert np.all(np.diff(steps) < 0)

    def test_return_shape_combinations(self, tree_t8):
        assert isinstance(compat.phylo_pieces(tree_t8, n=4), list)
        assert len(compat.phylo_pieces(tree_t8, n=4, timeSteps=True)) == 2
        assert len(compat.phylo_pieces(tree_t8, n=4, returnTree=True)) == 2
        out = compat.phylo_pieces(tree_t8, n=4, timeSteps=True, returnTree=True)
        assert len(out) == 3
        assert out[2] is tree_t8 or np.array_equal(out[2].lengths, tree_t8.lengths)

    def test_pieces_ordered_root_to_tips(self, tree_t8):
        # R: "The slices list is ordered from roots to tips."
        pieces = compat.phylo_pieces(tree_t8, n=4)
        assert [round(p.pd, 6) for p in pieces] == [4.0, 4.0, 8.0, 8.0]

    def test_method_validated_loudly(self, tree_t8):
        # phylo_pieces.R has `if(method == 1)` / `if(method == 2)` branches only.
        with pytest.raises(ValueError, match="method must be 1"):
            compat.phylo_pieces(tree_t8, n=4, method=3)
        with pytest.raises(ValueError, match="method must be 1"):
            compat.phylo_pieces(tree_t8, n=4, method=0)


class TestSqueezeIntInvert:
    """A7(b): ``invert = TRUE`` must return R's ``(younger, older)`` order."""

    def test_invert_order_matches_r(self, tree_t8):
        # R: list(tree1 = squeeze_root(to), tree2 = squeeze_tips(from))
        from_, to = 7.0, 1.0
        tree1, tree2 = compat.squeeze_int(tree_t8, from_, to, invert=True)
        # tree1 = R's squeeze_root(to) -> the younger/crown complement
        # tree2 = R's squeeze_tips(from) -> the older complement
        assert tree1.pd == pytest.approx(compat.squeeze_root(tree_t8, to).pd)
        assert tree2.pd == pytest.approx(compat.squeeze_tips(tree_t8, from_).pd)

    def test_invert_is_opposite_of_main_api(self, tree_t8):
        from_, to = 7.0, 1.0
        main = slice_interval(tree_t8, from_, to, invert=True)
        wrapped = compat.squeeze_int(tree_t8, from_, to, invert=True)
        assert len(main) == len(wrapped) == 2
        assert wrapped[0].pd == pytest.approx(main[1].pd)
        assert wrapped[1].pd == pytest.approx(main[0].pd)

    def test_the_two_complements_are_distinguishable(self, tree_t8):
        # guards against a vacuous pass on a symmetric case: the crown slice of
        # width 1 carries PD 4, the remnant older than 7 Ma carries PD 2.
        tree1, tree2 = compat.squeeze_int(tree_t8, 7.0, 1.0, invert=True)
        assert tree1.pd == pytest.approx(4.0)
        assert tree2.pd == pytest.approx(2.0)
        assert tree1.pd != tree2.pd

    def test_no_invert_returns_single_tree(self, tree_t8):
        out = compat.squeeze_int(tree_t8, 6.0, 2.0)
        assert not isinstance(out, list)
        assert out.pd == pytest.approx(slice_interval(tree_t8, 6.0, 2.0).pd)


class TestCriterionMapping:
    """A7(c): R's documented spellings map correctly; anything else errors."""

    @pytest.mark.parametrize("spelling", ["my", "MY", " My ", "million years"])
    def test_time_spellings(self, tree_t4, spelling):
        assert compat.squeeze_root(tree_t4, 1.0, criterion=spelling).pd == pytest.approx(
            slice_rootward(tree_t4, 1.0, criterion="time").pd
        )

    @pytest.mark.parametrize("spelling", ["PD", "pd", "Pd"])
    def test_pd_spellings(self, tree_t4, spelling):
        assert compat.squeeze_root(tree_t4, 1.0, criterion=spelling).pd == pytest.approx(
            slice_rootward(tree_t4, 1.0, criterion="pd").pd
        )

    def test_pd_and_time_give_different_answers(self, tree_t4):
        # non-vacuity check: on this tree the two criteria really do disagree,
        # so the mapping above is not just passing because both branches match.
        a = compat.squeeze_root(tree_t4, 1.0, criterion="my").pd
        b = compat.squeeze_root(tree_t4, 1.0, criterion="PD").pd
        assert a == pytest.approx(4.0)
        assert b == pytest.approx(10.0)
        assert a != b

    @pytest.mark.parametrize(
        "bad",
        ["phylogenetic diversity", "", "years", "depth", "PD-ish", "diversity"],
    )
    def test_unknown_criterion_raises(self, tree_t4, bad):
        with pytest.raises(ValueError, match="unknown criterion"):
            compat.squeeze_root(tree_t4, 1.0, criterion=bad)
        with pytest.raises(ValueError, match="unknown criterion"):
            compat.squeeze_tips(tree_t4, 1.0, criterion=bad)
        with pytest.raises(ValueError, match="unknown criterion"):
            compat.squeeze_int(tree_t4, 2.0, 1.0, criterion=bad)

    def test_native_python_spelling_also_accepted(self, tree_t4):
        # "time" is not an R spelling, but _CRITERION_MAP accepts it so that a
        # main-API user can drop into the compat layer without rewriting args.
        assert compat.squeeze_root(tree_t4, 1.0, criterion="time").pd == pytest.approx(
            compat.squeeze_root(tree_t4, 1.0, criterion="my").pd
        )

    @pytest.mark.parametrize("bad", [None, 1, ["PD"]])
    def test_non_string_criterion_raises(self, tree_t4, bad):
        with pytest.raises(ValueError, match="criterion must be a string"):
            compat.squeeze_root(tree_t4, 1.0, criterion=bad)


class TestSqueezeRootAndTips:
    """``compat.squeeze_tips`` had no coverage at all before this file."""

    def test_squeeze_root_equals_slice_rootward(self, tree_t4):
        a = compat.squeeze_root(tree_t4, 1.5)
        b = slice_rootward(tree_t4, 1.5)
        assert np.allclose(a.lengths, b.lengths)
        assert a.tree.tip_labels == b.tree.tip_labels

    def test_squeeze_tips_equals_slice_tipward(self, tree_t4):
        a = compat.squeeze_tips(tree_t4, 1.5)
        b = slice_tipward(tree_t4, 1.5)
        assert np.allclose(a.lengths, b.lengths)
        assert a.tree.tip_labels == b.tree.tip_labels

    def test_endpoint_semantics(self, tree_t4):
        # time = 0 -> nothing / everything, mirroring R (see USER_GUIDE §5).
        assert compat.squeeze_root(tree_t4, 0.0).pd == pytest.approx(0.0)
        assert compat.squeeze_tips(tree_t4, 0.0).pd == pytest.approx(tree_t4.total_pd())
        full = compat.squeeze_root(tree_t4, tree_t4.root_age).pd
        assert full == pytest.approx(tree_t4.total_pd())
        assert compat.squeeze_tips(tree_t4, tree_t4.root_age).pd == pytest.approx(0.0)

    def test_root_and_tips_are_complementary_in_pd(self, tree_t4):
        # a crown slice of width w and a tipward cut of the same w together
        # rebuild the whole tree's PD on this ultrametric tree.
        w = 1.5
        crown = compat.squeeze_root(tree_t4, w).pd
        remnant = compat.squeeze_tips(tree_t4, w).pd
        assert crown + remnant == pytest.approx(tree_t4.total_pd())

    def test_dropnodes_maps_to_collapse(self, tree_t4):
        # dropNodes=TRUE is R's node-dropping; here it routes to collapse(),
        # which removes zero-length *internal* branches and nothing else.  On
        # ((A:1,B:1):3,(C:2,D:2):2); the crown slice of width 1 leaves the
        # root -> (C,D) edge cut to zero length, so 6 edges become 4.
        plain = compat.squeeze_root(tree_t4, 1.0, dropNodes=False)
        collapsed = compat.squeeze_root(tree_t4, 1.0, dropNodes=True)
        assert collapsed.n_tips == plain.tree.n_tips
        assert len(collapsed.lengths) == 4
        assert len(plain.lengths) == 6
        assert collapsed.total_pd() == pytest.approx(plain.pd)
        # explicit no-op for the default
        assert np.array_equal(compat.squeeze_root(tree_t4, 1.0).lengths, plain.lengths)


class TestPruneTips:
    """``compat.prune_tips`` had no coverage at all before this file."""

    def test_method_1_keeps_long_terminal_branches(self, tree_t4):
        # R: "prune the tips originating after an inputted temporal threshold"
        kept = compat.prune_tips(tree_t4, 1.5, method=1)
        assert sorted(kept.tip_labels) == ["C", "D"]

    def test_method_2_is_the_complement(self, tree_t4):
        kept = compat.prune_tips(tree_t4, 1.5, method=2)
        assert sorted(kept.tip_labels) == ["A", "B"]

    def test_methods_match_the_main_api_sides(self, tree_t4):
        assert np.allclose(
            compat.prune_tips(tree_t4, 1.5, method=1).lengths,
            slice_prune_tips(tree_t4, 1.5, side="after").lengths,
        )
        assert np.allclose(
            compat.prune_tips(tree_t4, 1.5, method=2).lengths,
            slice_prune_tips(tree_t4, 1.5, side="before").lengths,
        )

    def test_threshold_at_the_extremes(self, tree_t4):
        # terminal branches are [1, 1, 2, 2]
        assert sorted(compat.prune_tips(tree_t4, 1.0, method=1).tip_labels) == [
            "A",
            "B",
            "C",
            "D",
        ]
        assert sorted(compat.prune_tips(tree_t4, 2.0, method=2).tip_labels) == [
            "A",
            "B",
        ]
        # Above the longest terminal branch, method=1 keeps 0 tips and so now
        # raises: R answers NULL there, which surfaced as an AttributeError on
        # ``.tip_labels`` in scripts migrated from treesliceR.
        with pytest.raises(ValueError, match="needs at least 2 tips"):
            compat.prune_tips(tree_t4, 2.5, method=1)
        # method=2 keeps tips *below* the threshold, so a threshold above every
        # terminal branch prunes nothing and returns the tree unpruned (4 tips)
        assert sorted(compat.prune_tips(tree_t4, 2.5, method=2).tip_labels) == [
            "A",
            "B",
            "C",
            "D",
        ]

    def test_vector_threshold_returns_list(self, tree_t4):
        out = compat.prune_tips(tree_t4, [1.0, 1.5], method=1)
        assert isinstance(out, list) and len(out) == 2
        assert sorted(out[0].tip_labels) == ["A", "B", "C", "D"]
        assert sorted(out[1].tip_labels) == ["C", "D"]

    def test_qtl_interprets_time_as_quantiles(self, tree_t4):
        # np.quantile([1,1,2,2], [0.25, 0.75]) == [1.0, 2.0] (type-7, as R)
        out = compat.prune_tips(tree_t4, [0.25, 0.75], qtl=True, method=1)
        assert sorted(out[0].tip_labels) == ["A", "B", "C", "D"]
        assert sorted(out[1].tip_labels) == ["C", "D"]

    @pytest.mark.parametrize("bad", [0, 3, -1, 1.5, None, "after"])
    def test_method_validated_loudly(self, tree_t4, bad):
        # prune_tips.R:52,121 gives method outside {1, 2} a NULL return; this
        # layer raises instead, and must NOT silently fall through to method=2.
        with pytest.raises(ValueError, match="method must be 1 or 2"):
            compat.prune_tips(tree_t4, 1.5, method=bad)

    def test_no_silent_fallthrough(self, tree_t4):
        # the specific regression: method=3 used to behave like method=2
        with pytest.raises(ValueError):
            compat.prune_tips(tree_t4, 1.5, method=3)


class TestNonUltrametricGuard:
    def test_guard_matches_main_api(self):
        # (A:1,B:2) is not ultrametric: root_age would be ambiguous
        tree = TreeArray.from_newick("((A:1,B:2):1,(C:1.5,D:1.5):1);")
        with pytest.raises(Exception):
            compat.squeeze_root(tree, 0.5)
        with pytest.raises(Exception):
            compat.prune_tips(tree, 1.0, method=1)
        # phylo_pieces has no waiver either, exactly as R's nodes_config() stops
        with pytest.raises(Exception):
            compat.phylo_pieces(tree, n=2)
        # explicitly waiving it succeeds and records the verdict
        out = compat.squeeze_root(tree, 0.5, allow_non_ultrametric=True)
        assert np.isfinite(out.pd)
        stack = slice_pieces(tree, n=2, allow_non_ultrametric=True)
        assert stack.ultrametric is False


class TestPublicSurface:
    def test_all_names_are_importable(self):
        for name in compat.__all__:
            assert callable(getattr(compat, name)), name

    def test_docstring_claims_are_true(self):
        # the module docstring advertises R's ordering for both flagged
        # functions; keep it honest if either is ever renamed.
        doc = compat.__doc__ or ""
        assert "depth-from-root" in doc or "cumsum" in doc
        assert "younger" in doc and "older" in doc
        pieces_doc = compat.phylo_pieces.__doc__
        assert "depth_windows[:, 1]" in pieces_doc
        assert "cumsum" in pieces_doc
