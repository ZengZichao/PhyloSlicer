"""Core TreeArray behaviour."""

from __future__ import annotations

import numpy as np
import pytest

from phyloslicer.core import TreeArray
from phyloslicer.simulate import birth_death_tree, yule_tree


class TestBasicProperties:
    def test_balanced4_structure(self, balanced4):
        assert balanced4.n_tips == 4
        assert balanced4.n_edges == 6
        assert balanced4.root == 4
        assert balanced4.root_age == pytest.approx(2.0)
        assert balanced4.total_pd() == pytest.approx(6.0)
        assert balanced4.is_ultrametric

    def test_caterpillar5_structure(self, caterpillar5):
        assert caterpillar5.root_age == pytest.approx(4.0)
        assert caterpillar5.total_pd() == pytest.approx(14.0)
        assert caterpillar5.is_ultrametric

    def test_node_depths(self, caterpillar5):
        d = caterpillar5.node_depths()
        assert d[caterpillar5.root] == 0.0
        assert d[:5].min() == pytest.approx(4.0)  # all tips at the root age

    def test_non_ultrametric_flag(self):
        edges = np.array([[2, 0], [2, 1]])
        tree = TreeArray(tip_labels=["A", "B"], edges=edges, lengths=[1.0, 2.0])
        assert not tree.is_ultrametric

    def test_validation(self):
        with pytest.raises(ValueError):
            TreeArray(tip_labels=["A"], edges=np.array([[1, 0]]), lengths=[1.0])
        with pytest.raises(ValueError):
            TreeArray(tip_labels=["A", "B"], edges=np.array([[2, 0], [2, 1]]), lengths=[-1.0, 1.0])

    def test_invalid_node_numbering_rejected(self):
        # tip 1 never appears as a child (id 1 used as internal instead)
        with pytest.raises(ValueError, match="numbering"):
            TreeArray(
                tip_labels=["A", "B", "C"],
                edges=np.array([[3, 0], [3, 2], [2, 4]]),
                lengths=[1.0, 1.0, 1.0],
            )
        # root id n_tips must not have a parent edge
        with pytest.raises(ValueError, match="root"):
            TreeArray(
                tip_labels=["A", "B", "C"],
                edges=np.array([[3, 2], [3, 0], [2, 1], [1, 3]]),
                lengths=[1.0, 1.0, 1.0, 1.0],
            )
        # duplicate incoming edges for one tip
        with pytest.raises(ValueError, match="incoming"):
            TreeArray(
                tip_labels=["A", "B", "C"],
                edges=np.array([[3, 0], [3, 1], [3, 2], [2, 2]]),
                lengths=[1.0, 1.0, 1.0, 1.0],
            )
        # disconnected component
        with pytest.raises(ValueError, match="unreachable"):
            TreeArray(
                tip_labels=["A", "B", "C"],
                edges=np.array([[3, 0], [3, 1], [5, 2]]),
                lengths=[1.0, 1.0, 1.0],
            )


class TestTipIntervals:
    def test_intervals_partition_tips(self, random_tree_20):
        intervals = random_tree_20.tip_intervals()
        # root edges together must cover all tips
        root = random_tree_20.root
        root_mask = random_tree_20.edges[:, 0] == root
        assert intervals[root_mask, 0].min() == 0
        assert intervals[root_mask, 1].max() == random_tree_20.n_tips
        # tip edges are singletons
        tip_edges = random_tree_20.edges[:, 1] < random_tree_20.n_tips
        assert (np.diff(intervals[tip_edges], axis=1) == 1).all()

    def test_membership_via_intervals(self, random_tree_20):
        """Path membership computed from intervals == explicit ancestor walk."""
        tree = random_tree_20
        intervals = tree.tip_intervals()
        parent = tree.parent_of()
        for tip in (0, 7, 19):
            path_edges = set()
            node = tip
            while node in parent:
                for e, (p, c) in enumerate(tree.edges.tolist()):
                    if c == node:
                        path_edges.add(e)
                        node = p
                        break
            edge_arr = np.fromiter(path_edges, dtype=np.int64)
            lo, hi = intervals[edge_arr, 0], intervals[edge_arr, 1]
            assert (lo <= tip).all() and (tip < hi).all()


class TestSubarray:
    def test_prune_keeps_pd_relations(self, caterpillar5):
        sub = caterpillar5.subarray(["A", "B", "C"])
        assert sub.n_tips == 3
        # MRCA of A,B,C is node 7; the root→MRCA stem ([5→6], [6→7]) is
        # discarded (ape::keep.tip, root.edge=0).  Remaining edges:
        # [7→8]=1, [8→A]=1, [8→B]=1, [7→C]=2 → PD = 5
        assert sub.total_pd() == pytest.approx(1 + 1 + 1 + 2)
        assert sub.tip_labels == ["A", "B", "C"]
        assert sub.is_ultrametric

    def test_prune_pair(self, caterpillar5):
        sub = caterpillar5.subarray(["D", "E"])
        # MRCA is the root itself → no stem to discard
        assert sub.total_pd() == pytest.approx(1 + 3 + 4)  # 5→6, 6→D, 5→E


class TestSimulation:
    def test_yule_ultrametric(self):
        tree = yule_tree(30, age=5.0, seed=1)
        assert tree.n_tips == 30
        assert tree.is_ultrametric
        assert tree.root_age == pytest.approx(5.0)

    def test_birth_death_conditioned(self):
        tree = birth_death_tree(12, age=1.0, birth_rate=2.0, death_rate=0.5, seed=3)
        assert tree.n_tips == 12
        assert tree.is_ultrametric
        assert tree.root_age == pytest.approx(1.0)

    def test_reproducible_with_seed(self):
        t1 = yule_tree(15, seed=11)
        t2 = yule_tree(15, seed=11)
        assert (t1.edges == t2.edges).all()
        assert np.allclose(t1.lengths, t2.lengths)
