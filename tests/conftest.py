"""Shared fixtures: hand-built trees with analytic properties."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phyloslicer.core import TreeArray
from phyloslicer.simulate import yule_tree


def balanced_4() -> TreeArray:
    """((A:1,B:1):1,(C:1,D:1):1) - root age 2, total PD 6, perfectly symmetric."""
    # nodes: tips 0..3 = A,B,C,D; root = 4; internal 5 = (A,B), 6 = (C,D)
    edges = np.array(
        [
            [4, 5],
            [4, 6],  # root edges, depth 0-1
            [5, 0],
            [5, 1],  # depth 1-2
            [6, 2],
            [6, 3],
        ]
    )
    lengths = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    return TreeArray(tip_labels=["A", "B", "C", "D"], edges=edges, lengths=lengths)


def caterpillar_5() -> TreeArray:
    """Ultrametric 5-tip caterpillar: ((((A:1,B:1):1,C:2):1,D:3):1,E:4) from the root.

    Root age 4 (all tips at depth 4), total PD 14, fully asymmetric.
    """
    # tips 0..4 = A..E; root = 5; internals 6 = (A,B,C,D), 7 = (A,B,C), 8 = (A,B)
    edges = np.array(
        [
            [5, 6],  # root -> n6, len 1 (depth 0-1)
            [5, 4],  # root -> E, len 4 (depth 0-4)
            [6, 7],  # n6 -> n7, len 1 (depth 1-2)
            [6, 3],  # n6 -> D, len 3 (depth 1-4)
            [7, 8],  # n7 -> n8, len 1 (depth 2-3)
            [7, 2],  # n7 -> C, len 2 (depth 2-4)
            [8, 0],  # n8 -> A, len 1 (depth 3-4)
            [8, 1],  # n8 -> B, len 1 (depth 3-4)
        ]
    )
    lengths = np.array([1.0, 4.0, 1.0, 3.0, 1.0, 2.0, 1.0, 1.0])
    return TreeArray(tip_labels=["A", "B", "C", "D", "E"], edges=edges, lengths=lengths)


@pytest.fixture
def balanced4() -> TreeArray:
    return balanced_4()


@pytest.fixture
def caterpillar5() -> TreeArray:
    return caterpillar_5()


@pytest.fixture
def random_tree_20() -> TreeArray:
    return yule_tree(20, age=10.0, seed=123)


@pytest.fixture
def random_mat_20(random_tree_20: TreeArray) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    mat = (rng.random((12, random_tree_20.n_tips)) < 0.3).astype(float)
    mat_df = pd.DataFrame(mat, columns=random_tree_20.tip_labels)
    mat_df.index = [f"s{i + 1}" for i in range(mat.shape[0])]
    # guarantee every species occurs at least once
    empty = mat_df.sum(axis=0) == 0
    if empty.any():
        cols = mat_df.columns[empty.to_numpy()]
        mat_df[cols] = 1.0
    return mat_df


@pytest.fixture
def span_x_extent():
    """Return ``(x0, x1)`` of the shading the legend carries under a label.

    ``Axes.axvspan`` hands back a ``Polygon`` before matplotlib 3.10 and a
    ``Rectangle`` from 3.10 on, and a legend handle can be a ``PolyCollection``
    wrapper around either, so the extent is read from whichever vertex accessor
    the artist offers.  Reading it through the legend, rather than positionally
    through ``ax.patches`` or ``ax.collections``, is what makes the assertion
    test the shading the reader actually sees.
    """

    def _extent(ax, label: str) -> tuple[float, float]:
        handles, labels = ax.get_legend_handles_labels()
        matched = [h for h, lab in zip(handles, labels, strict=True) if lab == label]
        assert len(matched) == 1, f"expected one legend entry {label!r}, got {labels}"
        artist = matched[0]
        from matplotlib.patches import Rectangle

        if isinstance(artist, Rectangle):
            return float(artist.get_x()), float(artist.get_x() + artist.get_width())
        if hasattr(artist, "get_paths"):
            paths = artist.get_paths()
            assert paths, f"{type(artist).__name__} carries no path"
            xy = paths[0].vertices
        elif hasattr(artist, "get_xy"):
            xy = artist.get_xy()
        else:
            xy = artist.get_path().vertices
        xs = np.asarray(xy, dtype=float)[:, 0]
        return float(xs.min()), float(xs.max())

    return _extent


@pytest.fixture
def register_reference():
    """Write a provenance manifest that registers one reference CSV.

    Mirrors ``validation/reference/manifest.txt`` for the shipped scenarios: the
    CLI ``validate`` command and ``validation/compare.py`` refuse any reference
    that is not registered by SHA256, shape, column pattern and orientation
    (an unregistered reference, and a self-supplied reference), so tests that
    exercise them have to register the
    files they invent.  The manifest is written *after* the CSV exists; a test
    that needs a mismatch is free to edit either side afterwards.
    """
    from phyloslicer.cli.main import format_reference_manifest, make_reference_record

    def _register(
        path,
        *,
        rate_column: str = "CpD",
        orientation: str = "site-by-index",
        columns_pattern: str = r"^.+$",
        manifest_name: str = "manifest.txt",
        **extra: str,
    ) -> Path:
        path = Path(path)
        record = make_reference_record(
            path,
            orientation=orientation,
            columns_pattern=columns_pattern,
            rate_column=rate_column,
            **extra,
        )
        text = format_reference_manifest(
            {"note": "provenance manifest written by the register_reference fixture"},
            {path.name: record},
        )
        manifest = path.parent / manifest_name
        manifest.write_text(text, encoding="utf-8")
        return manifest

    return _register
