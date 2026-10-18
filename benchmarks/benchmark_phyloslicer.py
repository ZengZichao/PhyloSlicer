#!/usr/bin/env python
"""Time phyloslicer on the same validation inputs used for the R benchmark
(benchmarks/benchmark_r.R) so the cross-language comparison refers to
byte-identical trees and matrices.

Every row of the output names the tree it was measured on (``source``) and the
matrix it was measured with (``structure``): the grid is 3 trees x 2 matrix
structures x 2 slice counts, and only ``source`` + ``structure`` make a row
unique, because ``bd100`` and ``coal100`` both have 100 tips.  ``reps`` records
how many timed runs each cell is the median of, and ``repeat_slice_s`` /
``repeat_cpd_s`` / ``repeat_cpe_s`` keep those runs.

Usage: python benchmarks/benchmark_phyloslicer.py
Writes benchmarks/results_py.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from machine_info import describe, machine_fields  # noqa: E402

from phyloslicer.io.matrices import align_tree_matrix, as_matrix
from phyloslicer.io.newick import read_newick
from phyloslicer.rates.api import cpd_rate, cpe_rate
from phyloslicer.slicing import slice_pieces

REPS = 5  # median of REPS timed runs, the protocol benchmarks/README.md documents

#: NV-08: read once, before any timing loop, and merged into every row this driver writes.
MACHINE = machine_fields()


def _vector(ts: list[float]) -> str:
    """The repeat vector as one CSV cell, so a median can be recomputed from it."""
    return ";".join(f"{x:.9g}" for x in ts)


def timed_runs(fn, *args, repeat: int = REPS) -> list[float]:
    """Every repeat, un-summarised, so a median can be recomputed from its inputs."""
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(*args)
        ts.append(time.perf_counter() - t0)
    return ts


def bench(fn, *args, repeat=REPS):
    """Median of `repeat` timed runs, matching the R side (benchmarks/README.md)."""
    return float(np.median(timed_runs(fn, *args, repeat=repeat)))


def main(ref_dir: str = "validation/reference", out_csv: str = "benchmarks/results_py.csv") -> None:
    print(describe(MACHINE))
    ref_dir = Path(ref_dir)
    rows = []
    for tree_path in sorted(ref_dir.glob("*__tree.tre")):
        tag = tree_path.name.replace("__tree.tre", "")
        tree = read_newick(tree_path)
        for mat_path in sorted(ref_dir.glob(f"{tag}__*__mat.csv")):
            structure = mat_path.name.removesuffix("__mat.csv").removeprefix(f"{tag}__")
            mat = pd.read_csv(mat_path, index_col=0)
            mat.columns = mat.columns.astype(str)
            mat = (mat > 0).astype(float)
            tree_a, mat_a = align_tree_matrix(tree, as_matrix(mat))
            for k in (25, 50):
                ts_slice = timed_runs(lambda: slice_pieces(tree_a, n=k).contribution_matrix())
                ts_cpd = timed_runs(lambda: cpd_rate(tree_a, mat_a, n_slices=k))
                ts_cpe = timed_runs(lambda: cpe_rate(tree_a, mat_a, n_slices=k))
                rows.append(
                    {
                        "tips": tree_a.n_tips,
                        "n_sites": mat_a.shape[0],
                        "n_slices": k,
                        "slice_s": float(np.median(ts_slice)),
                        "cpd_s": float(np.median(ts_cpd)),
                        "cpe_s": float(np.median(ts_cpe)),
                        "source": tag,
                        "structure": structure,
                        "reps": REPS,
                        "repeat_slice_s": _vector(ts_slice),
                        "repeat_cpd_s": _vector(ts_cpd),
                        "repeat_cpe_s": _vector(ts_cpe),
                        **MACHINE,
                    }
                )
                print(
                    f"{tag}/{structure} tips={tree_a.n_tips} sites={mat_a.shape[0]} k={k}: "
                    f"slice {rows[-1]['slice_s']:.4f}s cpd {rows[-1]['cpd_s']:.4f}s "
                    f"cpe {rows[-1]['cpe_s']:.4f}s"
                )
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
