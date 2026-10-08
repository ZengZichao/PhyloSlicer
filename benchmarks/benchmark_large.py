#!/usr/bin/env python
"""Python side of the large-input (1000- and 2000-tip, 250-assemblage)
cross-language benchmark, on the byte-identical inputs written by
``make_large_inputs.py`` and timed by ``benchmark_large_r.R``.

Reported times are the median of five repeats, the protocol documented in
``benchmarks/README.md`` (an earlier revision of this file recorded the fastest
of five repeats, which is not what that protocol says).  ``reps`` says so in the
table itself, and ``source`` names the input tree each row was measured on.

Usage (from the project root):
  python benchmarks/benchmark_large.py
Writes benchmarks/results_py_large.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from phyloslicer.io.matrices import align_tree_matrix, as_matrix
from phyloslicer.io.newick import read_newick
from phyloslicer.rates.api import cpd_rate
from phyloslicer.slicing import slice_pieces

from machine_info import describe, machine_fields  # noqa: E402

REPS = 5
N_SLICES = 100

#: NV-08: read once, before any timing loop, and merged into every row this driver writes.
MACHINE = machine_fields()


def median_time(fn, reps: int = REPS) -> float:
    ts = timed_runs(fn, reps)
    return float(np.median(ts))


def timed_runs(fn, reps: int = REPS) -> list[float]:
    """Every repeat, un-summarised."""
    ts = []
    fn()  # warm-up
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return ts


def main(
    in_dir: str = "benchmarks/large", out_csv: str = "benchmarks/results_py_large.csv"
) -> None:
    print(describe(MACHINE))
    root = Path(in_dir)
    rows = []
    for tp in sorted(root.glob("*__tree.tre")):
        tag = tp.name.replace("__tree.tre", "")
        mp = root / f"{tag}__mat.csv"
        if not mp.exists():
            continue
        tree = read_newick(tp)
        mat = pd.read_csv(mp, index_col=0)
        mat.columns = mat.columns.astype(str)
        mat = (mat > 0).astype(float)
        tree_a, mat_a = align_tree_matrix(tree, as_matrix(mat))
        ts_slice = timed_runs(lambda: slice_pieces(tree_a, n=N_SLICES).contribution_matrix())
        ts_cpd = timed_runs(lambda: cpd_rate(tree_a, mat_a, n_slices=N_SLICES))
        rows.append(
            {
                "tips": tree_a.n_tips,
                "n_sites": mat_a.shape[0],
                "n_slices": N_SLICES,
                "slice_s": float(np.median(ts_slice)),
                "cpd_s": float(np.median(ts_cpd)),
                "source": tag,
                "reps": REPS,
                "repeat_slice_s": ";".join(f"{x:.9g}" for x in ts_slice),
                "repeat_cpd_s": ";".join(f"{x:.9g}" for x in ts_cpd),
                **MACHINE,
            }
        )
        print(
            f"{tag} tips={tree_a.n_tips} sites={mat_a.shape[0]} k={N_SLICES}: "
            f"slice {rows[-1]['slice_s']:.6f}s cpd {rows[-1]['cpd_s']:.4f}s"
        )
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
