#!/usr/bin/env python
"""Python side of the large-input (1000- and 2000-tip, 250-assemblage)
cross-language benchmark, on the byte-identical inputs written by
``make_large_inputs.py`` and timed by ``benchmark_large_r.R``.

Reported times are the median of five repeats, the protocol documented in
``benchmarks/README.md`` (an earlier revision of this file recorded the fastest
of five repeats, which is not what that protocol says).

Usage (from the project root):
  python benchmarks/benchmark_large.py
Writes benchmarks/results_py_large.csv
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from phyloslicer.io.matrices import align_tree_matrix, as_matrix
from phyloslicer.io.newick import read_newick
from phyloslicer.rates.api import cpd_rate
from phyloslicer.slicing import slice_pieces

REPS = 5
N_SLICES = 100


def median_time(fn, reps: int = REPS) -> float:
    ts = []
    fn()  # warm-up
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def main(
    in_dir: str = "benchmarks/large", out_csv: str = "benchmarks/results_py_large.csv"
) -> None:
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
        t_slice = median_time(lambda: slice_pieces(tree_a, n=N_SLICES).contribution_matrix())
        t_cpd = median_time(lambda: cpd_rate(tree_a, mat_a, n_slices=N_SLICES))
        rows.append(
            {
                "tips": tree_a.n_tips,
                "n_sites": mat_a.shape[0],
                "n_slices": N_SLICES,
                "slice_s": t_slice,
                "cpd_s": t_cpd,
            }
        )
        print(
            f"{tag} tips={tree_a.n_tips} sites={mat_a.shape[0]} k={N_SLICES}: "
            f"slice {t_slice:.6f}s cpd {t_cpd:.4f}s"
        )
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
