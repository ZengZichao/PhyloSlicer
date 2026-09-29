#!/usr/bin/env python
"""Time phyloslicer on the same validation inputs used for the R benchmark
(benchmarks/benchmark_r.R) so the cross-language comparison refers to
byte-identical trees and matrices.

Usage: python benchmarks/benchmark_phyloslicer.py
Writes benchmarks/results_py.csv
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from phyloslicer.io.matrices import align_tree_matrix, as_matrix
from phyloslicer.io.newick import read_newick
from phyloslicer.rates.api import cpd_rate, cpe_rate
from phyloslicer.slicing import slice_pieces


def bench(fn, *args, repeat=5):
    """Median of `repeat` timed runs, matching the R side (benchmarks/README.md)."""
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(*args)
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def main(ref_dir: str = "validation/reference", out_csv: str = "benchmarks/results_py.csv") -> None:
    ref_dir = Path(ref_dir)
    rows = []
    for tree_path in sorted(ref_dir.glob("*__tree.tre")):
        tag = tree_path.name.replace("__tree.tre", "")
        tree = read_newick(tree_path)
        for mat_path in sorted(ref_dir.glob(f"{tag}__*__mat.csv")):
            mat = pd.read_csv(mat_path, index_col=0)
            mat.columns = mat.columns.astype(str)
            mat = (mat > 0).astype(float)
            tree_a, mat_a = align_tree_matrix(tree, as_matrix(mat))
            for k in (25, 50):
                t_slice = bench(lambda: slice_pieces(tree_a, n=k).contribution_matrix())
                t_cpd = bench(lambda: cpd_rate(tree_a, mat_a, n_slices=k))
                t_cpe = bench(lambda: cpe_rate(tree_a, mat_a, n_slices=k))
                rows.append(
                    {
                        "tips": tree_a.n_tips,
                        "n_sites": mat_a.shape[0],
                        "n_slices": k,
                        "slice_s": t_slice,
                        "cpd_s": t_cpd,
                        "cpe_s": t_cpe,
                    }
                )
                print(
                    f"{tag} tips={tree_a.n_tips} sites={mat_a.shape[0]} k={k}: "
                    f"slice {t_slice:.4f}s cpd {t_cpd:.4f}s cpe {t_cpe:.4f}s"
                )
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print("wrote", out_csv)


if __name__ == "__main__":
    main()
