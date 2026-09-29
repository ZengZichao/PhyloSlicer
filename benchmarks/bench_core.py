#!/usr/bin/env python
"""Micro-benchmarks for the PhyloSlicer kernels.

Times the two operations whose cost determines practical scale:
  1. slice_pieces (interval-intersection matrix) vs the R-style approach of
     materialising k sliced tree copies;
  2. cpd_rate (sparse-matrix pipeline) vs a naive per-site, per-slice loop.

Usage: python benchmarks/bench_core.py
Writes benchmarks/results.csv and prints a summary table.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from phyloslicer.rates.api import cpd_rate
from phyloslicer.simulate import make_dataset, yule_tree
from phyloslicer.slicing import slice_pieces


def naive_slice_pieces(tree, k):
    """R-style reference: materialise k sliced tree copies (O(k * E))."""
    C = slice_pieces(tree, n=k).contribution_matrix()
    out = []
    for j in range(k):
        out.append(tree.with_lengths(C[:, j]))  # full tree copy per slice
    return out


def naive_cpd(tree, mat, k):
    """Per-site, per-slice Python loops (O(n_sites * k * S))."""
    slices = slice_pieces(tree, n=k)
    C = slices.contribution_matrix()
    intervals = tree.tip_intervals()
    lo, hi = intervals[:, 0], intervals[:, 1]
    from phyloslicer.rates.fitting import fit_rate

    rates = []
    values = mat.to_numpy() > 0
    for s in range(values.shape[0]):
        tips = np.flatnonzero(values[s])
        mask = np.zeros(tree.n_edges, dtype=bool)
        for t in tips:
            mask |= (lo <= t) & (t < hi)
        prof = np.array([C[mask, j].sum() for j in range(k)])
        fit = fit_rate(prof, slices.ages)
        rates.append(fit.rate)
    return rates


def bench(fn, *args, repeat=3):
    best = np.inf
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn(*args)
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    rows = []
    for n_tips in (100, 1_000):
        for n_sites in (100, 1_000):
            tree = yule_tree(n_tips, age=1.0, seed=1)
            mat, _ = make_dataset(n_sites=n_sites, n_species=n_tips, structure="clustered", seed=2)
            mat.columns = tree.tip_labels
            mat.index = [f"s{i + 1}" for i in range(n_sites)]
            k = 100

            t_slice = bench(lambda: slice_pieces(tree, n=k).contribution_matrix())
            t_slice_naive = bench(lambda: naive_slice_pieces(tree, k), repeat=1)
            t_cpd = bench(lambda: cpd_rate(tree, mat, n_slices=k))
            t_cpd_naive = bench(lambda: naive_cpd(tree, mat, k), repeat=1)

            rows.append(
                {
                    "n_tips": n_tips,
                    "n_sites": n_sites,
                    "n_slices": k,
                    "slice_vectorized_s": t_slice,
                    "slice_naive_s": t_slice_naive,
                    "slice_speedup": t_slice_naive / t_slice,
                    "cpd_vectorized_s": t_cpd,
                    "cpd_naive_s": t_cpd_naive,
                    "cpd_speedup": t_cpd_naive / t_cpd,
                }
            )
            print(
                f"tips={n_tips:>6} sites={n_sites:>5}: slice {t_slice_naive / t_slice:8.1f}x, "
                f"cpd {t_cpd_naive / t_cpd:8.1f}x"
            )

    out = Path(__file__).parent / "results.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"results -> {out}")


if __name__ == "__main__":
    main()
