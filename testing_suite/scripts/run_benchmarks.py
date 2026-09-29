#!/usr/bin/env python
"""Benchmark PhyloSlicer on various data sizes.

Measures runtime and memory for the core operations across different
tree sizes and slice counts. Results are saved as CSV for citation.
"""

from __future__ import annotations

import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))

import phyloslicer as ps
from phyloslicer.simulate import make_dataset, yule_tree

BENCH_DIR = Path(__file__).resolve().parent.parent / "benchmarks"
BENCH_DIR.mkdir(parents=True, exist_ok=True)


def bench_one(n_tips, n_sites, n_slices, seed=42):
    """Benchmark a single configuration."""
    tree = yule_tree(n_tips, age=10.0, seed=seed)
    mat, coords = make_dataset(n_sites=n_sites, n_species=n_tips, structure="clustered", seed=seed)
    mat.columns = tree.tip_labels

    results = {}

    # Slicing
    tracemalloc.start()
    t0 = time.perf_counter()
    stack = ps.slice_pieces(tree, n=n_slices)
    t_slice = time.perf_counter() - t0
    _, peak_slice = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results["slice_s"] = t_slice
    results["slice_MB"] = peak_slice / 1e6

    # PD per slice
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree2, mat2 = ps.io.align_tree_matrix(tree, mat)

    tracemalloc.start()
    t0 = time.perf_counter()
    from phyloslicer.indices.api import pd_per_slice

    pd_per_slice(tree2, mat2, stack)
    t_pdperslice = time.perf_counter() - t0
    _, peak_pdperslice = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results["pd_per_slice_s"] = t_pdperslice
    results["pd_per_slice_MB"] = peak_pdperslice / 1e6

    # CpD rate
    tracemalloc.start()
    t0 = time.perf_counter()
    res = ps.cpd_rate(tree2, mat2, n_slices=n_slices)
    t_cpd = time.perf_counter() - t0
    _, peak_cpd = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results["cpd_rate_s"] = t_cpd
    results["cpd_rate_MB"] = peak_cpd / 1e6
    results["n_converged"] = int(res["converged"].sum())

    # Posterior (if trees available)
    trees = [yule_tree(n_tips, age=10.0, seed=seed + i) for i in range(3)]
    tracemalloc.start()
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ps.posterior_rate(trees, mat2, n_slices=n_slices)
    t_post = time.perf_counter() - t0
    _, peak_post = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    results["posterior_s"] = t_post
    results["posterior_MB"] = peak_post / 1e6

    return results


def main():
    print("=" * 70)
    print("PhyloSlicer Benchmark Suite")
    print("=" * 70)
    print(f"Python: {sys.executable}")
    print(f"Version: {ps.__version__}")
    print()

    configs = [
        # (n_tips, n_sites, n_slices)
        (20, 10, 50),
        (50, 20, 50),
        (100, 30, 50),
        (100, 30, 100),
        (100, 30, 200),
        (200, 50, 100),
        (500, 100, 100),
    ]

    all_results = []
    for n_tips, n_sites, n_slices in configs:
        print(f"Benchmarking: {n_tips} tips, {n_sites} sites, {n_slices} slices...")
        # Run 3 times, take median
        runs = []
        for _ in range(3):
            r = bench_one(n_tips, n_sites, n_slices)
            runs.append(r)

        # Median
        med = {}
        for key in runs[0]:
            vals = [r[key] for r in runs]
            med[key] = np.median(vals) if isinstance(vals[0], float) else vals[0]

        row = {
            "n_tips": n_tips,
            "n_sites": n_sites,
            "n_slices": n_slices,
            **med,
        }
        all_results.append(row)
        print(
            f"  slice: {med['slice_s']:.4f}s, cpd_rate: {med['cpd_rate_s']:.4f}s, "
            f"posterior: {med['posterior_s']:.4f}s"
        )

    df = pd.DataFrame(all_results)
    out_path = BENCH_DIR / "benchmark_results.csv"
    df.to_csv(out_path, index=False)
    print(f"\nResults saved to: {out_path}")
    print("\nBenchmark summary:")
    print(df.to_string(index=False))

    # Also write a markdown summary
    md_path = BENCH_DIR / "BENCHMARK_REPORT.md"
    with open(md_path, "w") as f:
        f.write("# PhyloSlicer Benchmark Report\n\n")
        f.write(f"**Version:** {ps.__version__}  \n")
        f.write(f"**Date:** {time.strftime('%Y-%m-%d')}  \n")
        f.write(f"**Python:** {sys.version.split()[0]}  \n\n")
        f.write("## Methodology\n\n")
        f.write("Each configuration is run 3 times; the table reports the median.\n")
        f.write("Memory is peak traced memory (tracemalloc).\n\n")
        f.write("## Results\n\n")
        f.write(
            "| n_tips | n_sites | n_slices | slice (s) | pd_per_slice (s) "
            "| cpd_rate (s) | cpd_MB | posterior (s) | n_converged |\n"
        )
        f.write(
            "|--------|---------|----------|-----------|------------------|--------------|--------|---------------|-------------|\n"
        )
        for _, row in df.iterrows():
            f.write(
                f"| {int(row['n_tips'])} | {int(row['n_sites'])} | {int(row['n_slices'])} "
                f"| {row['slice_s']:.4f} | {row['pd_per_slice_s']:.4f} "
                f"| {row['cpd_rate_s']:.4f} | {row['cpd_rate_MB']:.1f} "
                f"| {row['posterior_s']:.4f} | {int(row['n_converged'])} |\n"
            )
        f.write("\n## Scaling\n\n")
        f.write("The contribution matrix is `E x k` (edges x slices), so:\n")
        f.write("- Slicing cost is O(E * k) — linear in both tree size and slice count.\n")
        f.write("- PD per slice is O(S * E * k) via sparse matrix product.\n")
        f.write("- Rate fitting is vectorised: all sites solved simultaneously.\n")
    print(f"Benchmark report saved to: {md_path}")


if __name__ == "__main__":
    main()
