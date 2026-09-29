#!/usr/bin/env python
"""Fixed-k scaling benchmark, Python side: slicing-kernel time for the same
number of slices across tree sizes, on the byte-identical trees used by
``benchmark_scaling_r.R``.

Usage (from the project root):
  python benchmarks/benchmark_scaling.py
Writes benchmarks/results_scaling_py.csv
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

from phyloslicer.io.newick import read_newick
from phyloslicer.slicing import slice_pieces

REPS = 15  # median of REPS timed runs (sub-millisecond calls need the extra reps)
KS = (25, 50, 100)
DIRS = ("validation/reference", "benchmarks/large", "benchmarks/scaling")


def timed_runs(fn, reps: int = REPS) -> list[float]:
    """Every repeat, un-summarised.

    The archived tables used to store only median and min, so a documented ratio
    could not be re-derived, nor its spread re-checked, from the archived CSVs.
    Callers now write the whole vector.
    """
    ts = []
    fn()  # warm-up (JIT-free, but fills import caches)
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return ts


def median_time(fn, reps: int = REPS) -> tuple[float, float]:
    ts = timed_runs(fn, reps)
    return float(np.median(ts)), float(min(ts))


def main(out_csv: str = "benchmarks/results_scaling_py.csv") -> None:
    paths: list[Path] = []
    for d in DIRS:
        root = Path(d)
        if root.is_dir():
            paths += sorted(root.glob("*__tree.tre"))
    rows = []
    for tp in sorted(set(paths), key=lambda p: (p.read_text().count(","), p.name)):
        tree = read_newick(tp)
        for k in KS:
            ts = timed_runs(lambda: slice_pieces(tree, n=k).contribution_matrix())
            med, mn = float(np.median(ts)), float(min(ts))
            rows.append(
                {
                    "tips": tree.n_tips,
                    "n_slices": k,
                    "source": tp.stem,
                    "slice_s": med,
                    "slice_min": mn,
                    "reps": REPS,
                    "repeat_times_s": ";".join(f"{x:.9g}" for x in ts),
                }
            )
            print(
                f"{tp.parent.name}/{tp.stem} tips={tree.n_tips} k={k}: "
                f"slice median {med:.6f}s (min {mn:.6f}s)"
            )
    df = pd.DataFrame(rows).drop_duplicates(subset=["tips", "n_slices"])
    df.to_csv(out_csv, index=False)
    print("wrote", out_csv, len(df), "rows")


if __name__ == "__main__":
    main()
