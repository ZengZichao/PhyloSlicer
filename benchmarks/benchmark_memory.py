"""Python-side peak memory: the single (E, k) array vs the copy-per-slice design.

    python benchmarks/benchmark_memory.py                # 2000 tips, k = 25/50/100
    python benchmarks/benchmark_memory.py --tips 500

Writes ``benchmarks/results_memory_py.csv`` with ``tracemalloc`` peaks for

* ``array_peak_bytes``  - PhyloSlicer's contribution matrix for all slices at once;
* ``copies_peak_bytes`` - the reference design's structural cost, materialised
  here as ``k`` independent copies of the tree's edge and length arrays.

Scope, stated plainly: this measures the **Python** side only.  The comparable
peak for treesliceR inside its own R process is not measured here, so the
documentation reports the memory advantage as structural rather than as a
cross-language measurement.  Run from the repository root.
"""

from __future__ import annotations

import argparse
import csv
import gc
import sys
import tracemalloc
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phyloslicer import slice_pieces  # noqa: E402
from phyloslicer.core import TreeArray  # noqa: E402
from phyloslicer.simulate import yule_tree  # noqa: E402


def peak_bytes(fn, *args) -> int:
    """Peak traced allocation while ``fn`` runs, in bytes."""
    gc.collect()
    tracemalloc.start()
    try:
        fn(*args)
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def build_array(tree: TreeArray, k: int) -> None:
    """PhyloSlicer's path: one (E, k) contribution matrix for every window."""
    slice_pieces(tree, n=k, criterion="time").contribution_matrix()


def build_copies(tree: TreeArray, k: int) -> None:
    """The reference design's cost: a whole modified tree per slice."""
    kept = [(tree.edges.copy(), tree.lengths.copy(), list(tree.tip_labels)) for _ in range(k)]
    if len(kept) != k:  # pragma: no cover - keeps the copies live for the tracer
        raise AssertionError("copies were collected during measurement")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tips", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--grid", default="25,50,100")
    ap.add_argument("--out", default=str(Path(__file__).parent / "results_memory_py.csv"))
    args = ap.parse_args(argv)

    tree = yule_tree(n=args.tips, age=1.0, seed=args.seed)
    rows = []
    for k in [int(x) for x in args.grid.split(",")]:
        array_peak = peak_bytes(build_array, tree, k)
        copies_peak = peak_bytes(build_copies, tree, k)
        ratio = round(copies_peak / array_peak, 2) if array_peak else ""
        rows.append(
            {
                "tips": args.tips,
                "n_slices": k,
                "array_peak_bytes": array_peak,
                "copies_peak_bytes": copies_peak,
                "copies_over_array": ratio,
                "source": "benchmarks/benchmark_memory.py",
            }
        )
        print(
            f"k={k:4d}  array {array_peak / 1e6:8.3f} MB  copies {copies_peak / 1e6:8.3f} MB"
            f"  copies/array {ratio}"
        )

    with open(args.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {args.out}")
    print("note: Python-side peaks only; the treesliceR/R peak is not measured here.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
