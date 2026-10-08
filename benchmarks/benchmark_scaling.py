#!/usr/bin/env python
"""Fixed-k scaling benchmark, Python side: slicing-kernel time for the same
number of slices across tree sizes, on the byte-identical trees used by
``benchmark_scaling_r.R``.

Usage (from the project root):
  python benchmarks/benchmark_scaling.py
Writes benchmarks/results_scaling_py.csv

Every row carries ``source``, the tag of the ``<tag>__tree.tre`` it was measured
on, spelled exactly as ``benchmark_scaling_r.R`` spells it, because
``recompute_ratios.py`` pairs the two tables on ``(source, tips, n_slices)``.
``tips = 100`` is timed on two trees; the archived table keeps the one the
documentation quotes (see ``PUBLISHED_DUPLICATE``).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from phyloslicer.io.newick import read_newick
from phyloslicer.slicing import slice_pieces

from machine_info import describe, machine_fields  # noqa: E402

REPS = 15  # median of REPS timed runs (sub-millisecond calls need the extra reps)
KS = (25, 50, 100)
DIRS = ("validation/reference", "benchmarks/large", "benchmarks/scaling")

#: NV-08: read once, before any timing loop, and merged into every row this driver writes.
MACHINE = machine_fields()

#: A tree size covered by more than one input tree, and the one the published
#: ladder quotes.  ``tips = 100`` is measured twice - ``bd100`` (birth-death) and
#: ``coal100`` (coalescent) - and the documentation (benchmarks/README.md, caveat
#: 2, and the SI figure) declares ``bd100`` as the pairing on both sides, so the
#: archived table keeps exactly that row and names the discarded scenario in its
#: output instead of letting row order decide.
PUBLISHED_DUPLICATE = {100: "bd100"}


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


def scenario_tag(path: Path) -> str:
    """The input tag of a ``<tag>__tree.tre`` file, spelled as the R driver spells it.

    ``Path.stem`` would archive ``bd100__tree`` while ``benchmark_scaling_r.R``
    archives ``bd100``, and the two tables are paired by this value.
    """
    return path.name.removesuffix("__tree.tre")


def published_rows(rows: list[dict]) -> list[dict]:
    """One row per (tips, n_slices), with a duplicated tree size resolved by name.

    Replaces ``drop_duplicates``, whose outcome depended on the enumeration order:
    a size claimed by two trees keeps the declared published scenario, and a size
    claimed by two trees with no declared preference is an error rather than a coin
    toss.
    """
    by_cell: dict[tuple[int, int], list[dict]] = {}
    for row in rows:
        by_cell.setdefault((row["tips"], row["n_slices"]), []).append(row)
    out: list[dict] = []
    for (tips, _k), group in by_cell.items():
        if len(group) == 1:
            out.append(group[0])
            continue
        preferred = PUBLISHED_DUPLICATE.get(tips)
        if preferred is None:
            raise SystemExit(
                f"tips={tips} was timed on {len(group)} input trees "
                f"{[g['source'] for g in group]} and no published source is declared "
                "for that tree size; add it to PUBLISHED_DUPLICATE rather than letting "
                "row order pick the pairing"
            )
        named = [g for g in group if g["source"] == preferred]
        if len(named) != 1:
            raise SystemExit(
                f"tips={tips}: the declared published source {preferred!r} matched {len(named)} "
                f"of the timed trees {[g['source'] for g in group]}"
            )
        print(
            f"note: tips={tips} timed on {[g['source'] for g in group]}; "
            f"archiving the declared published source {preferred}"
        )
        out.append(named[0])
    return out


def main(out_csv: str = "benchmarks/results_scaling_py.csv") -> None:
    print(describe(MACHINE))
    paths: list[Path] = []
    for d in DIRS:
        root = Path(d)
        if root.is_dir():
            paths += sorted(root.glob("*__tree.tre"))
    rows = []
    for tp in sorted(set(paths), key=lambda p: (p.read_text().count(","), p.name)):
        tree = read_newick(tp)
        tag = scenario_tag(tp)
        for k in KS:
            ts = timed_runs(lambda: slice_pieces(tree, n=k).contribution_matrix())
            med, mn = float(np.median(ts)), float(min(ts))
            rows.append(
                {
                    "tips": tree.n_tips,
                    "n_slices": k,
                    "source": tag,
                    "slice_s": med,
                    "slice_min": mn,
                    "reps": REPS,
                    "repeat_times_s": ";".join(f"{x:.9g}" for x in ts),
                    **MACHINE,
                }
            )
            print(
                f"{tp.parent.name}/{tag} tips={tree.n_tips} k={k}: "
                f"slice median {med:.6f}s (min {mn:.6f}s)"
            )
    kept = published_rows(rows)
    pd.DataFrame(kept).to_csv(out_csv, index=False)
    print("wrote", out_csv, len(kept), "rows")


if __name__ == "__main__":
    main()
