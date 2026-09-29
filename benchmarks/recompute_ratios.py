"""Recompute every documented cross-language speedup from the archived tables.

    python benchmarks/recompute_ratios.py

A sentence such as "148-302 times faster" is only checkable if the reader knows
which cells of which table it summarises.  This script reads the archived
medians, recomputes each quoted range, prints the cells behind it, and exits
nonzero if a range no longer follows from the tables - so it doubles as a guard
on the numbers the documentation quotes.

It also states what the archive cannot support: the reference driver recorded
median, q25, q75 and minimum per cell but not its repeat vector, so a ratio can
be audited to its quartiles, not to its individual repeats.  The Python driver
now writes ``repeat_times_s``; the archived tables predate that column.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# What the documentation asserts, as (label, files, column, cell filter, expected
# value). An expected value of None means "print the range, assert nothing".
QUOTED = (
    (
        "whole ladder (headline speed-up claim)",
        ("results_scaling_py.csv", "results_scaling_r.csv"),
        "slice_s",
        lambda key: True,
        "76-556",
    ),
    (
        "ladder at 100 slices (fixed-k scaling claim)",
        ("results_scaling_py.csv", "results_scaling_r.csv"),
        "slice_s",
        lambda key: dict(key)["n_slices"] == 100,
        "148-302",
    ),
    (
        "end-to-end rate fits on 250-assemblage inputs",
        ("results_py_large.csv", "results_r_large.csv"),
        "cpd_s",
        lambda key: True,
        "6.7-8.1",
    ),
    (
        "slicing on the fixed validation grid",
        ("results_py.csv", "results_r.csv"),
        "slice_s",
        lambda key: True,
        None,
    ),
)


def load(name: str) -> dict[tuple[tuple[str, int], ...], dict[str, str]]:
    """Index a table by its dimension columns, which differ between tables."""
    with open(HERE / name, newline="") as fh:
        rows = list(csv.DictReader(fh))
    out: dict[tuple[tuple[str, int], ...], dict[str, str]] = {}
    for row in rows:
        key = tuple((col, int(row[col])) for col in ("tips", "n_sites", "n_slices") if col in row)
        out[key] = row
    return out


def describe(key: tuple[tuple[str, int], ...]) -> str:
    return " ".join(f"{name}={value}" for name, value in key)


def speedups(py: str, r: str, column: str) -> list[tuple[float, tuple[tuple[str, int], ...]]]:
    """Median ratio per paired cell, slowest agreement first."""
    try:
        left, right = load(py), load(r)
    except FileNotFoundError as exc:
        # a table that has gone must fail the check loudly; if it were simply
        # skipped the surviving cells would look like a narrower ladder
        print(f"MISSING  {exc.filename}")
        raise SystemExit(f"archived table missing: {exc.filename}")

    pairs = sorted(set(left) & set(right))
    unpaired = set(left) ^ set(right)
    if unpaired:
        names = sorted(describe(k) for k in unpaired)
        print(f"  note: {len(unpaired)} unpaired cell(s) not scored: {names}")
    out = []
    for key in pairs:
        if column in left[key] and column in right[key]:
            out.append((float(right[key][column]) / float(left[key][column]), key))
    return sorted(out)


def render(values: list[float]) -> str:
    """The range as the documentation writes it: integers when >= 100, one decimal below."""
    fmt = (lambda v: f"{v:.0f}") if max(values) >= 100 else (lambda v: f"{v:.1f}")
    return f"{fmt(min(values))}-{fmt(max(values))}"


def main() -> int:
    print("Cross-language speedups recomputed from the archived medians")
    print("  ratio = treesliceR median / PhyloSlicer median; larger is better\n")
    failures = 0
    for label, (py, r), column, keep, expected in QUOTED:
        table = [(value, key) for value, key in speedups(py, r, column) if keep(key)]
        if not table:
            print(f"MISSING  {label}: no paired cell for {column} in {py}")
            failures += 1
            continue
        computed = render([value for value, _key in table])
        verdict = "ok" if expected is None or computed == expected else "MISMATCH"
        failures += verdict == "MISMATCH"
        print(f"{verdict:8s} {label}")
        print(f"         {len(table)} cells, quoted {expected or computed}, computed {computed}")
        for value, key in table:
            print(f"           {value:9.1f}x  {describe(key)}")

    widest = speedups("results_scaling_py.csv", "results_scaling_r.csv", "slice_s")[-1]
    row = load("results_scaling_r.csv")[widest[1]]
    spread = float(row["slice_q75"]) / float(row["slice_q25"])
    verdict = "ok" if f"{spread:.1f}" == "4.3" else "MISMATCH"
    failures += verdict == "MISMATCH"
    print(f"\n{verdict:8s} spread of the widest endpoint (quoted as an IQR factor 4.3)")
    print(f"         {spread:.2f} at {describe(widest[1])}")

    print(
        "\nReproducibility limit: the reference tables store median, q25, q75 and minimum per cell,"
        "\nnot the raw repeats, so each ratio is auditable to its quartiles but not to its"
        "\nindividual runs."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
