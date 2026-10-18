"""Recompute every documented cross-language speedup from the archived tables.

    python benchmarks/recompute_ratios.py

A sentence such as "148-302 times faster" is only checkable if the reader knows
which cells of *which table* it summarises, and which input tree each of those
cells was measured on.  This script reads the archived medians, recomputes each
quoted range, prints the table and the scenario behind every cell, and exits
nonzero when

* a table has gone missing (a shrunken ladder must not look like a narrower
  claim), or when its row set silently differs from the other side's;
* a row that feeds a published number carries no scenario tag, or the tag
  ``unknown`` - the pairing would then rest on a guess;
* a tree size is covered by two input trees and the pairing is not the declared
  one (``benchmarks/README.md`` caveat 2 and the SI quote ``bd100``, so that is
  what is selected here, by name, never by row order);
* a quoted range no longer follows from the tables;
* a table no longer records the repeat count its documented protocol claims.

What the archive does *not* support, and what this script therefore prints
rather than asserts: the archived R tables store median, q25, q75 and minimum
per cell but no repeat vector, so a ratio is auditable to its quartiles, not to
its individual runs.  ``results_py.csv`` / ``results_r.csv`` hold one row per
(tree, matrix-structure, k) triple - twelve rows for four (tips, n_sites,
n_slices) cells - and every one of those rows is itself already the median of
five timed runs.  They are *not* five repeats of one cell, so collapsing them
into a single median would be a statistic nobody measured; the pairing is
therefore one-to-one on the full key, and nothing is chosen by row order.

Session drift affects both engines, not only the R one, so the script also
prints, for the tree sizes measured by two protocols, the two archived Python
and R values side by side next to the ratio each table yields.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

#: columns that identify a measured cell; a table may carry any subset
KEY_COLUMNS = ("tips", "n_sites", "n_slices", "source", "structure")
#: tags meaning "nobody can say which input produced this row"
UNKNOWN_TAGS = {"", "unknown", "-", "none"}
#: a tree size covered by more than one input tree, and the scenario the
#: documentation quotes it against (benchmarks/README.md caveat 2, the SI figure)
PREFERRED_SOURCE = {100: "bd100"}
#: the pairing the documentation does *not* quote.  It is computed and printed as
#: an appendix, so the choice and its consequences are visible instead of implied.
ALTERNATE_SOURCE = {100: "coal100"}

#: the repeat count each protocol documents (benchmarks/README.md "Protocols")
DECLARED_REPS = {
    "results_scaling_py.csv": 15,
    "results_scaling_r.csv": 15,
    "results_py.csv": 5,
    "results_r.csv": 5,
    "results_py_large.csv": 5,
    "results_r_large.csv": 5,
}

# What the documentation asserts, as (label, (python table, R table), column,
# cell filter, expected value).  An expected value of None means "print the
# range, assert nothing"; the two `_large` tables are the end-to-end protocol,
# the two `_scaling` tables the fixed-k kernel ladder, `results_{py,r}.csv` the
# fixed validation grid.
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

#: the cells two protocols both measure, reported as the drift between them
DRIFT_CELLS = ((1000, 100), (2000, 100))

_PAIRINGS: dict[tuple[str, str, str], tuple[list, str, tuple[str, ...]]] = {}
_TABLES: dict[str, list[dict[str, str]]] = {}


def _read(name: str) -> list[dict[str, str]]:
    if name in _TABLES:
        return _TABLES[name]
    path = HERE / name
    try:
        with open(path, newline="") as fh:
            rows = list(csv.DictReader(fh))
    except FileNotFoundError:
        # a table that has gone must fail the check loudly; if it were simply
        # skipped the surviving cells would look like a narrower ladder
        print(f"MISSING  {path}")
        raise SystemExit(f"archived table missing: {path}")
    if not rows:
        raise SystemExit(f"archived table {name} is empty")
    _TABLES[name] = rows
    return rows


def key_of(row: dict[str, str]) -> tuple[tuple[str, int | str], ...]:
    """The identity of a row: the dimension columns it happens to carry.

    The tables differ in their dimensions, and ``tips`` alone is not a key - at
    100 tips the scaling table holds a birth-death and a coalescent tree, and
    the validation grid holds two matrix structures per tree.
    """
    out: list[tuple[str, int | str]] = []
    for col in KEY_COLUMNS:
        if col not in row or row[col] is None:
            continue
        out.append((col, int(row[col]) if row[col].isdigit() else row[col]))
    return tuple(out)


def describe(key: tuple[tuple[str, int | str], ...]) -> str:
    return " ".join(f"{name}={value}" for name, value in key)


def _check_provenance(name: str, rows: list[dict[str, str]]) -> None:
    """Refuse to pair a row whose input tree cannot be named."""
    if "source" not in rows[0]:
        raise SystemExit(
            f"{name}: no `source` column, so the pairing for this table would rest on row "
            "order; re-run the timing driver, which writes the input tag, or attribute the "
            "rows by hand from the generator's enumeration (benchmarks/README.md records how)"
        )
    blank = [i for i, row in enumerate(rows, start=2) if not (row.get("source") or "").strip()]
    if blank:
        raise SystemExit(
            f"{name}: rows {blank[:5]} have an empty `source`, which is an untagged row: it "
            "would pair by accident rather than by provenance"
        )
    unknown = sorted(
        {row["source"] for row in rows if row["source"].strip().lower() in UNKNOWN_TAGS}
    )
    if unknown:
        raise SystemExit(
            f"{name}: scenario tag(s) {unknown} cannot be resolved to an input tree; every row "
            "of a paired table feeds a published ratio, so an untagged row must be re-measured "
            "or re-attributed, never guessed at"
        )
    seen: dict[tuple[tuple[str, int | str], ...], int] = {}
    for i, row in enumerate(rows, start=2):
        key = key_of(row)
        if key in seen:
            raise SystemExit(
                f"{name}: {describe(key)} appears twice (rows {seen[key]} and {i}); a "
                "dict-indexed pairing would silently keep the last one"
            )
        seen[key] = i


def _check_reps(name: str, rows: list[dict[str, str]]) -> str:
    """Confirm the archived repeat count is the one the protocol documents."""
    expected = DECLARED_REPS.get(name)
    if expected is None:
        return "no documented repeat count for this table"
    absent = [i for i, row in enumerate(rows, start=2) if not (row.get("reps") or "").strip()]
    if absent:
        raise SystemExit(
            f"{name}: rows {absent[:5]} carry no `reps`, so the claim that this table holds a "
            f"median of {expected} repeats is not checkable from the archive"
        )

    got = sorted({row["reps"] for row in rows})
    if got != [str(expected)]:
        raise SystemExit(
            f"{name}: `reps` reads {got} but benchmarks/README.md documents {expected} "
            "repeats for this protocol"
        )
    return f"median of {expected} timed repeats"


def _resolve_sizes(
    py_name: str, r_name: str, py: dict, r: dict, preferred: dict[int, str]
) -> tuple[dict[int, str | None], list[str]]:
    """Which scenario tags may be paired, per tree size.

    Returns ``{tips: declared source or None}`` plus the lines explaining the
    decision.  ``None`` means both tables carry the same scenario set for that
    size, so every scenario pairs one-to-one with itself.  A set that differs -
    the 100-tip rung of the scaling tables, where only the R side archived
    ``coal100`` - has to be resolved by the name the documentation quotes, and
    an undeclared disagreement is an error rather than a coin toss.
    """

    def sources(side: dict) -> dict[int, set[str]]:
        out: dict[int, set[str]] = {}
        for key in side:
            item = dict(key)
            out.setdefault(item["tips"], set()).add(item["source"])
        return out

    py_sizes, r_sizes = sources(py), sources(r)
    if set(py_sizes) != set(r_sizes):
        only_py = sorted(set(py_sizes) - set(r_sizes))
        only_r = sorted(set(r_sizes) - set(py_sizes))
        raise SystemExit(
            f"{py_name} and {r_name} cover different tree sizes: {only_py} only on the Python "
            f"side, {only_r} only on the R side.  Pairing a partial ladder would print a "
            "narrower claim than the one the documentation makes."
        )
    policy: dict[int, str | None] = {}
    notes: list[str] = []
    for tips in sorted(py_sizes):
        left, right = py_sizes[tips], r_sizes[tips]
        if left == right:
            policy[tips] = None
            if len(left) > 1 and tips not in preferred:
                # both sides agree, so the pairing is still one-to-one; nothing
                # is chosen by row order, but say which scenarios are in play
                notes.append(
                    f"tips={tips} is timed on {sorted(left)} on both sides; "
                    "each scenario is paired with itself"
                )
            continue
        declared = preferred.get(tips)
        if declared is None:
            raise SystemExit(
                f"tips={tips}: {py_name} has {sorted(left)} but {r_name} has {sorted(right)}; "
                "the scenario sets disagree and no published source is declared for this "
                "tree size, so there is no honest pairing"
            )
        if declared not in left or declared not in right:
            raise SystemExit(
                f"tips={tips}: the published pairing quotes {declared!r} "
                f"(PREFERRED_SOURCE) but the tables carry {sorted(left)} (Python) and "
                f"{sorted(right)} (R)"
            )
        policy[tips] = declared
        notes.append(
            f"tips={tips}: {py_name} has {sorted(left)}, {r_name} also has {sorted(right)}; "
            f"pairing only the declared published scenario {declared!r} and leaving the rest "
            "out of the ratio"
        )
    return policy, notes


def pairing(
    py_name: str, r_name: str, column: str, preferred: dict[int, str] = PREFERRED_SOURCE
) -> tuple[list, str, tuple[str, ...]]:
    """Median ratio per paired cell (smallest first), paired by (source, dimensions).

    Returns the ratios, the statistic both sides summarise, and the provenance
    notes.  Memoised: one table pair backs several claims, and a note belongs to
    the claim being printed, not to every lookup that reaches it.
    """
    memo = (py_name, r_name, column, tuple(sorted(preferred.items())))
    if memo in _PAIRINGS:
        return _PAIRINGS[memo]
    py_rows, r_rows = _read(py_name), _read(r_name)
    _check_provenance(py_name, py_rows)
    _check_provenance(r_name, r_rows)
    py_stat, r_stat = _check_reps(py_name, py_rows), _check_reps(r_name, r_rows)
    if py_stat != r_stat:
        raise SystemExit(
            f"{py_name} is {py_stat} but {r_name} is {r_stat}: the two sides would not be "
            "summarising the same number of runs"
        )
    py, r = {}, {}
    for row in py_rows:
        py[key_of(row)] = row
    for row in r_rows:
        r[key_of(row)] = row
    policy, notes = _resolve_sizes(py_name, r_name, py, r, preferred)
    out: list[tuple[float, tuple]] = []
    dropped: list[str] = []
    orphans: list[str] = []
    for key in sorted(set(py) ^ set(r)):
        # a cell only one table measured: expected where a declared published
        # scenario excludes the other one, a pairing hole everywhere else
        item = dict(key)
        chosen = policy.get(item["tips"])
        side = py_name if key in py else r_name
        if chosen is None:
            orphans.append(f"{describe(key)} (only in {side})")
        else:
            dropped.append(
                f"{describe(key)}: only in {side}, and a published ratio quotes {chosen!r} "
                "at this rung"
            )
    if orphans:
        raise SystemExit(
            f"{py_name} and {r_name} disagree about which cells exist: {', '.join(orphans)}. "
            "Pairing a hole would silently shorten the ladder."
        )
    for key in sorted(set(py) & set(r)):
        item = dict(key)
        chosen = policy.get(item["tips"])
        if chosen is not None and item["source"] != chosen:
            dropped.append(f"{describe(key)}: a published ratio quotes {chosen!r} at this rung")
            continue
        left, right = py[key], r[key]
        if column in left and column in right and right[column] not in ("", "NA"):
            out.append((float(right[column]) / float(left[column]), key))
    if not out:
        raise SystemExit(f"{py_name} and {r_name} share no pairable cell for {column}")
    notes.extend(sorted(dropped))

    result = (sorted(out), py_stat, tuple(notes))
    _PAIRINGS[memo] = result
    return result


def _cells(
    name: str, preferred: dict[int, str] = PREFERRED_SOURCE
) -> dict[tuple[int, int], dict[str, str]]:
    """Rows keyed by (tips, n_slices), keeping only the published scenario per rung."""
    return {
        (int(row["tips"]), int(row["n_slices"])): row
        for row in _read(name)
        if row["source"] == preferred.get(int(row["tips"]), row["source"])
    }


def render(values: list[float]) -> str:
    """The range as the documentation writes it: integers when >= 100, one decimal below."""
    fmt = (lambda v: f"{v:.0f}") if max(values) >= 100 else (lambda v: f"{v:.1f}")
    return f"{fmt(min(values))}-{fmt(max(values))}"


def _average_ranks(values: list) -> list[float]:
    ranks: dict = {}
    positions: dict = {}
    for i, value in enumerate(sorted(values), start=1):
        positions.setdefault(value, []).append(i)
    for value, group in positions.items():
        ranks[value] = sum(group) / len(group)
    return [ranks[v] for v in values]


def _spearman(xs: list, ys: list) -> float:
    """Rank correlation without a scipy dependency, ties averaged as scipy does."""
    return statistics.correlation(_average_ranks(xs), _average_ranks(ys))


def _summarise(title: str, cells: list[tuple[int, int, float]]) -> None:
    """Print one pairing: its ranges, its 100-tip triple, its trend across tree size."""
    values = [v for _t, _k, v in cells]
    at_100 = [v for _t, k, v in cells if k == 100]
    tips_100 = [(k, round(v, 1)) for t, k, v in cells if t == 100]
    rung_50 = sorted((t, v) for t, k, v in cells if k == 50)
    rho = _spearman([t for t, _v in rung_50], [v for _t, v in rung_50])
    print(f"         {title}")
    print(
        f"           ladder {render(values)} over {len(values)} cells, "
        f"at n_slices=100 {render(at_100)}"
    )
    print(f"           tips=100 ratios by (k, ratio): {tips_100}")
    print(
        f"           ratio vs tree size at k=50: rho {rho:+.2f} over "
        f"{[(t, round(v, 1)) for t, v in rung_50]}"
    )


def report_alternative_pairings() -> None:
    """Print the published pairing next to the one the archive used to produce.

    The archived R scaling table holds two trees at ``tips = 100`` - ``bd100``
    and ``coal100`` - while the archived Python table holds one row for that
    tree size, ``bd100`` (the driver keeps the declared published scenario; see
    ``benchmark_scaling.py``).  Indexing the old tables by ``(tips, n_slices)``
    alone therefore paired Python ``bd100`` against whichever R row came last,
    which was ``coal100``: a pairing across two different trees, quoted in a
    document that says ``bd100``.  Both readings are printed here so the choice
    is auditable rather than asserted - the three quoted ranges survive either,
    while the 100-tip triple and the trend across tree size do not, which is
    exactly why the pairing is now selected by name and the cross-tree reading is
    reported but never used.
    """
    published = [
        (dict(key)["tips"], dict(key)["n_slices"], value)
        for value, key in pairing("results_scaling_py.csv", "results_scaling_r.csv", "slice_s")[0]
    ]
    _summarise(f"published pairing {PREFERRED_SOURCE}", published)

    # the same ladder with the 100-tip rung read off the other R tree, against
    # the one Python row that exists at that size
    py = {(int(r["tips"]), int(r["n_slices"])): r for r in _read("results_scaling_py.csv")}
    r_rows = {}
    for row in _read("results_scaling_r.csv"):
        if row["source"] == ALTERNATE_SOURCE[100]:
            r_rows[(int(row["tips"]), int(row["n_slices"]))] = row
    alternate = [
        (dict(key)["tips"], dict(key)["n_slices"], value)
        for value, key in pairing("results_scaling_py.csv", "results_scaling_r.csv", "slice_s")[0]
        if dict(key)["tips"] != 100
    ]
    for cell, r_row in sorted(r_rows.items()):
        if cell in py:
            alternate.append(
                (cell[0], cell[1], float(r_row["slice_s"]) / float(py[cell]["slice_s"]))
            )
    _summarise(
        f"cross-tree reading the old row order produced "
        f"(Python bd100 vs R {ALTERNATE_SOURCE[100]}) - not a valid pairing",
        sorted(alternate),
    )


def report_drift() -> None:
    """Print the same operation as two archived protocols measured it.

    Session drift is usually attributed to treesliceR, but the Python side moves
    between sessions too - at 2000 tips and 100 slices the scaling protocol
    archived 0.95 ms while the end-to-end protocol archived 1.45 ms for the same
    call.  Both tables are real; they are simply different sessions, which is why
    every published ratio names the table it came from.
    """
    scaling, scaling_r = _cells("results_scaling_py.csv"), _cells("results_scaling_r.csv")
    large, large_r = _cells("results_py_large.csv"), _cells("results_r_large.csv")
    for tips, k in DRIFT_CELLS:
        cell = (tips, k)
        if cell not in large or cell not in scaling:
            continue
        a, b = scaling[cell], large[cell]
        ar, br = scaling_r[cell], large_r[cell]
        print(
            f"         tips={tips} n_slices={k}: Python {float(a['slice_s']) * 1e3:.4f} ms "
            f"[results_scaling_py.csv] vs {float(b['slice_s']) * 1e3:.4f} ms "
            f"[results_py_large.csv]; R {float(ar['slice_s']) * 1e3:.2f} ms vs "
            f"{float(br['slice_s']) * 1e3:.2f} ms [results_scaling_r.csv vs results_r_large.csv];"
            f" ratio {float(ar['slice_s']) / float(a['slice_s']):.1f}x vs "
            f"{float(br['slice_s']) / float(b['slice_s']):.1f}x"
        )


def main() -> int:
    print("Cross-language speedups recomputed from the archived medians")
    print("  ratio = treesliceR median / PhyloSlicer median; larger is better")
    print("  every cell below names the table pair and the input scenario it came from\n")
    failures = 0
    for label, (py, r), column, keep, expected in QUOTED:
        paired, statistic, notes = pairing(py, r, column)
        table = [(value, key) for value, key in paired if keep(key)]
        if not table:
            print(f"MISSING  {label}: no paired cell for {column} in {py}")
            failures += 1
            continue
        computed = render([value for value, _key in table])
        verdict = "ok" if expected is None or computed == expected else "MISMATCH"
        failures += verdict == "MISMATCH"
        print(f"{verdict:8s} {label}")
        print(f"         {len(table)} cells, quoted {expected or computed}, computed {computed}")
        print(f"         source:     {py} / {r}  (column {column}, {statistic} per row)")
        for line in notes:
            print(f"         provenance: {line}")
        for value, key in table:
            print(f"           {value:9.1f}x  {describe(key)}  [{py} / {r}]")

    print("\nSession drift is a property of the session, not of one engine:")
    report_drift()

    print("\nThe 100-tip rung paired both ways (the published one is bd100):")
    report_alternative_pairings()

    widest = pairing("results_scaling_py.csv", "results_scaling_r.csv", "slice_s")[0][-1]
    match = [entry for entry in _read("results_scaling_r.csv") if key_of(entry) == widest[1]]
    if not match:
        raise SystemExit(f"the widest cell {describe(widest[1])} is not in results_scaling_r.csv")
    spread = float(match[0]["slice_q75"]) / float(match[0]["slice_q25"])
    verdict = "ok" if f"{spread:.1f}" == "4.3" else "MISMATCH"
    failures += verdict == "MISMATCH"
    print(
        f"\n{verdict:8s} spread of the widest endpoint (quoted as an IQR factor 4.3)"
        f"\n         {spread:.2f} at {describe(widest[1])}  [results_scaling_r.csv]"
    )

    print(
        "\nReproducibility limit: the archived tables store median, q25, q75, minimum and the"
        "\nrepeat count per cell, but not the raw repeats, so each ratio is auditable to its"
        "\nquartiles and not to its individual runs.  Both timing drivers now write the repeat"
        "\nvector, so a re-run removes that limit; the archived tables do not contain it, and no"
        "\nratio here re-derives a median from repeats nobody archived."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
