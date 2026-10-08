"""The documented speed-up ranges must follow from the archived tables.

These are the only tests that read the quoted numbers rather than the package's.
They exist because the documented claims are ranges over a benchmark ladder, and
a range is easy to leave behind when a table is regenerated: the guard is the
same recomputation CI runs, exercised as a test.

The second half of the file guards the *provenance* of those ranges, which is the
part that used to be implicit.  ``results_scaling_r.csv`` measures the 100-tip
rung twice (``bd100`` and ``coal100``) and ``results_py.csv`` / ``results_r.csv``
hold three trees x two matrix structures per grid point, so a pairing that keys
on ``(tips, n_slices)`` alone quietly picks one of them by row order.  The
documented convention is ``bd100`` (benchmarks/README.md caveat 2, the SI
figure), so these tests assert that the pairing is made by that *name*, and that
a table which loses its tags, its repeat count, or a cell the other side still
has fails loudly instead of publishing a different number.
"""

from __future__ import annotations

import csv
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "recompute_ratios.py"
BENCH = ROOT / "benchmarks"

QUOTED = ("quoted 76-556", "quoted 148-302", "quoted 6.7-8.1")

#: every table the guard reads; a stub run needs all of them present unless a
#: test deliberately removes one
TABLES = (
    "results_scaling_py.csv",
    "results_scaling_r.csv",
    "results_py.csv",
    "results_r.csv",
    "results_py_large.csv",
    "results_r_large.csv",
)


def _run(script: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, check=False
    )


def _stub_run(
    tmp_path: Path,
    mutate: Callable[[Path], None] | None = None,
    omit: tuple[str, ...] = (),
) -> subprocess.CompletedProcess[str]:
    """Run the guard against a copy of the archived tables, optionally edited.

    Only ``HERE`` is redirected, so the pairing logic under test is the one CI
    and the documentation point at.
    """
    stage = tmp_path / "benchmarks"
    stage.mkdir(exist_ok=True)
    for name in TABLES:
        if name in omit:
            continue
        (stage / name).write_text((BENCH / name).read_text(encoding="utf-8"), encoding="utf-8")
    if mutate is not None:
        mutate(stage)
    source = SCRIPT.read_text(encoding="utf-8").replace(
        "HERE = Path(__file__).resolve().parent", f"HERE = Path({str(stage)!r})"
    )
    script = tmp_path / "recompute.py"
    script.write_text(source, encoding="utf-8")
    return _run(script)


def test_recompute_ratios_agrees_with_the_note():
    proc = _run(SCRIPT)
    assert proc.returncode == 0, f"published ranges drifted:\n{proc.stdout}\n{proc.stderr}"
    for claim in QUOTED:
        assert claim in proc.stdout, f"{claim} no longer follows from the tables:\n{proc.stdout}"


def test_missing_reference_table_is_reported_not_silently_skipped(tmp_path):
    """A table that vanishes must fail the check rather than shrink the ladder."""
    stub = tmp_path / "benchmarks"
    stub.mkdir()
    for name in ("results_scaling_py.csv", "results_scaling_r.csv"):
        (stub / name).write_text((ROOT / "benchmarks" / name).read_text(encoding="utf-8"))
    broken = SCRIPT.read_text(encoding="utf-8").replace(
        "HERE = Path(__file__).resolve().parent", f"HERE = Path({str(stub)!r})"
    )
    script = tmp_path / "recompute.py"
    script.write_text(broken, encoding="utf-8")
    proc = _run(script)
    assert proc.returncode != 0, "tables are missing here, so it must fail:\n" + proc.stdout
    assert "MISSING" in proc.stdout, proc.stdout


def test_the_100_tip_rung_is_paired_with_bd100_by_name(tmp_path):
    """The scaling ladder quotes bd100, and says so, at every 100-tip cell.

    Before the ``source`` column existed the pairing was ``(tips, n_slices)`` and
    the R side's *last* 100-tip row won - ``coal100`` - while the prose and the
    figure both said ``bd100``.  The two readings give the same quoted ranges, so
    a range-only guard could not tell them apart; this one can.
    """
    proc = _run(SCRIPT)
    assert proc.returncode == 0, proc.stdout
    ladder = proc.stdout.split("whole ladder", 1)[1].split("ladder at 100 slices", 1)[0]
    hundred = [line for line in ladder.splitlines() if re.search(r"^\s+[\d.]+x\s+tips=100 ", line)]
    assert len(hundred) == 3, hundred
    assert all("source=bd100" in line for line in hundred), hundred
    # coal100 may only show up as a named exclusion, never as a paired cell
    coal = [line for line in ladder.splitlines() if "source=coal100" in line]
    assert len(coal) == 3, coal
    assert all(line.lstrip().startswith("provenance:") for line in coal), coal
    assert "pairing only the declared published scenario 'bd100'" in proc.stdout
    assert "at this rung" in proc.stdout, proc.stdout  # the excluded coal100 rows are named


def test_the_alternate_pairing_is_reported_but_never_used(tmp_path):
    """The cross-tree reading is printed as an appendix, so the choice is visible.

    Both readings preserve the three quoted ranges, which is why the ranges alone
    were never evidence for the pairing; the 100-tip triple and the trend of the
    ratio across tree size are what differ.
    """
    proc = _run(SCRIPT)
    assert proc.returncode == 0, proc.stdout
    assert "not a valid pairing" in proc.stdout, proc.stdout
    assert "tips=100 ratios by (k, ratio): [(25, 93.8), (50, 148.9), (100, 210.8)]" in proc.stdout
    assert "tips=100 ratios by (k, ratio): [(25, 92.1), (50, 144.9), (100, 209.4)]" in proc.stdout
    assert "rho +0.43" in proc.stdout and "rho +0.54" in proc.stdout, proc.stdout


def test_every_published_ratio_names_its_table(tmp_path):
    """The mapping from claim to table is in the output, not only in the README."""
    proc = _run(SCRIPT)
    assert proc.returncode == 0, proc.stdout
    for pair in (
        "results_scaling_py.csv / results_scaling_r.csv",
        "results_py_large.csv / results_r_large.csv",
        "results_py.csv / results_r.csv",
    ):
        assert f"source:     {pair}" in proc.stdout, pair
    # the validation grid is twelve one-to-one cells, not four collapsed ones
    assert "12 cells, quoted 45-161" in proc.stdout, proc.stdout


def test_drift_is_attributed_to_both_engines(tmp_path):
    """The same Python call is archived at two values, so drift is not an R-only caveat."""
    proc = _run(SCRIPT)
    assert proc.returncode == 0, proc.stdout
    assert "Session drift is a property of the session, not of one engine" in proc.stdout
    assert "0.9516 ms [results_scaling_py.csv] vs 1.4513 ms [results_py_large.csv]" in proc.stdout
    assert "140.98 ms vs 156.17 ms" in proc.stdout
    assert "176.6x vs 104.7x" in proc.stdout and "148.2x vs 107.6x" in proc.stdout


def test_unknown_scenario_tag_fails_loudly(tmp_path):
    def mutate(stage: Path) -> None:
        path = stage / "results_scaling_py.csv"
        path.write_text(path.read_text(encoding="utf-8").replace(",bd100", ",unknown"))

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "cannot be resolved to an input tree" in proc.stderr, proc.stderr
    assert "unknown" in proc.stderr


def test_undeclared_scenario_disagreement_fails_loudly(tmp_path):
    """A rung whose two sides disagree without a declared pairing has no honest ratio."""

    def mutate(stage: Path) -> None:
        path = stage / "results_scaling_r.csv"
        path.write_text(path.read_text(encoding="utf-8").replace('"bd100"', '"bd101"'))

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "quotes 'bd100'" in proc.stderr, proc.stderr


def test_unrecorded_repeat_count_fails_instead_of_borrowing_one(tmp_path):
    """The claim "median of 5" must be checkable from the table that carries the number."""

    def mutate(stage: Path) -> None:
        path = stage / "results_py.csv"
        rows = list(csv.reader(path.read_text(encoding="utf-8").splitlines()))
        index = rows[0].index("reps")
        rows = [[cell for i, cell in enumerate(row) if i != index] for row in rows]
        path.write_text("\n".join(",".join(row) for row in rows) + "\n")

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "carry no `reps`" in proc.stderr, proc.stderr


def test_wrong_repeat_count_fails(tmp_path):
    """A table that says it summarised five runs is not the 15-run ladder."""

    def mutate(stage: Path) -> None:
        path = stage / "results_scaling_r.csv"
        rows = list(csv.reader(path.read_text(encoding="utf-8").splitlines()))
        index = rows[0].index("reps")
        for row in rows[1:]:
            row[index] = "5"
        path.write_text("\n".join(",".join(row) for row in rows) + "\n")

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "documents 15 repeats" in proc.stderr, proc.stderr


def test_a_row_the_other_table_does_not_have_is_a_hole_not_a_shrunken_ladder(tmp_path):
    """results_py.csv / results_r.csv pair one-to-one, or they fail."""

    def mutate(stage: Path) -> None:
        path = stage / "results_r.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        del lines[3]
        path.write_text("\n".join(lines) + "\n")

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "disagree about which cells exist" in proc.stderr, proc.stderr


def test_a_duplicated_key_is_not_resolved_by_row_order(tmp_path):
    def mutate(stage: Path) -> None:
        path = stage / "results_scaling_py.csv"
        lines = path.read_text(encoding="utf-8").splitlines()
        lines.append(lines[4])  # a second row claiming tips=100 n_slices=25 source=bd100
        path.write_text("\n".join(lines) + "\n")

    proc = _stub_run(tmp_path, mutate)
    assert proc.returncode != 0
    assert "appears twice" in proc.stderr, proc.stderr
