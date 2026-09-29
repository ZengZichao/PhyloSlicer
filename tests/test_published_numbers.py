"""The documented speed-up ranges must follow from the archived tables.

These are the only tests that read the quoted numbers rather than the package's.
They exist because the documented claims are ranges over a benchmark ladder, and
a range is easy to leave behind when a table is regenerated: the guard is the
same recomputation CI runs, exercised as a test.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks" / "recompute_ratios.py"

QUOTED = ("quoted 76-556", "quoted 148-302", "quoted 6.7-8.1")


def _run(script: Path):
    return subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True, check=False
    )


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
