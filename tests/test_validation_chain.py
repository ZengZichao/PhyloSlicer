"""The validation chain must not be able to validate itself.

``validation/compare.py`` and the CLI ``validate`` subcommand used to accept any
CSV as "the R reference": a transposed file was silently flipped into agreement
(at max_rel_err ~3e-14, i.e. a perfect-looking match of the wrong data), and any
single column - including phyloslicer's own output - counted as a reference.
NA handling was "verified" by comparing the *number* of non-finite values, which
two files with completely different NA patterns also pass.

These tests pin the replacements: registered hash + shape + column pattern +
orientation, refusal to run without them, refusal to transpose, and an
element-wise NA-pattern comparison.
"""

from __future__ import annotations

import importlib.util
import io
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phyloslicer.cli.main import (
    MANIFEST_FORMAT,
    ReferenceProvenanceError,
    format_reference_manifest,
    make_reference_record,
    parse_reference_manifest,
    sha256_file,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_DIR = REPO_ROOT / "validation"
REFERENCE_DIR = VALIDATION_DIR / "reference"

_spec = importlib.util.spec_from_file_location(
    "phyloslicer_validation_compare", VALIDATION_DIR / "compare.py"
)
compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(compare)

SCENARIO_FILES = [
    "bd20__tree.tre",
    "bd20__random__mat.csv",
    "bd20__random__CpD.csv",
    "bd20__random__CpE.csv",
    "bd20__random__rPD.csv",
    "bd20__random__rPE.csv",
]


def _shipped_records() -> tuple[dict, dict]:
    header, records, _ = compare.load_reference_manifest(
        REFERENCE_DIR, REFERENCE_DIR / "manifest.txt"
    )
    return header, records


@pytest.fixture
def full_scenario_dir(tmp_path):
    """The whole ``bd20__random`` scenario, beta included, as shipped.

    Every file the shipped manifest registers for that scenario is copied and
    re-registered under its shipped record, so the beta path of ``compare.py``
    runs on the real adjacency inputs and the real reference CSVs.  ``drop``
    removes a file *and* its entry, which is how a family treesliceR cannot
    produce is simulated; ``rehash`` re-pins a file's hash after a layout change,
    which is how a lying-about-its-orientation reference is simulated.
    """

    def _build(
        *, drop: tuple[str, ...] = (), tamper: str | None = None, rehash: tuple[str, ...] = ()
    ):
        target = tmp_path / "reference"
        target.mkdir(parents=True, exist_ok=True)
        _header, shipped = _shipped_records()
        names = sorted(
            n for n in shipped if n.startswith("bd20__random__") or n == "bd20__tree.tre"
        )
        records = {}
        for name in names:
            if name in drop:
                continue
            (target / name).write_bytes((REFERENCE_DIR / name).read_bytes())
            records[name] = dict(shipped[name])
        for name in rehash:
            records[name]["sha256"] = compare.sha256_file(target / name)
        _write_manifest(target / "manifest.txt", records)
        if tamper:
            victim = target / tamper
            victim.write_bytes(victim.read_bytes() + b"# tampered\n")
        return target

    return _build


def _write_manifest(path: Path, records: dict, header: dict | None = None) -> Path:
    text = format_reference_manifest(
        {"written": "test", "generator": "tests/test_validation_chain.py", **(header or {})},
        records,
    )
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def scenario_dir(tmp_path):
    """A minimal reference directory: the bd20 random scenario, fully registered."""

    def _build(
        *, omit: str | None = None, tamper: str | None = None, declare_unverified: bool = False
    ):
        target = tmp_path / "reference"
        target.mkdir(parents=True, exist_ok=True)
        for name in SCENARIO_FILES:
            if name == omit:
                continue
            (target / name).write_bytes((REFERENCE_DIR / name).read_bytes())
        records = {}
        for name in SCENARIO_FILES:
            if name == omit:
                continue
            path = target / name
            if name.endswith(".tre"):
                records[name] = {
                    "sha256": sha256_file(path),
                    "orientation": "newick-topology",
                    "columns_pattern": "-",
                    "rows": "-",
                    "cols": "-",
                }
            elif "__Cp" in name:
                letter = re.fullmatch(r".*__Cp([DE])\.csv", name).group(1)
                records[name] = make_reference_record(
                    path,
                    orientation="site-by-index",
                    columns_pattern=rf"^(Cp{letter}|P{letter}|p{letter}O)$",
                    rate_column=f"Cp{letter}",
                )
            else:
                is_profile = not name.endswith("mat.csv")
                pattern = r"^slice_[0-9]+$" if is_profile else r"^t[0-9]+$"
                if is_profile and declare_unverified:
                    # Simulate the state that prompted this gate: a
                    # profile the manifest declines to score.  The shipped files are
                    # site-by-slice now (regenerated by
                    # validation/r/generate_reference.R under R 4.5.3), so the
                    # declaration is made here rather than inherited from the
                    # reference directory - the point being that a declared
                    # reference is skipped, never quietly compared.
                    records[name] = make_reference_record(
                        path,
                        orientation="slice-by-site",
                        columns_pattern=pattern,
                        status="unverified",
                        reason="declared unverified for this test",
                    )
                    continue
                records[name] = make_reference_record(
                    path,
                    orientation="site-by-slice" if is_profile else "site-by-species",
                    columns_pattern=pattern,
                )
        _write_manifest(target / "manifest.txt", records)
        if tamper:  # after registration, so the hash no longer matches
            victim = target / tamper
            victim.write_bytes(victim.read_bytes() + b"# tampered\n")
        return target

    return _build


class TestShippedManifest:
    def test_every_shipped_reference_file_is_registered_and_untouched(self):
        header, records, manifest = compare.load_reference_manifest(
            REFERENCE_DIR, REFERENCE_DIR / "manifest.txt"
        )
        assert header["format"] == MANIFEST_FORMAT
        on_disk = {
            p.name
            for p in REFERENCE_DIR.iterdir()
            if p.suffix in (".csv", ".tre") and p.name != "equivalence_report.csv"
        }
        assert on_disk == set(records), "the manifest and the data directory disagree"
        # check_registration re-hashes every file: it must find nothing to report
        assert compare.check_registration(REFERENCE_DIR, records) == []
        assert manifest.is_file()

    def test_manifest_registers_the_regenerated_profiles(self):
        """The per-slice references are certified, and the manifest keeps the history.

        This could not be closed without R installed: the committed profiles
        contradicted their generator, so they were declared unverified rather than
        silently transposed.  They have since been regenerated by
        ``validation/r/generate_reference.R`` under R 4.5.3 / ape 5.8.1 /
        treesliceR 1.1.0, which reproduced every other reference file byte-for-byte
        - the proof that these profiles share their tree and matrix with the rest of
        the directory.  The manifest must register them now, and must still say what
        it used to record, because "this was unverifiable and here is what settled
        it" is the difference between an audit trail and a clean edit.
        """
        _header, records, _m = compare.load_reference_manifest(
            REFERENCE_DIR, REFERENCE_DIR / "manifest.txt"
        )
        header, _r = parse_reference_manifest((REFERENCE_DIR / "manifest.txt").read_text())
        profile_entries = [name for name in records if re.search(r"__rP[DE]\.csv$", name)]
        assert len(profile_entries) == 12
        for name in profile_entries:
            entry = records[name]
            assert entry["orientation"] == "site-by-slice"
            assert entry["status"] == "registered"
            assert (int(entry["rows"]), int(entry["cols"])) == (30, 50)
            assert "generate_reference.R" in entry["note"]
        assert not any(r.get("status") == "unverified" for r in records.values())
        rate_entries = [name for name in records if re.search(r"__Cp[DE]\.csv$", name)]
        assert len(rate_entries) == 12
        for name in rate_entries:
            assert records[name]["orientation"] == "site-by-index"
            assert records[name].get("rate_column") in ("CpD", "CpE")
        blob = header["limitations"]
        assert "RESOLVED 2026-09-21" in blob
        assert "status = unverified" in blob  # the superseded state, still described
        assert "byte-for-byte" in blob  # and the evidence that replaced it
        # what R did NOT buy: equal-PD slicing fails inside the reference itself
        assert "dimnames" in blob and "criterion = 'PD'" in blob
        assert not list(REFERENCE_DIR.glob("*__PD.csv")), (
            "the criterion='PD' variants must not appear as references while "
            "treesliceR 1.1.0 cannot produce them"
        )


class TestReportDeterminism:
    def test_two_fresh_runs_write_the_same_bytes(self, tmp_path):
        """The report is a function of the references, not of the machine's mood.

        Every agreement figure the documentation quotes comes out of this file.  A
        comparator that re-runs to slightly different last digits makes those
        quotes unfalsifiable, so the report has to be byte-reproducible - which
        it is only after ``validation/compare.py`` pins its BLAS threads: on a
        default-threaded interpreter 7 of the 162 rows jitter between consecutive
        runs (statuses unchanged, but a per-slice error moving 9.161e-14 ->
        9.141e-14), because multi-threaded reductions sum in a
        scheduling-dependent order.
        """
        work = tmp_path / "reference"
        shutil.copytree(REFERENCE_DIR, work)
        (work / "equivalence_report.csv").unlink()
        reports = []
        for _ in range(2):
            proc = subprocess.run(
                [sys.executable, str(VALIDATION_DIR / "compare.py"), str(work)],
                capture_output=True,
                text=True,
                env={**os.environ, "PYTHONPATH": str(REPO_ROOT / "src")},
            )
            assert proc.returncode == 0, proc.stderr[-2000:]
            reports.append((work / "equivalence_report.csv").read_bytes())
            (work / "equivalence_report.csv").unlink()
        assert reports[0] == reports[1], "the comparator is not reproducible"
        # Compared against the *shipped* report, the statuses must agree row for
        # row - and the figures only to their own tolerance.  Byte equality across
        # interpreters is not a claim this project can make: measured on
        # numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3 against the 2.4.6 run that
        # wrote the archived report, every one of the 162 statuses is identical
        # and 7 rows differ in their last digits (9.161e-14 -> 9.141e-14 style),
        # which is a library difference, not a data difference.
        shipped = pd.read_csv(REFERENCE_DIR / "equivalence_report.csv")
        fresh = pd.read_csv(io.BytesIO(reports[0]))
        assert list(fresh["scenario"]) == list(shipped["scenario"])
        assert list(fresh["status"]) == list(shipped["status"])
        scored = fresh[fresh["status"] == "verified"]
        assert (scored.max_rel_err <= scored.tol).all(), "a verified row is over tolerance"

    def test_the_thread_pin_runs_before_numpy_is_imported(self):
        """The pin only works if it precedes the numpy import - so say so."""
        source = (VALIDATION_DIR / "compare.py").read_text(encoding="utf-8")
        assert source.index("OMP_NUM_THREADS") < source.index("import numpy"), (
            "numpy has already loaded its BLAS back end by then"
        )
        for var in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            assert f'"{var}"' in source, var
        assert 'os.environ.setdefault(_var, "1")' in source


class TestRefusalToRun:
    def test_missing_manifest_is_fatal(self, tmp_path):
        target = tmp_path / "reference"
        target.mkdir()
        for name in SCENARIO_FILES:
            (target / name).write_bytes((REFERENCE_DIR / name).read_bytes())
        assert compare.main(str(target)) == 2
        assert not (target / "equivalence_report.csv").exists()

    def test_unregistered_reference_file_is_fatal(self, scenario_dir, capsys):
        target = scenario_dir(omit="bd20__random__CpD.csv")
        # a stray manifest is still needed to get past the first gate, so drop the
        # entry instead of the whole file: remove one block from the manifest
        manifest = target / "manifest.txt"
        text = manifest.read_text().replace("[file bd20__random__CpE.csv]", "[file renamed.csv]")
        manifest.write_text(text)
        assert compare.main(str(target)) == 2
        out = capsys.readouterr().err
        assert "bd20__random__CpE.csv: no entry" in out

    def test_tampered_reference_hash_is_fatal(self, scenario_dir, capsys):
        target = scenario_dir(tamper="bd20__random__mat.csv")
        assert compare.main(str(target)) == 2
        assert "bd20__random__mat.csv: bytes" in capsys.readouterr().err

    def test_cli_validate_refuses_reference_without_manifest_entry(
        self, tmp_path, random_tree_20, register_reference
    ):
        from phyloslicer.cli.main import main

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(_newick(random_tree_20))
        mat = pd.DataFrame(
            np.eye(4, random_tree_20.n_tips),
            columns=random_tree_20.tip_labels,
            index=[f"s{i}" for i in range(4)],
        )
        sites_path = tmp_path / "sites.csv"
        mat.to_csv(sites_path)
        ref_path = tmp_path / "ref.csv"
        # written with an index column, like every real treesliceR export (site_id)
        pd.DataFrame({"CpD_ref": [1.0, 2.0, 3.0, 4.0]}).to_csv(ref_path)
        register_reference(ref_path, rate_column="CpD_ref")
        # registered: accepted.  --no-strict separates the two verdicts this
        # command now reports: provenance (refused -> exit 2) versus numeric
        # agreement (a mismatch -> exit 1).  The values above are fabricated, so
        # this asserts the file was taken as a reference at all, not that it agrees.
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "4",
                    "--no-strict",
                    "--out",
                    str(tmp_path / "ok.csv"),
                ]
            )
            == 0
        )
        # same file, now governed by an empty manifest: refused
        (tmp_path / "empty.txt").write_text(
            format_reference_manifest({"note": "no records here"}, {})
        )
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "4",
                    "--out",
                    str(tmp_path / "bad.csv"),
                    "--manifest",
                    str(tmp_path / "empty.txt"),
                ]
            )
            == 2
        )

    def test_cli_validate_needs_a_declared_rate_column_and_orientation(
        self, tmp_path, random_tree_20
    ):
        from phyloslicer.cli.main import main

        tree_path = tmp_path / "in.tre"
        tree_path.write_text(_newick(random_tree_20))
        mat = pd.DataFrame(
            np.eye(4, random_tree_20.n_tips),
            columns=random_tree_20.tip_labels,
            index=[f"s{i}" for i in range(4)],
        )
        sites_path = tmp_path / "sites.csv"
        mat.to_csv(sites_path)
        ref_path = tmp_path / "ref.csv"
        pd.DataFrame({"CpD_ref": [1.0, 2.0, 3.0, 4.0]}).to_csv(ref_path)
        # registered but with no rate_column: the CLI must not fall back to iloc[:,0]
        record = make_reference_record(
            ref_path, orientation="site-by-index", columns_pattern=r"^CpD_ref$"
        )
        record.pop("rate_column", None)
        _write_manifest(tmp_path / "manifest.txt", {ref_path.name: record})
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "4",
                    "--manifest",
                    str(tmp_path / "manifest.txt"),
                    "--out",
                    str(tmp_path / "bad.csv"),
                ]
            )
            == 2
        )
        # registered as a per-slice profile: right file, wrong kind of reference
        record["rate_column"] = "CpD_ref"
        record["orientation"] = "slice-by-site"
        _write_manifest(tmp_path / "manifest.txt", {ref_path.name: record})
        assert (
            main(
                [
                    "validate",
                    str(tree_path),
                    str(sites_path),
                    "--reference",
                    str(ref_path),
                    "--n-slices",
                    "4",
                    "--manifest",
                    str(tmp_path / "manifest.txt"),
                    "--out",
                    str(tmp_path / "bad2.csv"),
                ]
            )
            == 2
        )


class TestOrientationIsEnforced:
    def test_transposed_profile_reference_raises(self, tmp_path):
        """A 30x50 file registered as 50x30 is refused, naming both shapes."""
        frame = pd.DataFrame(
            np.arange(30 * 50, dtype=float).reshape(30, 50),
            columns=[f"slice_{j + 1}" for j in range(50)],
            index=[f"site_{i + 1}" for i in range(30)],
        )
        path = tmp_path / "bd20__random__rPD.csv"
        frame.to_csv(path)
        # registered as the reference *should* look like: 50 slices x 30 sites
        record = make_reference_record(
            path,
            orientation="slice-by-site",
            columns_pattern=r"^slice_[0-9]+$",
            rows=50,
            cols=30,
        )
        with pytest.raises(ReferenceProvenanceError) as exc:
            compare.read_profile_reference(path.name, tmp_path, record)
        message = str(exc.value)
        assert "bd20__random__rPD.csv" in message
        assert "30x50" in message and "50x30" in message
        assert "slice-by-site" in message
        assert "transpos" in message

    def test_wrongly_declared_column_names_are_refused(self, tmp_path):
        frame = pd.DataFrame(
            np.zeros((50, 30)),
            columns=[f"site_{j + 1}" for j in range(30)],
            index=range(50),
        )
        path = tmp_path / "bd20__random__rPD.csv"
        frame.to_csv(path)
        record = make_reference_record(
            path, orientation="slice-by-site", columns_pattern=r"^slice_[0-9]+$"
        )
        with pytest.raises(ReferenceProvenanceError, match=r"registered pattern"):
            compare.read_profile_reference(path.name, tmp_path, record)

    def test_the_auto_transpose_is_really_gone(self):
        source = (VALIDATION_DIR / "compare.py").read_text()
        assert not re.search(r"\.T\b", source), "compare.py transposes again"
        assert "legacy column-bound layout" not in source

    def test_mismatch_survives_into_the_report(self, scenario_dir):
        """The registered orientation is what the comparison is checked against.

        A file transposed relative to its own registration must fail and name both
        shapes, instead of being flipped into a 1e-14 agreement.  This used to be
        demonstrated by re-registering the shipped 50x30 rPD file as the generator's
        30x50 layout; that file has since been regenerated by
        ``validation/r/generate_reference.R`` under R 4.5.3 and is genuinely
        site-by-slice now, so the transposition is manufactured on disk here rather
        than inherited from the reference directory.
        """
        target = scenario_dir()
        manifest = target / "manifest.txt"
        header, records = parse_reference_manifest(manifest.read_text())
        entry = records["bd20__random__rPD.csv"]
        assert (int(entry["rows"]), int(entry["cols"])) == (30, 50)
        assert entry["orientation"] == "site-by-slice"
        path = target / "bd20__random__rPD.csv"
        frame = pd.read_csv(path, index_col=0)
        assert frame.shape == (30, 50), frame.shape
        frame.T.to_csv(path)  # 50 rows x 30 columns: a transposition
        entry["sha256"] = sha256_file(path)  # the hash is honest; only the layout lies
        _write_manifest(
            manifest,
            records,
            header={k: v for k, v in header.items() if k != "format"},
        )
        assert compare.main(str(target)) == 1
        report = pd.read_csv(target / "equivalence_report.csv")
        row = report[report["scenario"] == "bd20__random:PD"].iloc[0]
        assert row["status"] == "provenance failure"
        assert "50x30" in row["note"] and "30x50" in row["note"]


class TestNaPatternsAndTotals:
    def test_equal_na_counts_are_not_equal_na_patterns(self):
        """The check the review called out: counting NAs passes, matching fails."""
        mine = np.array([np.nan, 1.0, 2.0, 3.0])
        ref = np.array([4.0, np.nan, 6.0, 7.0])
        assert int(np.sum(~np.isfinite(mine))) == int(np.sum(~np.isfinite(ref)))  # old "pass"
        note, agree = compare.na_pattern_report(mine, ref, label="pDO")
        assert agree is False
        assert "pDO: NA pattern differs at 2 row(s)" in note
        # identical patterns (both NA at the same row) do agree
        note2, agree2 = compare.na_pattern_report(mine, np.array([np.nan, 1, 2, 3]), label="pDO")
        assert agree2 is True and "1 NA rows" in note2

    def test_report_carries_reference_totals_and_na_patterns(self, scenario_dir):
        """The report must separate the three kinds of reference it compared.

        A3 is closed: the per-slice ``*__rPD/rPE`` profiles were regenerated by
        ``validation/r/generate_reference.R`` under R 4.5.3 - which reproduced every
        other reference file byte-for-byte - and are hash-registered, so profile
        agreement is now a scored row rather than "not evaluated".  The whole-tree
        ``PD`` / ``PE`` totals that treesliceR prints beside its rates remain their
        own ``:PD_total`` / ``:PE_total`` rows, because they are a different
        quantity with a different reference column.  The row count is whatever those
        kinds add up to, so this asserts the kinds and their ``status`` /
        ``ref_basis`` instead of a number that changes whenever a kind is added.
        """
        target = scenario_dir()
        assert compare.main(str(target)) == 0
        report = pd.read_csv(target / "equivalence_report.csv")
        assert set(report["status"]) == {"verified"}, report[report["status"] != "verified"][
            ["scenario", "status"]
        ].to_string()
        # "<scenario>:<kind>" - every row belongs to exactly one of the three kinds
        kind = report["scenario"].str.rsplit(":", n=1, expand=True)[1]
        rates = report[kind.isin(["CpD", "CpE"])]
        profiles = report[kind.isin(["PD", "PE"])]
        totals = report[kind.isin(["PD_total", "PE_total"])]
        assert len(rates) + len(profiles) + len(totals) == len(report)
        assert len(rates) == len(profiles) == len(totals) == 2  # D and E families

        for _, row in rates.iterrows():
            letter = row["scenario"].rsplit(":", 1)[1][-1]  # D / E
            # the reference stores the rate as CpD/CpE, the whole-tree value as
            # PD/PE and the origin time as pDO/pEO
            columns = (f"Cp{letter}", f"P{letter}", f"p{letter}O")
            assert row["status"] == "verified"
            assert row["ref_basis"] == f"bd20__random__{columns[0]}.csv:{columns[0]}"
            assert "NA patterns identical" in row["note"]
            # every compared column has its NA pattern checked element by
            # element - counting NAs was the check this replaced
            for column in columns:
                assert f"{column}: NA patterns identical" in row["note"]
            # the PD / PE total columns of the reference are really compared now
            assert f"P{letter} max rel err" in row["note"]
            assert float(row["max_rel_err"]) < 1e-4

        # the certified basis: the whole-tree totals, sourced from the *rate*
        # reference (never from the unverified profile files)
        for _, row in totals.iterrows():
            column = row["scenario"].rsplit(":", 1)[1].removesuffix("_total")  # PD / PE
            assert row["status"] == "verified"
            assert row["ref_basis"] == f"bd20__random__Cp{column[-1]}.csv:{column}"
            assert float(row["max_rel_err"]) < 1e-8  # deterministic tolerance
            assert float(row["tol"]) == compare.TOL_PROFILE
            assert float(row["pass_rate"]) == 1.0
            assert "not a per-slice profile" in row["note"]

        # ... and the per-slice profiles, now scored against a registered reference
        for _, row in profiles.iterrows():
            column = row["scenario"].rsplit(":", 1)[1]  # PD / PE
            assert row["status"] == "verified"
            assert row["ref_basis"] == f"bd20__random__r{column}.csv"
            assert float(row["max_rel_err"]) < compare.TOL_PROFILE
            assert float(row["pass_rate"]) == 1.0
        # a total is not a profile, even though both pass at the same tolerance
        assert "not a per-slice profile" in totals.iloc[0]["note"]

    def test_declared_unverified_is_skipped_not_passed(self, scenario_dir):
        """A profile the manifest declines to score is reported as such.

        The shipped reference no longer triggers this path (A3 is closed), but the
        behaviour is the whole reason the manifest exists: a skipped comparison must
        read "not evaluated" with no error figure, never a silent pass.
        """
        target = scenario_dir(declare_unverified=True)
        assert compare.main(str(target)) == 0
        report = pd.read_csv(target / "equivalence_report.csv")
        rows = report[report["scenario"].isin(["bd20__random:PD", "bd20__random:PE"])]
        assert len(rows) == 2 and set(rows["status"]) == {"not evaluated"}
        for _, row in rows.iterrows():
            assert "declared unverified for this test" in row["note"]
            assert "declared unverified in the manifest" in row["ref_basis"]
            # a skipped row must never masquerade as a pass
            assert pd.isna(row["pass_rate"]) and pd.isna(row["max_rel_err"])

    def test_strict_mode_fails_on_declared_unverified(self, scenario_dir):
        target = scenario_dir(declare_unverified=True)
        assert compare.main(str(target), strict=True) == 1


def _newick(tree) -> str:
    from phyloslicer.io import write_newick

    return write_newick(tree) + "\n"


# ---------------------------------------------------------------------- #
# the beta (CpB / CpB_RW / r_phylo(index="PB")) reference chain
# ---------------------------------------------------------------------- #
BETA_RATE_SCORING = {
    # family -> status compare.py must give it on the shipped references
    "CpB__sorensen_multisite": "verified",
    "CpB__turnover_multisite": "documented divergence",
    "CpB__nestedness_multisite": "documented divergence",
    "CpB__sorensen_pairwise": "no reference",
    "CpB__turnover_pairwise": "documented divergence",
    "CpB__nestedness_pairwise": "documented divergence",
    "CpB_RW__multisite": "verified",
    "CpB_RW__pairwise": "verified",
}
BETA_PROFILE_SCORING = {
    "rPB__sorensen_multisite": "verified",
    "rPB__turnover_multisite": "verified",
    "rPB__nestedness_multisite": "verified",
    "rPB__sorensen_pairwise": "no reference",
    "rPB__turnover_pairwise": "no reference",
    "rPB__nestedness_pairwise": "verified",
}


class TestBetaAdjacencySemantics:
    """The exported adjacency must be *the* neighbourhood PhyloSlicer computes.

    A beta reference is only comparable if both sides read the same focal
    neighbourhoods, so the queen relation the R generator used
    (``queen_adj(30)``, ``validation/r/generate_reference.R``) is checked against
    :func:`phyloslicer.spatial.adjacency_from_grid` here rather than asserted in a
    comment.
    """

    @staticmethod
    def _grid(n_sites: int) -> pd.DataFrame:
        side = int(np.ceil(np.sqrt(n_sites)))
        # same cell order as the generator: x = rep(1:side, times=side),
        # y = rep(1:side, each=side), truncated to n_sites
        cells = [(x, y) for y in range(1, side + 1) for x in range(1, side + 1)][:n_sites]
        return pd.DataFrame(
            {
                "minx": [c[0] - 1 for c in cells],
                "miny": [c[1] - 1 for c in cells],
                "maxx": [c[0] for c in cells],
                "maxy": [c[1] for c in cells],
            }
        )

    def test_queen_grid_reproduces_every_exported_adjacency(self):
        from phyloslicer.spatial import adjacency_from_grid

        _header, shipped = _shipped_records()
        adj_names = sorted(n for n in shipped if n.endswith("__adj.csv"))
        assert len(adj_names) == 6, adj_names
        for name in adj_names:
            ref = pd.read_csv(REFERENCE_DIR / name, index_col=0)
            assert (int(shipped[name]["rows"]), int(shipped[name]["cols"])) == ref.shape
            mine = adjacency_from_grid(self._grid(ref.shape[0]), method="queen")
            assert np.array_equal(mine, ref.to_numpy()), f"{name}: queen adjacency differs"
            # rook must NOT equal it, or the check above is vacuous
            rook = adjacency_from_grid(self._grid(ref.shape[0]), method="rook")
            assert not np.array_equal(rook, ref.to_numpy()), f"{name}: corner contacts missing"

    def test_exported_adjacency_is_reflexive_symmetric_and_nontrivial(self):
        _header, shipped = _shipped_records()
        for name in sorted(n for n in shipped if n.endswith("__adj.csv")):
            adj = pd.read_csv(REFERENCE_DIR / name, index_col=0).to_numpy()
            assert (adj == adj.T).all(), f"{name}: not symmetric"
            assert np.diag(adj).all(), f"{name}: a focal site is not its own neighbour"
            sizes = adj.sum(axis=1)
            assert sizes.min() >= 2, f"{name}: a neighbourhood has fewer than two sites"
            assert sizes.max() > sizes.min(), f"{name}: every neighbourhood the same size"


class TestShippedBetaReferences:
    def test_manifest_registers_the_beta_families_and_inputs(self):
        header, shipped = _shipped_records()
        beta = {n: r for n, r in shipped.items() if "__CpB" in n or "__rPB" in n}
        assert len(beta) == 66, sorted(beta)  # 7 rate + 4 profile families x 6 scenarios
        assert header.get("beta_families"), "the manifest must name the beta set"
        for name, rec in beta.items():
            for field in ("sha256", "rows", "cols", "orientation", "columns_pattern", "status"):
                assert rec.get(field), f"{name}: manifest entry lacks {field}"
            assert rec["status"] == "registered", name
            assert "generate_reference.R" in rec["note"], name
        for name, rec in shipped.items():
            if not name.endswith("__adj.csv"):
                continue
            assert rec["orientation"] == "focal-site-by-site", name
            assert (rec["rows"], rec["cols"]) == ("30", "30"), name
            assert rec.get("rate_column") is None, f"{name}: an input has no rate column"
        for name, rec in beta.items():
            if "__CpB_RW__" in name:
                assert rec["rate_column"] == "CpB_RW", name
                assert rec["total_column"] == "PB_RW", name
            elif "__CpB__" in name:
                assert rec["rate_column"] == "CpB", name
                assert rec["total_column"] == "PB", name
            else:
                assert rec["orientation"] == "site-by-slice", name
                assert (rec["rows"], rec["cols"]) == ("30", "50"), name

    def test_unsupported_families_have_no_reference_and_no_entry(self):
        """Families treesliceR cannot produce stay out of the registered set.

        A silently-forced file - a ragged k x P block squeezed into 30 x 50, or a
        per-pair vector collapsed into one number - would be worse than nothing,
        because it would look like a comparison.
        """
        _header, shipped = _shipped_records()
        for name in (
            "bd20__random__CpB__sorensen_pairwise.csv",
            "bd20__random__rPB__sorensen_pairwise.csv",
            "bd20__random__rPB__turnover_pairwise.csv",
        ):
            assert not (REFERENCE_DIR / name).exists(), name
            assert name not in shipped, name
        assert not list(REFERENCE_DIR.glob("*__CpB__sorensen_pairwise.csv"))
        assert compare.NO_REFERENCE.keys() == {
            "CpB__sorensen_pairwise",
            "rPB__sorensen_pairwise",
            "rPB__turnover_pairwise",
        }
        for suffix, why in compare.NO_REFERENCE.items():
            assert "treesliceR" in why, suffix


class TestBetaScoring:
    def test_no_adjacency_means_no_beta_rows(self, scenario_dir):
        """A scenario without the neighbourhood input contributes no beta row.

        This is the guard against a vacuous pass: without ``__adj.csv`` there is
        nothing to compare, and the report must say so by having no rows at all
        rather than by scoring an empty agreement.
        """
        target = scenario_dir()
        assert compare.main(str(target)) == 0
        report = pd.read_csv(target / "equivalence_report.csv")
        assert not report["scenario"].str.contains("CpB|rPB").any()

    def test_beta_rows_are_scored_as_expected(self, full_scenario_dir):
        target = full_scenario_dir()
        assert compare.main(str(target)) == 0, "a scored divergence must not fail the run"
        report = pd.read_csv(target / "equivalence_report.csv")
        beta = report[
            report["scenario"].str.split(":", n=1).str[1].str.startswith(tuple(compare.BETA_LABELS))
        ]
        # derive the family column from the subset, not from the full report: a
        # boolean Series taken off `report` and applied to `beta` is reindexed by
        # label, which silently keeps rows it should have dropped
        fam = beta["scenario"].str.split(":", n=1).str[1]
        # 8 rate families + 7 whole-tree totals (the family with no reference has
        # no total either) + 6 profile families
        assert len(beta) == 21, len(beta)

        for suffix, expected in BETA_RATE_SCORING.items():
            row = beta[fam == suffix]
            assert len(row) == 1, suffix
            assert row["status"].iloc[0] == expected, (
                suffix,
                row[["status", "max_rel_err", "note"]].to_string(),
            )
            if expected == "no reference":
                assert pd.isna(row["pass_rate"].iloc[0])
                assert row["ref_basis"].iloc[0] == "none"
                assert "treesliceR" in row["note"].iloc[0]
                continue
            assert row["n"].iloc[0] == 30
            assert f"bd20__random__{suffix}.csv:" in row["ref_basis"].iloc[0]
            # every scored rate row carries the deterministic whole-tree index and
            # the control fit, so a divergence can be read off the report itself
            assert "NA patterns identical" in row["note"].iloc[0]
            assert "control fit" in row["note"].iloc[0]
            assert "max rel err" in row["note"].iloc[0]

        # the deterministic half of the beta claim holds for every family,
        # including the four whose *rate* is an accepted divergence
        for suffix in BETA_RATE_SCORING:
            stem = "CpB_RW" if suffix.startswith("CpB_RW") else "CpB"
            total = "PB_RW" if stem == "CpB_RW" else "PB"
            if BETA_RATE_SCORING[suffix] == "no reference":
                assert not (fam == f"{suffix}:{total}").any(), suffix
                continue
            row = beta[fam == f"{suffix}:{total}"]
            assert len(row) == 1, suffix
            assert row["status"].iloc[0] == "verified", suffix
            assert float(row["tol"].iloc[0]) == compare.TOL_PROFILE
            assert float(row["max_rel_err"].iloc[0]) < compare.TOL_PROFILE

        for suffix, expected in BETA_PROFILE_SCORING.items():
            row = beta[fam == suffix]
            assert len(row) == 1, suffix
            assert row["status"].iloc[0] == expected, suffix
            if expected == "verified":
                assert row["n"].iloc[0] == 1500
                assert float(row["max_rel_err"].iloc[0]) < compare.TOL_PROFILE
                assert row["ref_basis"].iloc[0] == f"bd20__random__{suffix}.csv"

    def test_documented_divergence_is_not_counted_as_verified(self, full_scenario_dir):
        target = full_scenario_dir()
        assert compare.main(str(target)) == 0
        report = pd.read_csv(target / "equivalence_report.csv")
        div = report[report["status"] == "documented divergence"]
        assert len(div) == 4, div[["scenario", "max_rel_err"]].to_string()
        # admitted only because the error is large, not because it is small
        assert (div["max_rel_err"] > 1e-2).all()
        assert (div["pass_rate"] < 0.95).all()
        # ... and only when the reference-convention control closes it
        assert (div["note"].str.contains("control fit")).all()
        assert (div["note"].str.contains("accepted design divergence")).all()
        # no divergence row may also be a verified row
        assert set(div["scenario"]).isdisjoint(
            set(report[report["status"] == "verified"]["scenario"])
        )

    def test_strict_mode_fails_on_beta_gaps_and_divergences(self, full_scenario_dir):
        target = full_scenario_dir()
        assert compare.main(str(target)) == 0
        assert compare.main(str(target), strict=True) == 1

    def test_tampered_beta_reference_is_fatal(self, full_scenario_dir, capsys):
        target = full_scenario_dir(tamper="bd20__random__CpB__sorensen_multisite.csv")
        assert compare.main(str(target)) == 2
        assert "CpB__sorensen_multisite.csv: bytes" in capsys.readouterr().err

    def test_unregistered_adjacency_is_fatal(self, full_scenario_dir, capsys):
        target = full_scenario_dir()
        (target / "manifest.txt").write_text(
            (target / "manifest.txt")
            .read_text()
            .replace("[file bd20__random__adj.csv]", "[file renamed_adj.csv]")
        )
        assert compare.main(str(target)) == 2
        assert "bd20__random__adj.csv: no entry" in capsys.readouterr().err

    def test_rpb_reference_with_a_lying_orientation_is_refused(self, full_scenario_dir):
        """A transposed per-slice beta profile cannot be scored as if it fit.

        The hash is re-pinned over the transposed bytes, so the only thing left
        for the comparator to complain about is the *layout*.  Without that
        re-pin the run exits 2 on a provenance error before it ever looks at the
        shape, which would leave hash enforcement and shape enforcement
        indistinguishable.
        """
        target = full_scenario_dir()
        path = target / "bd20__random__rPB__sorensen_multisite.csv"
        pd.read_csv(path, index_col=0).T.to_csv(path)
        _header, shipped = _shipped_records()
        records = {
            n: dict(r)
            for n, r in shipped.items()
            if n.startswith("bd20__random__") or n == "bd20__tree.tre"
        }
        # keep the shipped 30x50 site-by-slice registration and re-pin only the
        # bytes: the file now claims to be that layout while being its transpose
        records[path.name]["sha256"] = compare.sha256_file(path)
        _write_manifest(target / "manifest.txt", records)
        assert compare.main(str(target)) == 1
        report = pd.read_csv(target / "equivalence_report.csv")
        row = report[report["scenario"] == "bd20__random:rPB__sorensen_multisite"].iloc[0]
        assert row["status"] == "provenance failure"
        assert "bytes" not in row["note"], row["note"]
        assert "50x30" in row["note"] and "30x50" in row["note"], row["note"]

    def test_missing_reference_degrades_to_a_named_gap(self, full_scenario_dir):
        """Drop a family's file: the row reports a gap, never a silent pass."""
        target = full_scenario_dir(drop=("bd20__random__CpB_RW__pairwise.csv",))
        assert compare.main(str(target)) == 0
        report = pd.read_csv(target / "equivalence_report.csv")
        row = report[report["scenario"] == "bd20__random:CpB_RW__pairwise"].iloc[0]
        assert row["status"] == "no reference"
        assert row["ref_basis"] == "none"
        assert pd.isna(row["max_rel_err"]) and pd.isna(row["pass_rate"])
        assert "no file for this family" in row["note"]


class TestReferenceConventionControl:
    """The comparator's own control fit, on data with a known answer."""

    def test_recovers_the_rate_of_a_scaled_curve(self):
        ages = np.linspace(1.0, 0.02, 50)
        truth = np.array([2.0, 7.5])
        # Two rows whose curves are scaled differently.  A row's total is *not*
        # the sum of its profile, so the reference convention (cumsum/whole-tree
        # index, CpB.R:392) and the shipped one (renormalise the profile to sum
        # to 1) are different objectives here; only the former returns
        # the rate the curve was built from.
        series, total = [], []
        for rate, amplitude in zip(truth, (0.9, 1.3)):
            curve = amplitude * np.exp(-rate * ages)
            profile = np.empty_like(curve)
            profile[0] = curve[0]
            profile[1:] = np.diff(curve)  # cumsum(profile) = curve
            series.append(profile)
            total.append(amplitude)  # cumsum/total = exp(-rate*ages)
        series = np.vstack(series)
        got = compare.fit_reference_convention(series, np.asarray(total), ages)
        assert np.allclose(got, truth, rtol=1e-4), got
        # The control must not be a no-op: scoring the same curves under the
        # profile-sum normalisation biases the rate by 4.7% and 13.4%, which is
        # the documented divergence this function exists to attribute.
        naive = compare.fit_reference_convention(series, series.sum(axis=1), ages)
        assert np.all(np.abs(naive / truth - 1) > 0.02), naive

    def test_degenerate_rows_come_back_nan_not_zero(self):
        ages = np.linspace(1.0, 0.02, 20)
        series = np.vstack([np.exp(-3.0 * ages), np.full(20, np.nan), np.zeros(20)])
        got = compare.fit_reference_convention(series, np.array([1.0, 1.0, 0.0]), ages)
        assert np.isfinite(got[0])
        assert np.isnan(got[1]) and np.isnan(got[2])

    def test_pairwise_turnover_uses_the_reference_curve_control(self):
        """The substitution is recorded in the row, so it cannot be mistaken."""
        assert compare.CURVE_CONTROL == {"CpB__turnover_pairwise": "multisite"}
        # CpB.R:605 (pairwise) and CpB.R:609 (multisite) are mean/mean and
        # sum/sum of the same quantity, so treesliceR's own turnover rate must be
        # method-independent - measured on the shipped references:
        for tag in ("bd20__random", "bd20__clustered", "coal100__clustered"):
            ms = pd.read_csv(REFERENCE_DIR / f"{tag}__CpB__turnover_multisite.csv", index_col=0)
            pw = pd.read_csv(REFERENCE_DIR / f"{tag}__CpB__turnover_pairwise.csv", index_col=0)
            both = ms["CpB"].notna() & pw["CpB"].notna()
            rel = (ms["CpB"][both] - pw["CpB"][both]).abs() / ms["CpB"][both].abs()
            assert float(rel.max()) < 1e-6, (tag, float(rel.max()))
            assert not np.allclose(ms["PB"], pw["PB"]), f"{tag}: PB should still differ"
