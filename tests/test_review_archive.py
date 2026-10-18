"""Tests for the anonymised peer-review archive built by ``tools/make_review_archive.py``.

The manuscript's submission checklist claims that CI "builds the review archive
through the packager, which fails the job if any personal path survives".  These
tests are what makes that claim true.  They check five things:

* the packager runs and the archive it writes passes its own identity guard;
* no author-identifying string survives in *any* archived entry, text or binary;
* no build/VCS/system artefact and no author-side tooling is archived;
* the provenance sidecar is complete - version, commit, and a SHA-256 plus a
  scrub record for every file it claims - and never quotes a removed literal;
* two builds are byte-identical, and the guard actually fails when identity is
  planted in an archive (so the passing case is not vacuous).

Finally, and this is the point of the whole exercise, the archive is unpacked
and *run*: its own test suite, the documented example and the cross-validation
comparator must all work for a reviewer who has nothing but this file.

Every identifying literal below is split into fragments, because this file is
itself archived and must be neither a hit for the guard it tests nor a victim of
the substitution it checks: no fragment is a scrub rule's literal, and no two
fragments sit adjacent in the source, so the joined string is never contiguous
in the bytes the packager reads.  The last test pins that this file ships
byte-for-byte unchanged.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PACKAGER = ROOT / "tools" / "make_review_archive.py"
PY = sys.executable

#: Set this to skip the (slow) in-archive execution tests, e.g. on a machine
#: without the optional extras.  The build/scrub/guard tests always run.
SKIP_RUN_ENV = "PHYLOSLICER_SKIP_ARCHIVE_RUN"

needs_packager = pytest.mark.skipif(
    not PACKAGER.is_file(),
    reason=(
        "tools/ is build-time tooling and is deliberately excluded from the reviewer "
        "archive (its scrub table quotes the literals it removes), so this module "
        "only runs from the source checkout"
    ),
)


#: name -> literal.  Split so that no single fragment equals a scrub rule's
#: literal and no two fragments are adjacent in these bytes: the joined strings
#: therefore never occur contiguously in this file, and neither the guard nor the
#: substitution can see them here.


def _identity_patterns() -> dict[str, str]:
    """The guard literals, read from build-time tooling that never ships.

    Raises ImportError when tools/ is absent, which is why every caller sits
    behind @needs_packager.
    """
    import importlib.util

    path = ROOT / "tools" / "identity_literals.py"
    spec = importlib.util.spec_from_file_location("phyloslicer_identity_literals", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.IDENTITY_PATTERNS


#: Entries that must never appear in a review archive.
FORBIDDEN_NAME_PARTS = (
    ".git",
    ".github",
    "tools",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    "build",
    "dist",
    ".egg-info",
)
FORBIDDEN_NAMES = (".DS_Store", "Thumbs.db")
FORBIDDEN_SUFFIXES = (".pyc", ".pyo")

#: RI-18, spelled in pieces for the same reason the identity literals are: this module is
#: archived, and a file that names the marker the guard rejects *is* a file carrying the
#: marker.  Joining two fragments leaves nothing contiguous for `marker_violations()` to
#: find, and `test_this_file_ships_unscrubbed` pins that the fragments are not scrub rules.
UNRESOLVED_MARKERS = tuple(
    "".join(parts)
    for parts in (
        ("AUTHOR ", "ACTION"),
        ("confirm", " which mode"),
        ("TODO:", " confirm"),
        ("TBD:", " confirm"),
    )
)

#: Documents that carry author identity in the source tree and must therefore
#: show up in the sidecar's scrub list.
MUST_BE_SCRUBBED = (
    "pyproject.toml",
    "CITATION.cff",
    ".zenodo.json",
    "LICENSE",
    "README.md",
    "README_CN.md",
    "SECURITY.md",
    "SECURITY_CN.md",
    "CODE_OF_CONDUCT.md",
    "CODE_OF_CONDUCT_CN.md",
    "CONTRIBUTING.md",
    "CONTRIBUTING_CN.md",
    "RELEASE.md",
    "RELEASE_CN.md",
    "tests/test_release_hygiene.py",
)

#: What a reviewer needs in order to install, run and re-verify the package.
MUST_BE_PRESENT = (
    "pyproject.toml",
    "MANIFEST.in",
    ".gitignore",
    "LICENSE",
    "README.md",
    "README_CN.md",
    "CHANGELOG.md",
    "CHANGELOG_CN.md",
    "src/phyloslicer/__init__.py",
    "tests/test_release_hygiene.py",
    "examples/quickstart.py",
    "docs/USER_GUIDE.md",
    "docs/USER_GUIDE_CN.md",
    "benchmarks/recompute_ratios.py",
    "validation/compare.py",
    "validation/reference/manifest.txt",
    "testing_suite/scripts/run_all_tests.py",
    "REVIEW_ARCHIVE_PROVENANCE.json",
    "PEER_REVIEW_NOTICE.md",
    "PEER_REVIEW_NOTICE_CN.md",
)


def build(
    out_dir: Path,
    *,
    fmt: str = "zip",
    flags: tuple[str, ...] = (),
    strict: bool = False,
) -> Path:
    """Run the packager and return the archive it wrote.

    ``strict=True`` passes no provenance flag at all, which is how the gate is
    demonstrated (the build must then *refuse*).  Otherwise the two provenance flags
    are supplied for the checkout in hand, because the packager is allowed - and now
    obliged - to name its own provenance state: ``--allow-dirty`` when the working tree
    carries uncommitted edits, ``--allow-no-git`` when there is no ``.git`` at all.
    Either way the sidecar has to record the truth about that choice, which
    :func:`test_provenance_sidecar_is_complete_and_self_consistent` then asserts.
    """
    if not strict:
        flags = (*_provenance_flags(), *flags)
    proc = subprocess.run(
        [PY, str(PACKAGER), "--out", str(out_dir), "--format", fmt, *flags],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, f"packager failed:\n{proc.stdout}\n{proc.stderr}"
    assert "identity guard PASSED" in proc.stdout, proc.stdout
    suffix = ".tar.gz" if fmt == "tar.gz" else ".zip"
    archives = sorted(out_dir.glob(f"*{suffix}"))
    assert len(archives) == 1, f"expected one archive, got {archives}"
    return archives[0]


def run_packager(*args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [PY, str(PACKAGER), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def _git(*args: str) -> str | None:
    """Read-only git query: stdout on success, ``None`` when git or the repo is absent."""
    try:
        proc = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    except OSError:
        return None
    return proc.stdout if proc.returncode == 0 else None


def checkout_state() -> dict[str, bool]:
    """Does this checkout have git metadata, and does it carry uncommitted edits?

    ``git status --porcelain`` answers with an empty *string* when the tree is clean and
    fails outright when there is no repository, so the two cases are told apart here
    rather than by ``.strip()``, which cannot.
    """
    head = _git("rev-parse", "HEAD")
    porcelain = _git("status", "--porcelain")
    return {
        "git": head is not None,
        "dirty": porcelain is not None and any(line.strip() for line in porcelain.splitlines()),
    }


def _provenance_flags() -> tuple[str, ...]:
    """The provenance flags this checkout needs before the packager will build at all."""
    state = checkout_state()
    flags = []
    if state["dirty"]:
        flags.append("--allow-dirty")
    if not state["git"]:
        flags.append("--allow-no-git")
    return tuple(flags)


def entries_of(path: Path) -> dict[str, bytes]:
    """Archive entries keyed by path *inside* the unpacked root."""
    with zipfile.ZipFile(path) as zf:
        raw = {i.filename: zf.read(i) for i in zf.infolist() if not i.is_dir()}
    prefix = {name.split("/")[0] for name in raw}
    assert len(prefix) == 1, f"archive must unpack to a single root, got {prefix}"
    return {"/".join(name.split("/")[1:]): data for name, data in raw.items()}


@pytest.fixture(scope="module")
def archive(tmp_path_factory) -> Path:
    """One shared build; the packager takes ~2 s and the run tests need bytes."""
    out = tmp_path_factory.mktemp("review-archive")
    return build(out)


@pytest.fixture(scope="module")
def entries(archive: Path) -> dict[str, bytes]:
    return entries_of(archive)


@pytest.fixture(scope="module")
def provenance(entries: dict[str, bytes]) -> dict:
    return json.loads(entries["REVIEW_ARCHIVE_PROVENANCE.json"])


@needs_packager
def test_packager_reports_and_exits_zero(tmp_path: Path) -> None:
    out = tmp_path / "out"
    proc = run_packager("--out", str(out), "--verbose", *_provenance_flags())
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "identity guard PASSED" in proc.stdout
    # the build summary must state the provenance it recorded and the language census
    for token in ("language mix       :", "provenance         :"):
        assert token in proc.stdout, f"summary is missing {token!r}"
    # summary requirements: counts, bytes, per-directory, exclusions, scrubs
    for token in (
        "entries            :",
        "uncompressed bytes",
        "per directory",
        "scrubbed files",
        "excluded entries",
    ):
        assert token in proc.stdout, f"summary is missing {token!r}"


@needs_packager
def test_no_identifying_string_survives_anywhere(entries: dict[str, bytes]) -> None:
    """Every entry, text or binary, is scanned as bytes for every pattern."""
    hits = []
    for name, data in sorted(entries.items()):
        blob = data
        for pattern, literal in _identity_patterns().items():
            needle = literal.encode("utf-8")
            if needle in blob or needle in name.encode("utf-8"):
                # a grant number inside a float column of a reference table is a
                # digit run, not a citation; the packager scopes numeric patterns
                # to documents and so do we
                if pattern in ("grant_1", "grant_2", "grant_3") and name.endswith(
                    (".csv", ".tre", ".nex", ".png")
                ):
                    continue
                hits.append((name, pattern))
    assert not hits, f"author identity survived in: {hits}"


@needs_packager
def test_artefacts_and_author_tooling_are_excluded(entries: dict[str, bytes]) -> None:
    offenders = []
    for name in entries:
        parts = name.split("/")
        if any(p in FORBIDDEN_NAME_PARTS for p in parts):
            offenders.append(name)
        elif parts[-1] in FORBIDDEN_NAMES or name.endswith(FORBIDDEN_SUFFIXES):
            offenders.append(name)
        elif name in {"REVIEW_ARCHIVE.md", "REVIEW_ARCHIVE_CN.md"}:
            offenders.append(name)
    assert not offenders, f"artefacts/author tooling shipped to reviewers: {offenders}"
    # and the things reviewers do need are there
    missing = [name for name in MUST_BE_PRESENT if name not in entries]
    assert not missing, f"reviewer archive is missing: {missing}"


@needs_packager
def test_this_file_ships_unscrubbed(entries: dict[str, bytes]) -> None:
    """The archived copy of this module must equal the source copy, byte for byte.

    If a fragment above ever became a contiguous rule literal, the packager would
    rewrite the archived copy and the assembled patterns would silently name the
    wrong strings - a test that can no longer fail.  Byte equality is the guard
    against that, and it is what makes the fragment discipline above verifiable
    rather than a comment.
    """
    source = (ROOT / "tests" / "test_review_archive.py").read_bytes()
    assert entries["tests/test_review_archive.py"] == source, (
        "this file was rewritten by the scrub: some fragment here is a rule literal"
    )


@needs_packager
def test_english_documents_keep_their_chinese_counterparts(
    entries: dict[str, bytes], provenance: dict
) -> None:
    """The project's own hygiene rule must not be broken by the selection.

    CM-19: the note promises "an English user guide, a Chinese user manual", so the
    bilingual mix is recorded in the sidecar as a census read off the shipped bytes -
    and this checks the census against the entries, the pairing, and the scrubbing: the
    Chinese documents are the mirrored half of an already-anonymised tree, not a second
    identity surface.
    """
    names = set(entries)
    missing = [
        name
        for name in names
        if name.endswith(".md")
        and not name.endswith("_CN.md")
        and "/node_modules/" not in name
        and f"{name[:-3]}_CN.md" not in names
    ]
    assert not missing, f"English documents without a Chinese counterpart: {missing}"

    mix = provenance["language_mix"]
    cjk = re.compile(r"[㐀-䶿一-鿿]")
    counted = {}
    for name, data in sorted(entries.items()):
        if Path(name).suffix.lower() not in {
            ".py",
            ".md",
            ".toml",
            ".yml",
            ".yaml",
            ".R",
            ".cff",
            ".json",
            ".txt",
            ".in",
            ".cfg",
            ".sh",
        }:
            continue
        n = len(cjk.findall(data.decode("utf-8", "ignore")))
        if n:
            counted[name] = n
    assert mix["documents_with_chinese"] == len(counted), (
        f"sidecar census {mix['documents_with_chinese']} != shipped {len(counted)}"
    )
    assert mix["total_chinese_characters"] == sum(counted.values())
    assert mix["markdown_documents_with_chinese"] == sum(1 for k in counted if k.endswith(".md"))
    assert mix["chinese_counterpart_documents"] == sum(1 for k in counted if "_CN." in k)
    # the census is not decoration: each counted Chinese manual must sit beside the
    # English document of record it translates.  English files may themselves quote a few
    # CJK characters (sample output, the `_CN` naming rule) without being Chinese
    # documents, so the pairing is checked from the Chinese side only.
    orphan = [
        k for k in counted if k.endswith("_CN.md") and k[: -len("_CN.md")] + ".md" not in entries
    ]
    assert not orphan, f"Chinese documents with no English counterpart in the archive: {orphan}"
    assert mix["chinese_counterpart_documents"] >= 1
    assert "docs/USER_GUIDE_CN.md" in counted, "the note promises a Chinese user manual"
    assert "docs/USER_GUIDE.md" in entries and "docs/USER_GUIDE_CN.md" in entries


@needs_packager
def test_provenance_sidecar_is_complete_and_self_consistent(
    entries: dict[str, bytes], provenance: dict
) -> None:
    version = provenance["package"]["version"]
    source = (ROOT / "src" / "phyloslicer" / "__init__.py").read_text(encoding="utf-8")
    assert f'__version__ = "{version}"' in source, "sidecar version is not the package version"

    control = provenance["source_control"]
    if control.get("git_unavailable"):
        # NV-10, the reviewer's situation: a copy without `.git`.  The packager then
        # makes no commit claim at all, so there is nothing of the right length to
        # assert - and saying so is the point, rather than failing on `len("unknown")`.
        assert control["commit"] is None, "a no-git build must not write a commit string"
        assert "No git metadata" in control["note"]
        assert "b151683" in control["note"], "the release must still be named for editors"
    else:
        assert re.fullmatch(r"[0-9a-f]{40}", str(control["commit"])), (
            f"the sidecar names something that is not a commit: {control['commit']!r}"
        )
        assert "post-release maintenance commit" in control["note"]
        assert "b151683" in control["note"], "the release relationship must be stated for editors"
        # NV-04/RI-08: a dirty build is allowed only when it advertises the dirt.
        if control.get("git_dirty"):
            assert control["git_dirty_paths"], "git_dirty without a path list is a half-truth"
            assert len(control["git_dirty_paths"]) == control["git_dirty_count"]
            assert "uncommitted" in control["note"], "HEAD must be explicitly disclaimed"
            assert _provenance_flags() == ("--allow-dirty",), (
                "the suite built with --allow-dirty but the checkout says otherwise"
            )
        else:
            assert control["git_dirty_paths"] == [], "a clean build must not list paths"
            assert "uncommitted" not in control["note"]

    decision = provenance["anonymisation"]["licence_decision"]
    assert decision["resolved"] is True, "the licence mode may not be shipped as an open question"
    assert decision["decided_mode"] == "placeholder", decision
    assert decision["this_build"] == provenance["anonymisation"]["licence_mode"]
    assert decision["licence"] == "BSD-3-Clause"

    hashed = {item["path"]: item["sha256"] for item in provenance["files"]}
    unhashed = set(entries) - {"REVIEW_ARCHIVE_PROVENANCE.json"}
    assert set(hashed) == unhashed, "every archived file must be hashed in the sidecar"
    for name, digest in hashed.items():
        assert hashlib.sha256(entries[name]).hexdigest() == digest, f"{name} hash drift"

    scrubbed = {item["path"] for item in provenance["anonymisation"]["scrubbed_files"]}
    unreported = set(MUST_BE_SCRUBBED) - scrubbed
    assert not unreported, f"identity-bearing files not reported as scrubbed: {sorted(unreported)}"
    for record in provenance["anonymisation"]["scrubbed_files"]:
        assert record["rules"] and record["occurrences"] > 0
        assert record["path"] in entries
    rules = {item["rule"] for item in provenance["anonymisation"]["rules"]}
    used = {r for item in provenance["anonymisation"]["scrubbed_files"] for r in item["rules"]}
    assert used <= rules, f"scrubbed files cite unknown rules: {used - rules}"
    # CM-20: the ORCID policy is recorded, and no synthetic iD was used for it
    assert "ORCID" in provenance["anonymisation"]["orcid_policy"]
    assert ("0000" + "-0000-0000-000").encode() not in entries["REVIEW_ARCHIVE_PROVENANCE.json"], (
        "the synthetic iD token came back"
    )


@needs_packager
def test_provenance_never_quotes_a_removed_literal(entries: dict[str, bytes]) -> None:
    """A leak in the audit trail would undo the whole exercise."""
    sidecar = entries["REVIEW_ARCHIVE_PROVENANCE.json"].decode("utf-8")
    for pattern, literal in _identity_patterns().items():
        assert literal not in sidecar, f"the sidecar quotes {pattern}"


@needs_packager
def test_two_builds_are_byte_identical(tmp_path: Path) -> None:
    """Determinism, diagnosed: a hash difference names the entry that moved.

    The archive is built from a live working tree, so a mismatch is either the
    packager being non-deterministic or somebody editing the tree between the two
    builds.  Reporting the differing entries tells those apart immediately.
    """
    paths = [build(tmp_path / sub) for sub in ("a", "b")]
    digests = [hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
    if digests[0] != digests[1]:
        first, second = (entries_of(p) for p in paths)
        moved = sorted(
            name for name in set(first) | set(second) if first.get(name) != second.get(name)
        )
        pytest.fail(
            f"builds differ ({digests[0][:12]} vs {digests[1][:12]}); entries that "
            f"changed between the two reads of the tree: {moved}"
        )


@needs_packager
def test_tar_gz_builds_and_passes_the_guard(tmp_path: Path) -> None:
    path = build(tmp_path / "targz", fmt="tar.gz")
    proc = subprocess.run(
        [PY, str(PACKAGER), "--check-only", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr[-2000:]
    assert "identity guard PASSED" in proc.stdout


@needs_packager
def test_tar_gz_is_reproducible_too(tmp_path: Path) -> None:
    digests = [
        hashlib.sha256(build(tmp_path / sub, fmt="tar.gz").read_bytes()).hexdigest()
        for sub in ("a", "b")
    ]
    assert digests[0] == digests[1], f"tar.gz build is not reproducible: {digests}"


@needs_packager
def test_check_only_passes_on_the_real_archive(archive: Path) -> None:
    proc = subprocess.run(
        [PY, str(PACKAGER), "--check-only", str(archive)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


@needs_packager
def test_guard_fails_when_identity_is_planted(tmp_path: Path, archive: Path) -> None:
    """Without this, 'the guard passed' would prove nothing."""
    entries = entries_of(archive)
    entries["pyproject.toml"] = (
        entries["pyproject.toml"]
        + f'\nauthors = [{{ name = "{_identity_patterns()["author_github_login"]}" }}]'.encode()
    )
    tampered = tmp_path / "tampered.zip"
    with zipfile.ZipFile(tampered, "w") as zf:
        for name, data in sorted(entries.items()):
            zf.writestr(f"PhyloSlicer_review_0.1.0/{name}", data)
    proc = subprocess.run(
        [PY, str(PACKAGER), "--check-only", str(tampered)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0, "the guard accepted an archive with an author name in it"
    assert "IDENTITY GUARD FAILED" in proc.stderr
    assert "pyproject.toml" in proc.stderr, proc.stderr[-1500:]
    assert "hash_mismatch" in proc.stderr, "the sidecar hashes must also be checked"


@needs_packager
def test_guard_fails_when_an_artefact_is_planted(tmp_path: Path, archive: Path) -> None:
    entries = entries_of(archive)
    entries["src/phyloslicer/__pycache__/x.cpython-311.pyc"] = b"cached"
    entries[".DS_Store"] = b"metadata"
    tampered = tmp_path / "artefacts.zip"
    with zipfile.ZipFile(tampered, "w") as zf:
        for name, data in sorted(entries.items()):
            zf.writestr(f"PhyloSlicer_review_0.1.0/{name}", data)
    proc = subprocess.run(
        [PY, str(PACKAGER), "--check-only", str(tampered)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert "artefact" in proc.stderr, proc.stderr[-1500:]


@needs_packager
def test_packager_refuses_an_unverified_provenance(tmp_path: Path) -> None:
    """NV-04/RI-08 and NV-10: a false commit claim now stops the build instead of shipping.

    Two states are refused without an explicit flag - uncommitted edits, and no `.git` at
    all - and each becomes legal only by a flag that also makes the sidecar tell the truth
    about it.  Whichever state this checkout is in, the refusal is exercised for it; the
    other is reported as not applicable rather than silently passing.
    """
    state = checkout_state()
    if not state["git"]:
        proc = run_packager("--out", str(tmp_path / "nogit"))
        assert proc.returncode != 0, "a build with no git metadata must not claim a commit"
        assert "PROVENANCE GATE" in proc.stderr, proc.stderr[-1200:]
        assert "--allow-no-git" in proc.stderr
        ok = run_packager("--out", str(tmp_path / "allowed"), "--allow-no-git")
        assert ok.returncode == 0, ok.stderr[-2000:]
        side = json.loads(
            entries_of(next(iter(sorted((tmp_path / "allowed").glob("*.zip")))))[
                "REVIEW_ARCHIVE_PROVENANCE.json"
            ]
        )
        assert side["source_control"]["commit"] is None
        assert side["source_control"]["git_unavailable"] is True
        pytest.skip("no .git here: the dirty-worktree half of the gate cannot be exercised")
    if not state["dirty"]:
        pytest.skip("the worktree is clean, so there is no dirt for the gate to refuse")
    proc = run_packager("--out", str(tmp_path / "refused"))
    assert proc.returncode != 0, "a dirty worktree built against HEAD would misdescribe the archive"
    assert "PROVENANCE GATE" in proc.stderr, proc.stderr[-1500:]
    assert "--allow-dirty" in proc.stderr
    assert not (tmp_path / "refused").exists(), "the refused build must not write an archive"
    allowed = run_packager("--out", str(tmp_path / "allowed"), "--allow-dirty")
    assert allowed.returncode == 0, allowed.stderr[-2000:]
    side = json.loads(
        entries_of(next(iter(sorted((tmp_path / "allowed").glob("*.zip")))))[
            "REVIEW_ARCHIVE_PROVENANCE.json"
        ]
    )["source_control"]
    assert side["commit"] and re.fullmatch(r"[0-9a-f]{40}", str(side["commit"]))
    assert side["git_dirty"] is True and side["git_dirty_count"] == len(side["git_dirty_paths"])
    assert side["git_dirty_paths"], "--allow-dirty must record the uncommitted paths"
    assert "uncommitted" in side["note"]


@needs_packager
def test_no_unresolved_action_marker_survives(entries: dict[str, bytes], provenance: dict) -> None:
    """RI-18: the archive must not advertise a decision its author has not made."""
    for name, data in sorted(entries.items()):
        for marker in UNRESOLVED_MARKERS:
            assert marker.encode() not in data, f"{name} still carries {marker!r}"
    decision = provenance["anonymisation"]["licence_decision"]
    assert decision["resolved"] is True and decision["decided_mode"] == "placeholder"
    assert decision["this_build"] == provenance["anonymisation"]["licence_mode"]
    assert decision["licence"] == "BSD-3-Clause"
    assert decision["deviation"] is None, decision["deviation"]
    # the maintainer-side notes are where the marker used to live.  They may *name* the
    # forbidden phrases - that is how they explain the rule - but no line of them may
    # *be* one, so the directive form (**marker**) and the Chinese action line are gone
    for name in ("REVIEW_ARCHIVE.md", "REVIEW_ARCHIVE_CN.md"):
        path = ROOT / name
        if path.is_file():
            text = path.read_text(encoding="utf-8")
            for marker in UNRESOLVED_MARKERS:
                assert f"**{marker}**" not in text, f"{name} still opens an action item"
            for directive in ("需作者确认", "提请作者确认"):
                assert directive not in text, f"{name} still asks the author to confirm"


@needs_packager
def test_guard_fails_when_invalid_orcid_is_planted(tmp_path: Path, archive: Path) -> None:
    """CM-20: a well-formed-looking iD that resolves nowhere is rejected, not shipped.

    The planted value is the synthetic token an earlier build used as its replacement
    (with a deliberately wrong check character), so this case is the regression that keeps
    the anonymisation from re-inventing fake identifiers.
    """
    entries = entries_of(archive)
    tampered = tmp_path / "orcid.zip"
    # assembled from fragments: this module ships, and a file that carries the fake iD it
    # plants would itself be the leak `orcid_violations()` is there to catch
    planted_id = "0000" + "-0000-0000" + "-0009"
    cff = entries["CITATION.cff"].decode("utf-8")
    assert "WITHHELD_FOR_PEER_REVIEW" in cff, "the shipped metadata must not carry an iD at all"
    planted = cff.replace(
        "orcid: https://orcid.org/WITHHELD_FOR_PEER_REVIEW",
        f"orcid: https://orcid.org/{planted_id}",
        1,
    )
    assert planted != cff, "the CITATION.cff orcid line moved: the plant found nothing to replace"
    with zipfile.ZipFile(tampered, "w") as zf:
        for name, data in entries.items():
            blob = planted.encode("utf-8") if name == "CITATION.cff" else data
            zf.writestr(f"PhyloSlicer_review_0.1.0/{name}", blob)
    proc = run_packager("--check-only", str(tampered))
    assert proc.returncode != 0, "the guard accepted a fabricated ORCID iD"
    assert "invalid_orcid" in proc.stderr, proc.stderr[-1500:]
    assert planted_id in proc.stderr, proc.stderr[-1500:]


@needs_packager
def test_guard_fails_when_resolvable_zenodo_doi_is_planted(tmp_path: Path, archive: Path) -> None:
    """RI-09/CM-02: the version DOI is the identifier that de-anonymises the archive.

    One ``doi.org`` click from inside a "double-blind" package reaches the Zenodo record
    that names both authors, their affiliation and their repository, so the guard now
    rejects *any* ``10.5281/zenodo.<digits>`` run rather than only the two known ids.
    """
    entries = entries_of(archive)
    cff = entries["CITATION.cff"].decode("utf-8")
    token_line = "10.5281/zenodo." + "WITHHELD_FOR_PEER_REVIEW"
    assert token_line in cff, "the shipped metadata lost its withheld DOI line"
    # the digits come from the excluded tools/ table, never spelled in this shipped file
    version_id = _identity_patterns()["zenodo_version_record"]
    planted = cff.replace(token_line, f"10.5281/zenodo.{version_id}")
    assert planted != cff, "the plant found nothing to replace"
    assert version_id not in cff, (
        "the version record id is already sitting in the shipped CITATION.cff"
    )
    tampered = tmp_path / "doi.zip"
    with zipfile.ZipFile(tampered, "w") as zf:
        for name, data in entries.items():
            blob = planted.encode("utf-8") if name == "CITATION.cff" else data
            zf.writestr(f"PhyloSlicer_review_0.1.0/{name}", blob)
    proc = run_packager("--check-only", str(tampered))
    assert proc.returncode != 0, "the guard shipped an archive with a live Zenodo DOI"
    assert "resolvable_zenodo_doi" in proc.stderr, proc.stderr[-1500:]
    assert "CITATION.cff" in proc.stderr, proc.stderr[-1500:]


@needs_packager
def test_guard_fails_when_an_identity_is_reassembled_from_fragments(
    tmp_path: Path, archive: Path
) -> None:
    """RI-17/CM-03: splitting a withheld literal no longer hides it - it fails.

    The archived test file used to carry the whole identity table spelled as
    ``_frag("Ze", "ng")`` so the guard could not see it, which also let a reader
    reassemble it.  The table lives in ``tools/`` now, and this case pins that the
    same trick in *any* shipped ``.py`` entry is a guard failure, not an evasion.
    """
    surname = _identity_patterns()["author_surname"]
    cut = len(surname) // 2
    planted_source = (
        "def _frag(*parts):\n"
        '    return "".join(parts)\n\n'
        f'TABLE = {{"who": _frag("{surname[:cut]}", "{surname[cut:]}")}}\n'
    )
    assert surname not in planted_source, "the plant spelled the literal contiguously"
    entries = entries_of(archive)
    entries["src/phyloslicer/_planted_identity.py"] = planted_source.encode("utf-8")
    tampered = tmp_path / "frag.zip"
    with zipfile.ZipFile(tampered, "w") as zf:
        for name, data in sorted(entries.items()):
            zf.writestr(f"PhyloSlicer_review_0.1.0/{name}", data)
    proc = run_packager("--check-only", str(tampered))
    assert proc.returncode != 0, "the guard accepted a fragment-joined author surname"
    assert "reassembled_identity" in proc.stderr, proc.stderr[-1500:]
    assert "_planted_identity.py" in proc.stderr, proc.stderr[-1500:]


def _extract(archive: Path, target: Path) -> Path:
    with zipfile.ZipFile(archive) as zf:
        zf.extractall(target)
    return next(iter(p for p in target.iterdir() if p.is_dir()))


@pytest.fixture(scope="module")
def unpacked(archive: Path, tmp_path_factory) -> Path:
    target = tmp_path_factory.mktemp("unpacked")
    return _extract(archive, target)


def _run(
    cwd: Path, argv: list[str], timeout: int, *, with_src_path: bool = False
) -> tuple[int, str, float]:
    started = time.monotonic()
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}
    env.pop("PYTHONPATH", None)
    if with_src_path:
        # what pyproject.toml's `pythonpath = ["src"]` does for pytest, spelled
        # out for the scripts that do not insert src/ themselves
        env["PYTHONPATH"] = str(cwd / "src")
    proc = subprocess.run(
        [PY, *argv],
        cwd=cwd,
        capture_output=True,
        text=True,
        timeout=timeout,
        env=env,
        check=False,
    )
    return proc.returncode, proc.stdout + proc.stderr, time.monotonic() - started


@needs_packager
@pytest.mark.skipif(os.environ.get(SKIP_RUN_ENV) == "1", reason=f"{SKIP_RUN_ENV}=1")
def test_reviewer_can_run_the_test_suite(unpacked: Path) -> None:
    """The 'a reviewer can run it' claim, machine-checked."""
    rc, output, elapsed = _run(unpacked, ["-m", "pytest", "-q", "-p", "no:cacheprovider"], 900)
    print(f"\nin-archive pytest: rc={rc} {elapsed:.1f}s\n{output[-2500:]}")
    assert rc == 0, f"the archived suite does not pass:\n{output[-4000:]}"


@needs_packager
@pytest.mark.skipif(os.environ.get(SKIP_RUN_ENV) == "1", reason=f"{SKIP_RUN_ENV}=1")
def test_reviewer_can_run_the_documented_example(unpacked: Path) -> None:
    rc, output, elapsed = _run(
        unpacked,
        ["examples/quickstart.py", "--sites", "24", "--species", "60", "--trees", "3"],
        300,
    )
    print(f"\nin-archive quickstart: rc={rc} {elapsed:.1f}s\n{output[:900]}")
    assert rc == 0, output[-3000:]
    assert "converged" in output and "rate_median" in output


@needs_packager
@pytest.mark.skipif(os.environ.get(SKIP_RUN_ENV) == "1", reason=f"{SKIP_RUN_ENV}=1")
def test_reviewer_can_rerun_the_cross_validation(unpacked: Path, tmp_path: Path) -> None:
    """The hash-registered-reference claim the manuscript rests on.

    ``validation/compare.py`` — unlike ``examples/quickstart.py`` and the
    benchmark drivers — does not put ``src/`` on ``sys.path`` itself, so a
    reviewer runs it either with ``PYTHONPATH=src`` (what the user guide's
    "run everything" block implies for a source checkout) or after installing
    the package.  Both are exercised: here the no-install form, and in
    :func:`test_reviewer_can_install_the_archived_package` the installed form.
    """
    manifest = unpacked / "validation" / "reference" / "manifest.txt"
    assert manifest.is_file(), "the reference manifest must ship for the comparator to run"
    rc, output, elapsed = _run(unpacked, ["validation/compare.py"], 900, with_src_path=True)
    print(f"\nin-archive validation/compare.py: rc={rc} {elapsed:.1f}s\n{output[-1200:]}")
    assert rc == 0, output[-4000:]
    assert "verdict: PASS" in output, output[-1500:]
    # the manifest is not decoration: a reference it does not register has to stop
    # the comparator, which is what makes the hashes a claim rather than a list.
    # Run against a copy so the shared unpacked tree stays pristine for later tests.
    work = tmp_path / "reference"
    shutil.copytree(unpacked / "validation" / "reference", work)
    (work / "bogus__unregistered.csv").write_text("site_id,value\ns1,1.0\n", encoding="utf-8")
    rc2, output2, _ = _run(unpacked, ["validation/compare.py", str(work)], 900, with_src_path=True)
    assert rc2 != 0, "an unregistered reference file must stop the comparator"
    print(f"unregistered reference refused as expected: {output2[-300:]}")


@needs_packager
@pytest.mark.skipif(os.environ.get(SKIP_RUN_ENV) == "1", reason=f"{SKIP_RUN_ENV}=1")
@pytest.mark.skipif(
    subprocess.run([PY, "-c", "import setuptools"], capture_output=True).returncode != 0,
    reason="setuptools is needed to build the archived package offline",
)
def test_reviewer_can_install_the_archived_package(unpacked: Path, tmp_path: Path) -> None:
    """``pip install`` of the archive alone must work, with no author metadata."""
    target = tmp_path / "site"
    proc = subprocess.run(
        [
            PY,
            "-m",
            "pip",
            "install",
            "--quiet",
            "--no-deps",
            "--no-build-isolation",
            "--target",
            str(target),
            ".",
        ],
        cwd=unpacked,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    env = {
        **os.environ,
        "PYTHONPATH": str(target),  # pip --target unpacks packages at this root
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    run = subprocess.run(
        [PY, "-c", "import phyloslicer; print(phyloslicer.__version__)"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    assert run.returncode == 0, run.stderr[-2000:]
    source = (ROOT / "src" / "phyloslicer" / "__init__.py").read_text(encoding="utf-8")
    assert f'__version__ = "{run.stdout.strip()}"' in source


@needs_packager
def test_license_preserved_mode_ships_the_notice_verbatim(tmp_path: Path) -> None:
    """The escape hatch for the licence decision, and its documented cost.

    ``--license preserved`` is the alternative the author may choose instead of
    tokenising the copyright holder (see REVIEW_ARCHIVE.md).  It must ship
    ``LICENSE`` byte-for-byte - a half-substituted notice would name nobody and
    license nothing - and it must still scrub everything else, so checking the
    result in the default mode has to fail rather than quietly pass.
    """
    proc = run_packager("--out", str(tmp_path), "--license", "preserved", *_provenance_flags())
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "copyright holder" in proc.stderr, "the exception must be announced, not silent"
    kept = entries_of(next(iter(sorted(tmp_path.glob("*.zip")))))
    assert kept["LICENSE"] == (ROOT / "LICENSE").read_bytes(), "the notice was altered"
    # everything else is still scrubbed, so the guard has to reject this archive
    # unless the caller declares the licence exception
    strict = subprocess.run(
        [PY, str(PACKAGER), "--check-only", str(next(iter(sorted(tmp_path.glob("*.zip")))))],
        capture_output=True,
        text=True,
        check=False,
    )
    assert strict.returncode != 0, "an archive that names the holder must not pass by default"
    assert "LICENSE" in strict.stderr, strict.stderr[-1200:]


@needs_packager
def test_source_files_were_not_modified() -> None:
    """Anonymisation writes into the archive, never into the working tree."""
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert _identity_patterns()["author_given_name"] in license_text, (
        "the source licence must keep naming the copyright holder"
    )
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert _identity_patterns()["author_email"] in pyproject
    assert (ROOT / "REVIEW_ARCHIVE.md").is_file(), "the archive recipe must be documented"
