#!/usr/bin/env python
"""Build the anonymised source archive that the Ecography editors forward to reviewers.

Ecography Software Notes allow a submission to be reviewed through a
supplementary archive instead of an anonymous repository mirror.  That archive
has to satisfy two properties at once, and they pull against each other:

* it must be *runnable* - a reviewer unpacks it and immediately runs the test
  suite, the documented example and the cross-validation comparator;
* it must carry *no author identity* - not in the documents, not in the
  metadata, not in a cache directory, not in an absolute local path.

This tool produces that archive and then proves it did: every entry it writes
is re-read and scanned, and if a single identifying pattern survives the build
it prints the offending entry, the pattern and the surrounding text and exits
non-zero.  ``--check-only`` runs the same guard over an archive that already
exists, so the guard can be exercised without rebuilding anything.

Two design decisions are load-bearing and explained in ``REVIEW_ARCHIVE.md``:

1. Scrubbing is a *uniform token substitution* applied to every archived text
   file, including the tests that assert on the metadata.  Because
   ``tests/test_release_hygiene.py`` demands that ``pyproject.toml``,
   ``CITATION.cff``, ``.zenodo.json``, ``README*.md`` and ``LICENSE`` all name
   the same author and the same repository, replacing the literals *everywhere*
   keeps those consistency assertions true in the anonymised copy: the reviewer
   gets a self-consistent package rather than a package whose own hygiene test
   fails.  The source files in the working tree are never modified.
2. ``tools/`` itself is *not* shipped inside the archive, because the scrub
   table here necessarily contains the literals it removes; shipping it would
   undo the anonymisation.  ``REVIEW_ARCHIVE.md`` is excluded for the same
   reason.  ``tests/test_review_archive.py`` ships, and skips itself when there
   is no packager to drive.  Its scrub table used to sit inside that shipped
   file, spelled in fragments so the guard could not see it - which also meant a
   reviewer could reassemble it (RI-17/CM-03); the literals now live in
   ``tools/identity_literals.py``, which is excluded with this module, and
   ``fragment_violations()`` rejects any archived source file that reassembles a
   withheld literal out of fragment joins, so the trick cannot be repeated.

Usage::

    python tools/make_review_archive.py --verbose
    python tools/make_review_archive.py --format tar.gz --out dist/review_archive
    python tools/make_review_archive.py --check-only \\
        dist/review_archive/PhyloSlicer_review_0.1.0.zip

No ``PYTHONPATH`` is needed: the tool reads the version out of the package
source instead of importing it, so it neither needs the runtime dependencies nor
a built extension.  Stdlib only, because build-time tooling must not add a
runtime dependency.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import getpass
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import tarfile
import time
from dataclasses import dataclass, field
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

TOOL_VERSION = "1.0.0"
ROOT = Path(__file__).resolve().parents[1]

#: Top-level documents that must ship: the reviewer needs the licence, the
#: changelog and the readme, and ``test_release_hygiene`` additionally reads
#: every file in its ``VERSION_BEARING`` list plus ``.gitignore`` itself, so an
#: archive that dropped one of them would fail its own test suite.
INCLUDE_FILES = (
    "pyproject.toml",
    "MANIFEST.in",
    ".gitignore",
    "LICENSE",
    "NOTICE.md",
    "NOTICE_CN.md",
    "README.md",
    "README_CN.md",
    "CHANGELOG.md",
    "CHANGELOG_CN.md",
    "CONTRIBUTING.md",
    "CONTRIBUTING_CN.md",
    "RELEASE.md",
    "RELEASE_CN.md",
    "CODE_OF_CONDUCT.md",
    "CODE_OF_CONDUCT_CN.md",
    "SECURITY.md",
    "SECURITY_CN.md",
    "CITATION.cff",
    ".zenodo.json",
)

#: Source trees a reviewer needs in order to run and re-verify the package.
INCLUDE_DIRS = ("src", "tests", "examples", "docs", "benchmarks", "validation", "testing_suite")

#: Top-level entries that are deliberately left out, with the reason recorded in
#: the build summary and in the provenance sidecar.
EXCLUDED_TOPLEVEL = {
    # build / VCS / editor artefacts
    ".git": "version-control metadata: commit author names and e-mails, local paths",
    ".github": "repository automation and .github/CODEOWNERS, which names the owner",
    ".pre-commit-config.yaml": "commit-time tooling; nothing in it is needed to run the package",
    ".ruff_cache": "tool cache",
    ".pytest_cache": "tool cache",
    ".mypy_cache": "tool cache",
    "build": "build output",
    "dist": "build output",
    # author-side documentation
    "tools": "build-time packager: its scrub table necessarily quotes the literals it removes",
    "REVIEW_ARCHIVE.md": "documents the anonymisation and therefore quotes the literals it removes",
    "REVIEW_ARCHIVE_CN.md": "Chinese counterpart of the above, same reason",
}

#: Directory names pruned wherever they appear.
PRUNE_DIRS = {
    ".git",
    ".github",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    ".tox",
    ".venv",
    "venv",
    "node_modules",
    "build",
    "dist",
    ".eggs",
    "tools",
}

#: File names pruned wherever they appear (macOS/Windows metadata, caches).
PRUNE_NAMES = {".DS_Store", "Thumbs.db", "desktop.ini"}

#: Byte-code and install-metadata suffixes.
PRUNE_SUFFIXES = (".pyc", ".pyo", ".egg-info")

#: Text suffixes that ``LICENSE_MODE`` and the numeric-pattern scope key off.
DOC_SUFFIXES = {
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
}

#: Timestamp used when neither ``SOURCE_DATE_EPOCH`` nor git is available.
FALLBACK_EPOCH = int(dt.datetime(2026, 1, 1, tzinfo=dt.timezone.utc).timestamp())

#: The tag the deposited Zenodo record was cut from, and the release-preparation
#: commit named in the submission checklist.  Both are recorded in the provenance
#: sidecar so an editor can reconcile this archive against the deposited version.
RELEASED_TAG = "v0.1.0"
RELEASE_PREP_COMMIT = "b1516835caa7cc227fbac5c407fda278efe5b499"

#: RI-18: the licence-mode question is closed, not shipped as an open action.  The
#: submitted archive is the ``placeholder`` build of the BSD-3-Clause licence: the
#: notice ships in full - all three conditions and the disclaimer verbatim - and only
#: the copyright-holder names (``Zichao Zeng and Yinzhao Wang``, per CM-21) are
#: tokenised.  ``preserved`` remains available as an explicit alternative, but a build
#: that is not in the recorded mode is labelled a deviation in the sidecar, and the
#: guard refuses any archive whose sidecar leaves the decision unresolved.
RESOLVED_LICENSE_MODE = "placeholder"
RESOLVED_LICENSE_STATEMENT = (
    "BSD-3-Clause, recorded decision: this archive is the 'placeholder' build, in which "
    "LICENSE carries the tokens 'Author One Withheld' / 'Author Two Withheld' in place of "
    "the two copyright holders and nothing else in the notice is altered.  The deposited "
    "and public versions carry both names verbatim."
)

#: Phrases that would make the archive advertise an unfinished authorial decision.
#: The guard rejects them wherever they appear, so the marker cannot come back.
UNRESOLVED_ACTION_MARKERS = (
    "AUTHOR ACTION",
    "confirm which mode",
    "TODO: confirm",
    "TBD: confirm",
)

#: NV-10 / NV-04: the two provenance states the packager refuses to build in unless the
#: author says otherwise, spelled out rather than degraded into a fake value.
GIT_ABSENT_MESSAGE = (
    "PROVENANCE GATE: no readable git metadata (no `.git`, or no `git` on PATH), so no "
    "commit can be named.  Refusing to write an archive that would claim one.  Re-run "
    "with --allow-no-git to build an archive whose sidecar records commit = null and "
    "git_unavailable = true (no commit claim at all), or run this tool from a checkout."
)
DIRTY_MESSAGE = (
    "PROVENANCE GATE: the working tree carries uncommitted edits, so naming HEAD alone "
    "would misdescribe the bytes being shipped.  Refusing to build.  Commit them, or "
    "re-run with --allow-dirty, which records every uncommitted path in the sidecar under "
    "source_control.git_dirty_paths and says so in source_control.note."
)

# --------------------------------------------------------------------------
# identity rules
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Rule:
    """One literal author-identifying string and the neutral token replacing it."""

    name: str
    description: str
    literal: str
    replacement: str
    #: Numeric identifiers (grant numbers) are only scanned in documents, code
    #: and configuration: an eight-digit run inside a reference table's float is
    #: a normal artefact, not a funding number.
    doc_only: bool = False

    def apply(self, text: str) -> tuple[str, int]:
        """Replace every occurrence of the literal; return the text and the count."""
        n = text.count(self.literal)
        if not n:
            return text, 0
        return text.replace(self.literal, self.replacement), n


#: CM-20: the token that replaces an ORCID iD.  Deliberately not iD-shaped, so a
#: metadata harvester cannot mistake it for a registered identifier.
WITHHELD_ORCID_TOKEN = "WITHHELD_FOR_PEER_REVIEW"

#: Static rules, longest literal first so that ``Zichao Zeng`` is never left as a
#: fragment of an earlier ``Zeng`` replacement.  These are the strings the audit
#: found in the tree; see REVIEW_ARCHIVE.md for the file:line evidence.
STATIC_RULES: tuple[Rule, ...] = (
    Rule(
        "author_name_cn_block",
        "first author, Chinese attribution line",
        "曾子超（Zichao Zeng）",
        "第一作者（同行评审期间隐去）",
    ),
    Rule(
        "author_name_en_block",
        "first author, English attribution line with Chinese name",
        "Zichao Zeng (曾子超)",
        "Author One Withheld",
    ),
    Rule(
        "affiliation_full",
        "institution string with school",
        "School of Life Sciences and Biotechnology, Shanghai Jiao Tong University",
        "Institution Withheld For Peer Review",
    ),
    Rule(
        "affiliation_cn",
        "institution, Chinese",
        "上海交通大学生命科学技术学院",
        "所属机构（同行评审期间隐去）",
    ),
    Rule(
        "affiliation_cn_short",
        "institution, Chinese, short form",
        "上海交通大学",
        "机构（同行评审期间隐去）",
    ),
    Rule(
        "author_family_given",
        "first author, family-given order used by .zenodo.json",
        "Zeng, Zichao",
        "Withheld, Author One",
    ),
    Rule(
        "coauthor_family_given",
        "co-author, family-given order used by .zenodo.json",
        "Wang, Yinzhao",
        "Withheld, Author Two",
    ),
    Rule(
        "author_name_given_family",
        "first author, given-family order used by pyproject/CITATION/LICENSE",
        "Zichao Zeng",
        "Author One Withheld",
    ),
    #: CM-21: ``LICENSE`` now names both paper authors as copyright holders
    #: (``Copyright (c) 2026, Zichao Zeng and Yinzhao Wang``), so the bare surname of
    #: the second one would otherwise survive the scrub as a scan-only hit
    #: (``surname_2``).  Tokenise the pair as a unit, in given-family order.
    Rule(
        "coauthor_name_given_family",
        "co-author, given-family order used by the LICENSE copyright line",
        "Yinzhao Wang",
        "Author Two Withheld",
    ),
    Rule(
        "author_github_login",
        "GitHub owner of the public repository",
        "ZengZichao",
        "AUTHOR_WITHHELD_FOR_PEER_REVIEW",
    ),
    Rule(
        "author_email",
        "first author e-mail",
        "zengzichao@sjtu.edu.cn",
        "author.one.withheld@example.invalid",
    ),
    Rule(
        "coauthor_email",
        "co-author e-mail",
        "wyz@sjtu.edu.cn",
        "author.two.withheld@example.invalid",
    ),
    Rule(
        "zenodo_doi",
        "Zenodo record DOI, which resolves to the named authors",
        "10.5281/zenodo.23057620",
        "10.5281/zenodo.WITHHELD_FOR_PEER_REVIEW",
    ),
    Rule(
        "zenodo_version_doi",
        "Zenodo version DOI of the release under review, which resolves to the "
        "named authors and their affiliation on the record page",
        "10.5281/zenodo.23058460",
        "10.5281/zenodo.WITHHELD_FOR_PEER_REVIEW",
    ),
    #: CM-20: the replacement must not look like an issued identifier.  An earlier
    #: build wrote the synthetic iDs ``0000-0000-0000-0001/2``, which are
    #: well-formed enough to be mistaken for a deposit error (and, since the all-zero
    #: stem happens to carry a valid ISO 7064 check digit, a checksum test would not
    #: even catch them).  The policy is now: no iD at all - the field keeps its key
    #: and carries the same neutral token as the rest of the scrub table, and
    #: ``orcid_violations()`` rejects any ORCID-shaped string in the archive that does
    #: not satisfy the checksum rule, so neither a synthetic nor a real iD can ship.
    Rule(
        "author_orcid",
        "first author ORCID",
        "0000-0001-6553-970X",
        WITHHELD_ORCID_TOKEN,
    ),
    Rule(
        "coauthor_orcid",
        "co-author ORCID",
        "0000-0001-8715-1386",
        WITHHELD_ORCID_TOKEN,
    ),
    Rule("given_name_1", "first author given name", "Zichao", "AuthorOne"),
    Rule("given_name_2", "co-author given name", "Yinzhao", "AuthorTwo"),
    Rule("cjk_name_1", "first author name in Chinese characters", "曾子超", "（已隐去）"),
    Rule(
        "family_name_1",
        "CITATION.cff family-names field, first author",
        "family-names: Zeng",
        "family-names: Withheld",
    ),
    Rule(
        "family_name_2",
        "CITATION.cff family-names field, co-author",
        "family-names: Wang",
        "family-names: Withheld",
    ),
    Rule(
        "institution",
        "institution name",
        "Shanghai Jiao Tong University",
        "Institution Withheld For Peer Review",
    ),
    # README.md wraps the institution across a line break, so the full-string rule
    # above cannot reach the second half; without these two the city and the tail
    # of the university name would have survived anonymisation in README.md:144.
    Rule(
        "institution_tail",
        "institution name, wrapped across lines",
        "Jiao Tong University",
        "Institution Withheld For Peer Review",
    ),
    Rule("institution_fragment", "institution name fragment", "Jiao Tong", "Withheld"),
    Rule("affiliation_city", "affiliation city", "Shanghai", "Withheld"),
    Rule(
        "school",
        "school name",
        "School of Life Sciences and Biotechnology",
        "School Withheld For Peer Review",
    ),
    Rule("email_domain", "institutional e-mail domain", "sjtu.edu.cn", "example.invalid"),
    Rule("grant_1", "NSFC grant number", "42422209", "GRANT_WITHHELD_FOR_PEER_REVIEW", True),
    Rule("grant_2", "NSFC grant number", "42272354", "GRANT_WITHHELD_FOR_PEER_REVIEW", True),
    Rule(
        "grant_3",
        "national R&D programme number",
        "2023YFC3108600",
        "GRANT_WITHHELD_FOR_PEER_REVIEW",
        True,
    ),
)

#: Patterns the guard rejects but no rule rewrites: they are either incidental
#: (a cache directory that should never have been collected) or third-party text
#: that must survive untouched, so the right fix is to exclude or to edit, not to
#: substitute.  ``Wang`` and ``Zeng`` are only safe as scan patterns because the
#: whole tree was grepped and each occurs solely in the author metadata.
#: The short literals are assembled from parts: this file is a ``.py`` in the
#: tree, so writing a home-directory prefix out in full would trip
#: ``test_release_hygiene.py`` and its personal-path sweep - the same reason that
#: test spells the prefix in pieces.  (An earlier revision of this comment
#: illustrated the point by committing the very prefix it warns about, and the
#: source-tree suite failed on it: the guard earns its keep in both directions.)
SCAN_ONLY_LITERALS: tuple[tuple[str, str, str, bool], ...] = (
    ("home_users_dir", "any macOS home path", "/" + "Users" + "/", False),
    ("home_windows_dir", "any Windows home path", "C:" + "\\" + "Users", True),
    ("windows_unc", "Windows home path, alternate spelling", "C:/" + "Users", False),
    ("github_hostname", "a GitHub URL of any kind", "github.com/" + "Zeng", False),
    ("orcid_scheme", "an ORCID iD of any kind", "orcid.org/" + "0000-0001", False),
    (
        "synthetic_orcid",
        "a synthetic all-zero ORCID iD used as an anonymisation token (CM-20)",
        "0000" + "-0000-0000-000",
        False,
    ),
    ("surname_1", "first author surname", "Ze" + "ng", False),
    ("surname_2", "co-author surname", "Wa" + "ng", False),
    ("zenodo_record", "the Zenodo record id", "2305" + "7620", False),
    ("zenodo_version_record", "the Zenodo version-record id", "2305" + "8460", False),
)


@dataclass
class Hit:
    """A surviving identifying pattern."""

    where: str
    pattern: str
    description: str
    context: str


def rules_for(root: Path) -> list[Rule]:
    """Static rules plus the machine-local paths, longest literal first.

    The home directory, this checkout's path and the local account name are
    identity for a reviewer: they name the author's machine.  They are derived
    at build time rather than hard-coded, so the tool works on any machine and
    still fails if such a path is ever committed.
    """
    home = str(Path.home())
    local = [
        Rule(
            "checkout_path", "absolute path of this working tree", str(root), "/review/PhyloSlicer"
        ),
        Rule("home_directory", "absolute home directory of the builder", home, "/home/reviewer"),
        Rule("local_account", "local account name", getpass.getuser(), "reviewer"),
    ]
    return sorted((*STATIC_RULES, *local), key=lambda r: -len(r.literal))


def scan_patterns(rules: list[Rule]) -> list[tuple[str, str, str, bool]]:
    """(name, description, literal, doc-only) tuples the guard rejects anywhere.

    Every scrub rule's own literal is re-added as a pattern, so a rule that did
    not fire - because a new file introduced the string in a form the rule does
    not cover - is caught rather than shipped.
    """
    from_rules = [(r.name, r.description, r.literal, r.doc_only) for r in rules]
    return from_rules + [(n, d, lit, doc) for n, d, lit, doc in SCAN_ONLY_LITERALS]


# --------------------------------------------------------------------------
# provenance helpers
# --------------------------------------------------------------------------


def git(*args: str) -> str | None:
    """Run a read-only git command; ``None`` when git or the repo is absent."""
    try:
        proc = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    except OSError:
        return None
    if proc.returncode:
        return None
    return proc.stdout.strip() or None


def git_raw(*args: str) -> str | None:
    """Like :func:`git`, but an empty success is ``''`` rather than ``None``.

    ``git status --porcelain`` answers with an empty string exactly when the tree is
    clean, which ``git()`` cannot tell apart from "no repository here".
    """
    try:
        proc = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=False)
    except OSError:
        return None
    return proc.stdout if proc.returncode == 0 else None


def worktree_state() -> dict[str, object]:
    """The raw git facts about the checkout the archive is read from."""
    porcelain = git_raw("status", "--porcelain")
    return {
        "head": git("rev-parse", "HEAD"),
        "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
        "cleanly_queried": porcelain is not None,
        "dirty_paths": sorted(
            {
                line[3:].strip().split(" -> ")[-1]
                for line in (porcelain or "").splitlines()
                if line.strip()
            }
        ),
    }


def enforce_provenance_policy(
    state: dict[str, object], *, allow_dirty: bool, allow_no_git: bool
) -> None:
    """Refuse to build on a provenance the archive could not tell truthfully (NV-04/RI-08).

    Two states are gated: no git metadata at all, and a worktree with uncommitted
    edits.  Both used to pass silently - the first by writing the literal ``"unknown"``
    into a field the suite then measured for length, the second by naming HEAD as if it
    described the shipped bytes while dozens of paths differed from it.  Both now stop
    the build, and both are survivable only by an explicit flag that also makes the
    sidecar say so.
    """
    if state["head"] is None:
        if not allow_no_git:
            raise SystemExit(GIT_ABSENT_MESSAGE)
        return
    dirty = state["dirty_paths"]
    if dirty and not allow_dirty:
        listed = "\n".join(f"    {p}" for p in list(dirty)[:10])  # type: ignore[index]
        more = len(dirty) - 10  # type: ignore[arg-type]
        raise SystemExit(
            f"{DIRTY_MESSAGE}\n  {len(dirty)} uncommitted path(s):"  # type: ignore[arg-type]
            f"\n{listed}" + (f"\n    … and {more} more" if more > 0 else "")
        )


def _ancestor(sha: str) -> bool | None:
    """Is ``sha`` an ancestor of HEAD?  ``None`` when git cannot say."""
    try:
        proc = subprocess.run(
            ["git", "merge-base", "--is-ancestor", sha, "HEAD"],
            cwd=ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    if proc.returncode not in (0, 1):
        return None
    return proc.returncode == 0


def package_version(root: Path) -> str:
    """The version string, read from the package source rather than imported.

    ``src/phyloslicer/__init__.py`` is the single source of truth declared in
    ``pyproject.toml``; importing it here would drag numpy/pandas into a
    build-time tool.
    """
    text = (root / "src" / "phyloslicer" / "__init__.py").read_text(encoding="utf-8")
    found = re.search(r'^__version__\s*=\s*["\']([^"\']+)["\']', text, re.M)
    if not found:
        raise SystemExit(f"cannot read __version__ from {root / 'src/phyloslicer/__init__.py'}")
    return found.group(1)


def fixed_epoch() -> tuple[int, str]:
    """The single timestamp every entry gets, so two builds are byte-identical."""
    env = os.environ.get("SOURCE_DATE_EPOCH")
    if env:
        return int(env), "SOURCE_DATE_EPOCH"
    commit_epoch = git("show", "-s", "--format=%ct", "HEAD")
    if commit_epoch:
        return int(commit_epoch), "commit date of HEAD"
    return FALLBACK_EPOCH, "fixed fallback timestamp (no git, no SOURCE_DATE_EPOCH)"


def source_control(state: dict[str, object]) -> dict[str, object]:
    """Commit identity of the working tree, and where it sits relative to the release.

    The released DOI in the public metadata points at the commit Zenodo archived,
    not at this checkout, so an editor reconciling the archive against the
    deposited version needs both named.  A commit hash identifies nobody on its
    own: while the repository is private, only someone who already has access to
    its history can resolve one to an author.

    NV-10: when there is no git metadata, the commit field is ``null`` and
    ``git_unavailable`` is true - an absent claim, not a fake 7-character one.
    NV-04/RI-08: when the build was allowed on a dirty tree, the uncommitted paths are
    listed and HEAD is explicitly disclaimed as a full description of the archive.
    """
    head = state["head"]
    branch = state["branch"]
    dirty = list(state["dirty_paths"] or [])  # type: ignore[arg-type]
    if head is None:
        return {
            "git_available": False,
            "git_unavailable": True,
            "commit": None,
            "commit_format": "none: no commit is claimed",
            "branch": None,
            "describe": None,
            "released_tag": RELEASED_TAG,
            "released_tag_commit": None,
            "commits_after_released_tag": None,
            "release_preparation_commit": RELEASE_PREP_COMMIT,
            "git_dirty": None,
            "git_dirty_count": None,
            "git_dirty_paths": [],
            "note": (
                "No git metadata was readable when this archive was built (no `.git` "
                "directory, or no `git` executable), and the build was explicitly allowed "
                "with --allow-no-git.  `commit` is therefore null: this archive makes no "
                "claim about a source commit, and the provenance cases of "
                "tests/test_review_archive.py skip the commit assertions rather than fail.  "
                "Reconcile against the released tag v0.1.0 and the release-preparation "
                f"commit {RELEASE_PREP_COMMIT[:7]} instead."
            ),
        }
    released = git("rev-parse", "v0.1.0^{commit}")
    after = git("rev-list", "--count", "v0.1.0..HEAD")
    known = _ancestor(RELEASE_PREP_COMMIT)
    reachable = {True: "an ancestor of", False: "not an ancestor of", None: "unresolvable in"}[
        known
    ]
    prep_subject = git("log", "-1", "--format=%s", str(RELEASE_PREP_COMMIT))
    prep_subject = prep_subject or "subject unavailable"
    note = (
        f"Built from post-release maintenance commit {str(head)[:7]} on '{branch}', "
        f"{after} commits after the released tag {RELEASED_TAG} ({str(released)[:7]}). "
        f"HEAD is therefore neither the tagged release commit nor the "
        f"release-preparation commit {RELEASE_PREP_COMMIT[:7]} "
        f"('{prep_subject}'), "
        f"which is {reachable} HEAD and is what the deposited Zenodo record was "
        "cut from."
    )
    if dirty:
        note += (
            f" WARNING: the build was allowed with --allow-dirty while {len(dirty)} path(s) "
            "carried uncommitted edits, so this commit hash does NOT describe every byte of "
            "this archive; see git_dirty_paths."
        )
    return {
        "git_available": True,
        "git_unavailable": False,
        "commit": head,
        "commit_format": "sha1, 40 hexadecimal characters",
        "branch": branch,
        "describe": git("describe", "--tags", "--always", "HEAD") or "no tag reachable",
        "released_tag": RELEASED_TAG,
        "released_tag_commit": released,
        "commits_after_released_tag": after,
        "release_preparation_commit": RELEASE_PREP_COMMIT,
        "git_dirty": bool(dirty),
        "git_dirty_count": len(dirty),
        "git_dirty_paths": dirty,
        "note": note,
    }


# --------------------------------------------------------------------------
# collection and scrubbing
# --------------------------------------------------------------------------


def collect(root: Path, out_dir: Path) -> tuple[list[Path], list[tuple[str, str]]]:
    """Return the files to archive and the (path, reason) pairs left behind."""
    excluded: list[tuple[str, str]] = []
    chosen: list[Path] = []

    def keep(path: Path) -> bool:
        rel = path.relative_to(root)
        parts = rel.parts
        if any(p in PRUNE_DIRS for p in parts[:-1]):
            return False
        if any(p.endswith(s) for p in parts for s in (".egg-info",)):
            return False
        if path.name in PRUNE_NAMES or path.suffix in PRUNE_SUFFIXES:
            return False
        try:
            resolved = path.resolve()
            if resolved == out_dir.resolve() or out_dir.resolve() in resolved.parents:
                return False
        except OSError:
            pass
        return True

    for name in INCLUDE_FILES:
        path = root / name
        if path.is_file():
            chosen.append(path)
        else:
            excluded.append((name, "listed for inclusion but absent from the working tree"))

    for name in INCLUDE_DIRS:
        directory = root / name
        if not directory.is_dir():
            excluded.append((name, "listed for inclusion but absent from the working tree"))
            continue
        for path in sorted(directory.rglob("*")):
            if path.is_file():
                if keep(path):
                    chosen.append(path)
                else:
                    excluded.append(
                        (str(path.relative_to(root)), "build/VCS/editor artefact or archive output")
                    )

    for entry in sorted(os.listdir(root)):
        if entry in INCLUDE_FILES or entry in INCLUDE_DIRS:
            continue
        reason = EXCLUDED_TOPLEVEL.get(entry)
        if reason is None:
            reason = "not on the reviewer include list"
        excluded.append((entry, reason))

    return sorted(set(chosen)), excluded


def scrub(text: str, rules: list[Rule]) -> tuple[str, list[str], int]:
    """Apply the whole token map; report which rules fired and how often."""
    used: list[str] = []
    total = 0
    for rule in rules:
        text, n = rule.apply(text)
        if n:
            used.append(rule.name)
            total += n
    return text, used, total


def read_entry(
    path: Path, rel: str, rules: list[Rule], *, license_mode: str
) -> tuple[bytes, list[dict[str, object]]]:
    """Bytes to archive plus one scrub record per file actually changed.

    ``--license preserved`` ships ``LICENSE`` byte-for-byte.  Substituting even
    only the holder's name would corrupt the notice - the bare given-name rule
    turns it into ``AuthorOne Zeng``, which names nobody and licenses nothing -
    and a half-scrubbed licence is worse than either complete choice.
    """
    raw = path.read_bytes()
    records: list[dict[str, object]] = []
    if license_mode == "preserved" and rel == "LICENSE":
        return raw, records
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw, records  # binary: PNG/CSV blobs are scrubbed by exclusion, not in place
    scrubbed, used, n = scrub(text, rules)
    if used:
        records.append({"path": rel, "rules": used, "occurrences": n})
    return scrubbed.encode("utf-8"), records


# --------------------------------------------------------------------------
# reviewer-facing notes written into the archive
# --------------------------------------------------------------------------

#: Replaced by the package version in each generated notice.
VERSION_MARKER = "VERSION_MARKER"

NOTICE_EN = """\
# PhyloSlicer VERSION_MARKER - note for peer reviewers

This copy of PhyloSlicer has been anonymised for double-blind review.  Every
author name, e-mail address, ORCID iD, institutional affiliation, funding
identifier, public repository URL and released DOI in the metadata has been
replaced by a neutral token such as `Author One Withheld`,
`AUTHOR_WITHHELD_FOR_PEER_REVIEW` or `Institution Withheld For Peer Review`.
The substitution is uniform across the whole tree, documents, code and tests
included, which is why the package's own metadata-consistency test still passes
here: the literals the test reads and the literals the metadata carries are the
same tokens.

Nothing else about the software has been changed.  The code, the reference data
and the documentation are the submitted versions.

`LICENSE` is carried in full, including all three BSD-3-Clause conditions and
the disclaimer; only the name of the copyright holder is withheld for review.
The licence grant itself is unchanged, and the working copy of the package -
which the editors can see on acceptance - names the holder.

## What you can run straight from this directory

```bash
python -m pytest -q                          # the full test suite
python examples/quickstart.py --sites 24 --species 60 --trees 3
python benchmarks/recompute_ratios.py        # re-derives the quoted speed-ups
```

`pyproject.toml` puts `src/` on the path for pytest, and the example and the
benchmark drivers insert it themselves, so those need no install step.  For the
cross-validation comparator and the CLI, use either form:

```bash
python validation/compare.py  # re-runs the cross-validation; no env setup needed
python -m pip install -e .                   # or install once, then run plainly
```

`REVIEW_ARCHIVE_PROVENANCE.json` next to this file records the package version,
the source commit this archive was built from, which files were scrubbed (by
rule, never by quoting the removed literal) and a SHA-256 for every entry.
"""

NOTICE_CN = """\
# PhyloSlicer VERSION_MARKER —— 致同行评审专家的说明

本副本已为双盲评审做匿名化处理。元数据中的作者姓名、邮箱、ORCID iD、所属机构、
资助编号、公开仓库地址与已发布的 DOI 均已替换为中性占位符，例如
`Author One Withheld`、`AUTHOR_WITHHELD_FOR_PEER_REVIEW`、
`Institution Withheld For Peer Review`。替换是全局统一的：文档、代码与测试使用
同一套占位符，因此软件包自身那条“元数据必须一致”的测试在这里依然成立——测试读取
的字面量与元数据里的字面量是同一个占位符。

软件本身没有做任何其他改动。代码、参考数据与文档都是送审版本。

`LICENSE` 全文保留，BSD-3-Clause 的三项条件与免责声明一字未删；仅为评审隐去了版权
持有人的姓名，授权条款本身不变，且工作目录中的原始文件仍然写明持有人。

## 可直接在本目录运行

```bash
python -m pytest -q                          # 完整测试套件
python examples/quickstart.py --sites 24 --species 60 --trees 3
python benchmarks/recompute_ratios.py        # 重算文中引用的加速比
```

`pyproject.toml` 已为 pytest 把 `src/` 加入路径，示例与基准脚本也自行插入该路径，
因此以上命令无需安装步骤。交叉验证比较器与命令行工具则可用以下任一形式：

```bash
python validation/compare.py  # 重跑交叉验证；无需设环境变量
python -m pip install -e .                   # 或先安装，再直接运行
```

同目录下的 `REVIEW_ARCHIVE_PROVENANCE.json` 记录了包版本号、构建所用的源码提交、
被脱敏的文件清单（按规则名记录，不复述被删除的字面量）以及每个条目的 SHA-256。
"""


def generated_entries(version: str) -> dict[str, bytes]:
    """Files the packager writes into the archive but keeps out of the tree.

    The version is stamped into each notice so a reviewer can tie the note to the
    copy in hand; ``REVIEW_ARCHIVE_PROVENANCE.json`` carries the full record.
    """
    return {
        "PEER_REVIEW_NOTICE.md": NOTICE_EN.replace(VERSION_MARKER, version).encode("utf-8"),
        "PEER_REVIEW_NOTICE_CN.md": NOTICE_CN.replace(VERSION_MARKER, version).encode("utf-8"),
    }


# --------------------------------------------------------------------------
# the guard
# --------------------------------------------------------------------------


def scan_blob(name: str, data: bytes, patterns: list[tuple[str, str, str, bool]]) -> list[Hit]:
    """Every pattern that survives in one archive entry (name included)."""
    hits: list[Hit] = []
    haystack = name.encode("utf-8")
    suffix = Path(name).suffix.lower()
    doc_like = suffix in DOC_SUFFIXES or suffix == ""
    for pname, desc, literal, doc_only in patterns:
        if doc_only and not doc_like:
            continue
        needle = literal.encode("utf-8")
        for where, blob in ((f"{name} (entry name)", haystack), (name, data)):
            start = blob.find(needle)
            if start < 0:
                continue
            window = blob[max(0, start - 60) : start + len(needle) + 60]
            context = window.decode("utf-8", "replace").replace("\n", "\\n")
            hits.append(Hit(where, pname, desc, context))
    return hits


#: An ORCID iD as it is written on the register: four dashed groups, the last one
#: ending in a digit or the checksum character ``X``.
ORCID_SHAPE = re.compile(r"\d{4}-\d{4}-\d{4}-\d{3}[0-9X]")

#: CJK ideographs (the block the archived Chinese documents are written in).
CJK_SHAPE = re.compile(r"[㐀-䶿一-鿿]")


def orcid_check_digit(first15: str) -> str:
    """The ISO 7064 Mod 11,2 check character for the leading 15 digits of an iD."""
    total = 0
    for ch in first15:
        total = (total + int(ch)) * 2
    remainder = (12 - (total % 11)) % 11
    return "X" if remainder == 10 else str(remainder)


def orcid_violations(name: str, data: bytes) -> list[Hit]:
    """CM-20: no malformed ORCID iD may ship, and no anonymisation token may look like one.

    The policy this tool enforces is *token, not identifier*: an iD field that has been
    withheld carries ``WITHHELD_FOR_PEER_REVIEW``.  A well-formed-looking string that
    resolves nowhere is worse than an obvious token, because a metadata harvester reads
    it as an error in the deposit - so every iD-shaped run in an archived text entry is
    re-checked against the checksum rule, and anything that fails is a violation.
    """
    if not data:
        return []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    hits: list[Hit] = []
    for found in set(ORCID_SHAPE.findall(text)) | set(ORCID_SHAPE.findall(name)):
        body = found.replace("-", "")
        want = orcid_check_digit(body[:15])
        if want == found[-1]:
            continue  # a genuine, checksum-valid iD: not this tool's business
        start = text.find(found)
        hits.append(
            Hit(
                name,
                "invalid_orcid",
                f"ORCID-shaped identifier {found} fails the ISO 7064 check digit (expected {want})",
                text[max(0, start - 60) : start + len(found) + 60].replace("\n", "\\n"),
            )
        )
    return hits


def marker_violations(name: str, data: bytes) -> list[Hit]:
    """RI-18: an unresolved authorial action must not survive into a shipped file."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    return [
        Hit(name, "unresolved_action", f"unresolved decision marker {marker!r}", marker)
        for marker in UNRESOLVED_ACTION_MARKERS
        if marker in text
    ]


#: RI-09/CM-02: a Zenodo DOI is the one identifier in the metadata that resolves
#: straight to the public deposit, and the deposit record names both authors, their
#: e-mails, their ORCID iDs and their institution.  The literal scrub rules catch the
#: two known digits; this shape catches *any* deposit id, including one added later,
#: which is why the guard is a regex and not another literal.  The token the scrub
#: writes (``zenodo.WITHHELD_FOR_PEER_REVIEW``) has no digits, so it never matches.
ZENODO_DOI_SHAPE = re.compile(r"10\.5281/zenodo\.([0-9]{6,})")


def doi_violations(name: str, data: bytes) -> list[Hit]:
    """CM-02: no resolvable Zenodo DOI may survive in the reviewer archive.

    ``10.5281/zenodo.<digits>`` is a click away from the authors' names, so an
    archived entry that still carries one de-anonymises the submission no matter how
    carefully the name literals were tokenised.  Only the digits are reported, so the
    message cannot itself widen the leak.
    """
    if not data:
        return []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    hits: list[Hit] = []
    for match in ZENODO_DOI_SHAPE.finditer(text):
        hits.append(
            Hit(
                name,
                "resolvable_zenodo_doi",
                "a Zenodo version/concept DOI with its record digits intact "
                f"({match.group(0)[:14]}…) survived: it resolves to the named deposit",
                "…/zenodo.[record id withheld here]",
            )
        )
    return hits


#: RI-17/CM-03: the shape of the bypass that defeated the guard once - a dict of
#: withheld literals, each spelled as fragments joined at runtime so neither the
#: scrubber nor the scanner could see it.  A file may not do that here.
#:
#: Exempt, on one test: a fragment that names nobody cannot deanonymise anything.
#: ``home_users_dir`` and its Windows spellings are generic path prefixes, and
#: ``tests/test_release_hygiene.py`` assembles them on purpose because it is the test
#: that sweeps the tree for personal paths - its own assertion would trip the sweep
#: otherwise.  ``synthetic_orcid`` is the fabricated all-zero iD an earlier build used
#: as its *token*; it is registered nowhere, and ``tests/test_review_archive.py`` has to
#: assemble it in order to plant it and prove the guard rejects it.  All four are still
#: rejected contiguously by ``scan_blob()`` - only the reassembly scan stands down.
FRAGMENT_EXEMPT_PATTERNS = frozenset(
    {
        "home_users_dir",
        "home_windows_dir",
        "windows_unc",
        "synthetic_orcid",
        "checkout_path",
    }
)

FRAGMENT_JOIN_FUNCS = frozenset({"_frag", "_join", "join", "concat"})


def _reassembled(node: ast.AST) -> str | None:
    """The string an expression evaluates to, if it is pure string-literal glue.

    Covers ``"a" + "b"``, ``"ab"`` constants, ``_frag("a", "b")``, ``"".join(("a",
    "b"))`` and any nesting of those - exactly the forms an anonymisation-bypassing
    table uses.  ``None`` for anything that is not a closed literal expression.
    """
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return node.value
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left, right = _reassembled(node.left), _reassembled(node.right)
        return left + right if left is not None and right is not None else None
    if isinstance(node, ast.Call):
        func = node.func
        args = node.args
        # "".join(<literal container>)
        if isinstance(func, ast.Attribute) and func.attr == "join":
            container = args[0] if args else None
            parts = _literal_sequence(container)
            sep = _reassembled(func.value) if isinstance(func.value, ast.Constant) else None
            if parts is not None and sep is not None:
                return sep.join(parts)
            return None
        # _frag("a", "b", ...) / any helper called by name
        if isinstance(func, ast.Name) and func.id in FRAGMENT_JOIN_FUNCS:
            parts = [_reassembled(a) for a in args]
            return "".join(parts) if all(p is not None for p in parts) else None
    return None


def _literal_sequence(node: ast.AST | None) -> list[str] | None:
    """Flatten a tuple/list/set/generator-of-literals into its string items."""
    if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
        out = [_reassembled(item) for item in node.elts]
        return out if all(o is not None for o in out) else None
    return None


def fragment_violations(name: str, data: bytes, literals: dict[str, str]) -> list[Hit]:
    """RI-17/CM-03: reject an archived file that reassembles a withheld literal.

    The guard compares every closed string-literal expression in an archived Python
    file against the scrub set, so splitting ``Zeng`` into ``"Ze" + "ng"`` - or into a
    helper call over fragments - is now a *failure* rather than an evasion.  Names of
    the offending literals are reported, never their content.
    """
    if Path(name).suffix != ".py" or not data:
        return []
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return []
    try:
        tree = ast.parse(text)
    except SyntaxError:
        return []
    hits: list[Hit] = []
    seen: set[str] = set()
    for node in ast.walk(tree):
        # only *reassembled* values: a plain contiguous literal is already caught by
        # scan_blob() (and would be scrubbed), so scanning constants here would just
        # double-report and punish files that legitimately spell an ordinary string
        if not isinstance(node, (ast.BinOp, ast.Call)):
            continue
        value = _reassembled(node)
        if not value or len(value) < 3:
            continue
        for pname, literal in literals.items():
            if literal and literal in value:
                key = f"{name}:{pname}:{node.lineno}"
                if key in seen:
                    continue
                seen.add(key)
                hits.append(
                    Hit(
                        name,
                        "reassembled_identity",
                        f"string expression at line {node.lineno} reassembles the withheld "
                        f"literal of rule '{pname}'",
                        f"{type(node).__name__} at line {node.lineno}",
                    )
                )
    return hits


def without_root(entries: dict[str, bytes], prefix: str) -> dict[str, bytes]:
    """Drop the single root directory the archive unpacks to, once and only once.

    ``_strip_prefix`` is written for archive member names; applying it twice - which is
    what happens when ``check()`` re-reads an archive that ``read_archive()`` has already
    normalised - would eat ``docs/`` off ``docs/USER_GUIDE_CN.md``.  The language census
    is keyed on paths inside the unpacked root in both directions, so it goes through
    here.
    """
    head = prefix + "/"
    return {(k[len(head) :] if k.startswith(head) else k): v for k, v in entries.items()}


def language_mix(entries: dict[str, bytes]) -> dict[str, object]:
    """CM-19: record which archived documents are Chinese, and how many.

    The note promises a reviewer "an English user guide, a Chinese user manual"; the
    promise needs a count that is read off the bytes actually shipped rather than a
    sentence in a readme, and this is that count.  It doubles as the identity surface
    census: every file here has been through the same scrub and the same guard.
    """
    per_file: dict[str, int] = {}
    for name, data in sorted(entries.items()):
        if Path(name).suffix.lower() not in DOC_SUFFIXES:
            continue
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            continue
        n = len(CJK_SHAPE.findall(text))
        if n:
            per_file[name] = n
    markdown = {k: v for k, v in per_file.items() if k.endswith(".md")}
    return {
        "detection": (
            "characters in the CJK Unified Ideographs blocks (U+3400-U+4DBF, U+4E00-U+9FFF) "
            "in the decoded text of every archived entry whose suffix is a document, code or "
            "configuration suffix"
        ),
        "documents_with_chinese": len(per_file),
        "markdown_documents_with_chinese": len(markdown),
        "chinese_counterpart_documents": sum(1 for k in per_file if "_CN." in k),
        "total_chinese_characters": sum(per_file.values()),
        "largest_documents": [
            {"path": k, "chinese_characters": v}
            for k, v in sorted(per_file.items(), key=lambda kv: (-kv[1], kv[0]))[:12]
        ],
        "note": (
            "Every English document of record ships beside a Chinese translation of it; the "
            "English text is authoritative and both were scrubbed and guard-scanned by the "
            "same pass, so the mirrored set adds no identity surface that the English set "
            "does not already have."
        ),
    }


def artefact_violations(names: list[str]) -> list[str]:
    """Entry names that must never appear in a review archive."""
    banned_dirs = {
        ".git",
        ".github",
        "__pycache__",
        ".pytest_cache",
        ".ruff_cache",
        ".mypy_cache",
        ".tox",
        ".venv",
        "venv",
        "build",
        "dist",
        "tools",
    }
    bad = []
    for name in names:
        parts = name.split("/")
        if any(p in banned_dirs for p in parts):
            bad.append(name)
        elif any(p.endswith(".egg-info") for p in parts[:-1]):
            bad.append(name)
        elif Path(name).name in PRUNE_NAMES or Path(name).suffix in {".pyc", ".pyo"}:
            bad.append(name)
    return sorted(set(bad))


def report_hits(hits: list[Hit], *, label: str) -> None:
    print(f"\nIDENTITY GUARD FAILED ({label}):", file=sys.stderr)
    for hit in hits:
        print(
            f"  {hit.where}\n    pattern: {hit.pattern} ({hit.description})\n"
            f"    context: …{hit.context}…",
            file=sys.stderr,
        )
    print(
        f"\n{len(hits)} surviving identifying match(es); the archive is not fit to send.",
        file=sys.stderr,
    )


# --------------------------------------------------------------------------
# writers
# --------------------------------------------------------------------------


def write_zip(path: Path, entries: list[tuple[str, bytes]], epoch: int) -> None:
    stamp = time.gmtime(epoch)[:6]
    with ZipFile(path, "w", ZIP_DEFLATED) as zf:
        for name, data in entries:
            info = ZipInfo(name, date_time=stamp)
            info.compress_type = ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o644 << 16
            zf.writestr(info, data)


def write_targz(path: Path, entries: list[tuple[str, bytes]], epoch: int) -> None:
    with open(path, "wb") as fh:
        with gzip.GzipFile(filename="", mode="wb", compresslevel=9, mtime=0, fileobj=fh) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.GNU_FORMAT) as tf:
                for name, data in entries:
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    info.mtime = epoch
                    info.mode = 0o644
                    info.type = tarfile.REGTYPE
                    info.uid = info.gid = 0
                    info.uname = info.gname = ""
                    tf.addfile(info, io.BytesIO(data))


def read_archive(path: Path) -> dict[str, bytes]:
    """Every entry of an existing archive, keyed by name without the root prefix."""
    out: dict[str, bytes] = {}
    if zipfile_name(path):
        with ZipFile(path) as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                out[_strip_prefix(info.filename)] = zf.read(info)
    else:
        with tarfile.open(path, "r:gz") as tf:
            for member in tf.getmembers():
                if not member.isfile():
                    continue
                handle = tf.extractfile(member)
                assert handle is not None
                out[_strip_prefix(member.name)] = handle.read()
    return out


def zipfile_name(path: Path) -> bool:
    return path.suffix == ".zip"


def _strip_prefix(name: str) -> str:
    parts = name.split("/")
    return "/".join(parts[1:]) if len(parts) > 1 and parts[0] else name


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# --------------------------------------------------------------------------
# build / check
# --------------------------------------------------------------------------


@dataclass
class Build:
    prefix: str
    entries: dict[str, bytes] = field(default_factory=dict)
    scrubbed: list[dict[str, object]] = field(default_factory=list)
    excluded: list[tuple[str, str]] = field(default_factory=list)
    notes: dict[str, str] = field(default_factory=dict)


def build(
    root: Path,
    out_dir: Path,
    *,
    license_mode: str,
    verbose: bool,
    state: dict[str, object] | None = None,
) -> Build:
    rules = rules_for(root)
    state = state if state is not None else worktree_state()
    version = package_version(root)
    archive = Build(prefix=f"PhyloSlicer_review_{version}")
    chosen, excluded = collect(root, out_dir)
    archive.excluded = excluded

    for path in chosen:
        rel = str(path.relative_to(root)).replace(os.sep, "/")
        data, records = read_entry(path, rel, rules, license_mode=license_mode)
        archive.entries[f"{archive.prefix}/{rel}"] = data
        archive.scrubbed.extend(records)

    for name, data in generated_entries(version).items():
        archive.entries[f"{archive.prefix}/{name}"] = data
        archive.notes[name] = "written by the packager; not present in the working tree"

    digest = {
        name[len(archive.prefix) + 1 :]: {
            "bytes": len(data),
            "sha256": sha256(data),
        }
        for name, data in sorted(archive.entries.items())
    }
    rules_used: dict[str, dict[str, object]] = {}
    for record in archive.scrubbed:
        for rule_name in record["rules"]:  # type: ignore[union-attr]
            entry = rules_used.setdefault(
                rule_name,
                {"description": _rule_description(rule_name), "files": [], "occurrences": 0},
            )
            entry["files"].append(record["path"])
            entry["occurrences"] = int(entry["occurrences"]) + 1
    provenance = {
        "tool": TOOL_VERSION,
        "tool_path": "tools/make_review_archive.py",
        "package": {"name": "phyloslicer", "version": version},
        "source_control": source_control(state),
        "build": {
            "python": sys.version.split()[0],
            "python_full": sys.version,
            "platform": sys.platform,
            "timestamp_utc": dt.datetime.fromtimestamp(
                fixed_epoch()[0], dt.timezone.utc
            ).isoformat(),
            "timestamp_source": fixed_epoch()[1],
            "deterministic": True,
        },
        "anonymisation": {
            "policy": (
                "Uniform literal-to-token substitution across every archived text file "
                "(documents, code, metadata and tests alike), plus exclusion of build, "
                "VCS and repository-automation artefacts. Source files in the working "
                "tree are never modified."
            ),
            "licence_mode": license_mode,
            #: RI-18: the licence-mode question, closed.  A reviewer (or an editor)
            #: reading the sidecar gets the decision and the statement of fact, never an
            #: open action, and `check()` refuses an archive whose record says otherwise.
            "licence_decision": {
                "resolved": True,
                "decided_mode": RESOLVED_LICENSE_MODE,
                "this_build": license_mode,
                "licence": "BSD-3-Clause",
                "statement": RESOLVED_LICENSE_STATEMENT,
                "deviation": (
                    None
                    if license_mode == RESOLVED_LICENSE_MODE
                    else (
                        f"this build used --license {license_mode}, which is not the recorded "
                        f"decision ({RESOLVED_LICENSE_MODE}); it names the copyright holders and "
                        "is a deliberate, sidecar-declared exception, not an unresolved question"
                    )
                ),
            },
            "orcid_policy": (
                "withheld as a non-identifier token ('"
                f"{WITHHELD_ORCID_TOKEN}'), never a synthetic iD; every ORCID-shaped string in "
                "an archived entry is re-checked against the ISO 7064 Mod 11,2 rule and "
                "rejected if it fails"
            ),
            "rules": [
                {"rule": name, **spec}  # never the removed literal, only its description
                for name, spec in sorted(rules_used.items())
            ],
            "scrubbed_files": archive.scrubbed,
            "declared_exceptions": (
                [
                    {
                        "path": "LICENSE",
                        "kind": "copyright holder name retained, notice shipped verbatim",
                        "reason": (
                            "BSD-3-Clause is granted by a named holder; the notice cannot be "
                            "stripped without altering the licence. The archive is supplied "
                            "privately to the editorial office only."
                        ),
                        "consequence": (
                            "tests/test_release_hygiene.py is scrubbed like every other text "
                            "file, so its assertion that the licence names the holder no longer "
                            "matches this verbatim LICENSE: in 'preserved' mode the archived "
                            "suite fails that one test. Use 'placeholder' for an archive that "
                            "passes its own suite."
                        ),
                    }
                ]
                if license_mode == "preserved"
                else []
            ),
            "third_party_names_kept": [
                "NOTICE.md / NOTICE_CN.md and the CITATION.cff reference entry name the "
                "treesliceR authors, who are third parties, not the authors of this submission"
            ],
        },
        "excluded": [{"path": p, "reason": r} for p, r in archive.excluded],
        "language_mix": language_mix(without_root(archive.entries, archive.prefix)),
        "generated_inside_archive": archive.notes,
        "files": [
            {"path": name, "bytes": info["bytes"], "sha256": info["sha256"]}
            for name, info in sorted(digest.items())
        ],
        "summary": {
            "file_count": len(archive.entries) + 1,
            "total_bytes": sum(len(d) for d in archive.entries.values()),
            "per_directory": per_directory(archive.entries),
        },
        "note_for_reviewers": (
            "See PEER_REVIEW_NOTICE.md / PEER_REVIEW_NOTICE_CN.md. Metadata fields reading "
            "'Withheld' are anonymisation tokens, not placeholders for unfinished work: the "
            "real attribution is carried by the manuscript title page and by the public "
            "release, and is restored verbatim on acceptance."
        ),
    }
    sidecar = json.dumps(provenance, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    archive.entries[f"{archive.prefix}/REVIEW_ARCHIVE_PROVENANCE.json"] = sidecar.encode("utf-8")
    if verbose:
        print(f"provenance sidecar: {len(sidecar)} bytes over {len(digest)} hashed entries")
    return archive


def _rule_description(name: str) -> str:
    for rule in STATIC_RULES:
        if rule.name == name:
            return rule.description
    return {
        "checkout_path": "absolute path of this working tree",
        "home_directory": "absolute home directory of the builder",
        "local_account": "local account name",
    }.get(name, name)


def per_directory(entries: dict[str, bytes]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for name, data in entries.items():
        rel = _strip_prefix(name)
        top = rel.split("/")[0] if "/" in rel else "(top level)"
        entry = counts.setdefault(top, {"files": 0, "bytes": 0})
        entry["files"] += 1
        entry["bytes"] += len(data)
    return {k: v for k, v in sorted(counts.items())}


def finalise(out_dir: Path, fmt: str, archive: Build, verbose: bool) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{archive.prefix}.{fmt}"
    epoch, source = fixed_epoch()
    ordered = sorted(archive.entries.items())
    if fmt == "zip":
        write_zip(path, ordered, epoch)
    else:
        write_targz(path, ordered, epoch)
    if verbose:
        print(f"timestamp source: {source} -> {epoch}")
    return path


def summary(
    path: Path,
    archive: Build,
    *,
    verbose: bool,
    state: dict[str, object] | None = None,
) -> None:
    data = path.read_bytes()
    entries = len(archive.entries)
    total = sum(len(v) for v in archive.entries.values())
    print(f"\narchive            : {path}")
    print(f"entries            : {entries}")
    print(f"uncompressed bytes : {total:,}")
    print(f"archive bytes      : {len(data):,}  (SHA-256 {sha256(data)})")
    print("per directory      :")
    for name, info in per_directory(archive.entries).items():
        print(f"  {name:<22} {info['files']:>4} files  {info['bytes']:>12,} bytes")
    mix = language_mix(without_root(archive.entries, archive.prefix))
    print(
        f"language mix       : {mix['documents_with_chinese']} archived documents carry "
        f"Chinese ({mix['markdown_documents_with_chinese']} of them .md, "
        f"{mix['chinese_counterpart_documents']} named *_CN.*), "
        f"{mix['total_chinese_characters']:,} Chinese characters in total"
    )
    control = source_control(state) if state is not None else None
    if control is not None:
        commit = control["commit"]
        if commit is None:
            provenance_line = "commit null (git_unavailable: no claim made)"
        elif control["git_dirty"]:
            provenance_line = (
                f"commit {str(commit)[:7]} + {control['git_dirty_count']} uncommitted "
                "path(s) recorded in the sidecar (--allow-dirty)"
            )
        else:
            provenance_line = f"commit {str(commit)[:7]} (worktree clean at build time)"
        print(f"provenance         : {provenance_line}")
    if verbose:
        print(f"scrubbed files     : {len(archive.scrubbed)}")
        for record in archive.scrubbed:
            print(f"  {record['path']:<34} {', '.join(record['rules'])}")
        print(f"excluded entries   : {len(archive.excluded)}")
        for name, reason in archive.excluded:
            print(f"  {name:<34} {reason}")


def check(path: Path, *, license_mode: str, verbose: bool) -> list[Hit]:
    """Re-read an existing archive and return every identity/artefact violation."""
    rules = rules_for(ROOT)
    patterns = scan_patterns(rules)
    entries = read_archive(path)
    #: rule name -> removed literal, for the fragment-reassembly scan below
    literals = {
        pname: literal
        for pname, _desc, literal, _doc_only in patterns
        if pname not in FRAGMENT_EXEMPT_PATTERNS
    }
    hits: list[Hit] = []
    for name, data in sorted(entries.items()):
        hits.extend(scan_blob(name, data, patterns))
        hits.extend(orcid_violations(name, data))
        hits.extend(marker_violations(name, data))
        hits.extend(doi_violations(name, data))
        hits.extend(fragment_violations(name, data, literals))
    for bad in artefact_violations(list(entries)):
        hits.append(Hit(bad, "artefact", "build/VCS artefact in the archive", "-"))
    sidecar_name = "REVIEW_ARCHIVE_PROVENANCE.json"
    if sidecar_name not in entries:
        hits.append(Hit(sidecar_name, "provenance", "the provenance sidecar is missing", "-"))
    else:
        provenance = json.loads(entries[sidecar_name])
        declared = {item["path"]: item["sha256"] for item in provenance["files"]}
        for name, data in sorted(entries.items()):
            if name == sidecar_name:
                continue
            want = declared.get(name)
            if want is None:
                hits.append(Hit(name, "undeclared", "entry absent from the sidecar", "-"))
            elif want != sha256(data):
                hits.append(Hit(name, "hash_mismatch", "SHA-256 differs from the sidecar", "-"))
        for name in declared:
            if name not in entries:
                hits.append(Hit(name, "missing", "declared by the sidecar but absent", "-"))
        hits.extend(sidecar_violations(provenance, entries, license_mode))
        if license_mode != "preserved":
            for exc in provenance["anonymisation"]["declared_exceptions"]:
                hits.append(Hit(exc["path"], "declared_exception", exc["reason"], "-"))
        if verbose:
            control = provenance["source_control"]
            commit = control["commit"]
            decided = provenance["anonymisation"]["licence_decision"]["decided_mode"]
            print(
                f"sidecar: {len(declared)} hashed entries, "
                f"{len(provenance['anonymisation']['scrubbed_files'])} scrubbed files, "
                f"version {provenance['package']['version']}, "
                f"commit {commit[:7] if commit else 'null (git_unavailable)'}, "
                f"licence mode {provenance['anonymisation']['licence_mode']} "
                f"(recorded decision {decided}), "
                f"{provenance['language_mix']['documents_with_chinese']} documents carry Chinese"
            )
    return [h for h in hits if not _allowed_by_mode(h, license_mode)]


def sidecar_violations(
    provenance: dict[str, object], entries: dict[str, bytes], license_mode: str
) -> list[Hit]:
    """NV-10 / NV-04 / RI-18 / CM-19: the self-claims the sidecar makes must hold.

    Every one of these is a claim the archive can check itself: a commit that is named
    has to be a real 40-character SHA-1 (or be absent, and say so), the licence mode has
    to be the recorded decision or an acknowledged deviation from it, and the language
    census has to match the bytes next to it.
    """
    hits: list[Hit] = []
    control = provenance["source_control"]
    commit = control["commit"]
    if control.get("git_unavailable"):
        if commit is not None:
            hits.append(
                Hit(
                    "REVIEW_ARCHIVE_PROVENANCE.json",
                    "commit_claim",
                    "git_unavailable is set but a commit string was written",
                    str(commit),
                )
            )
    elif commit is None:
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "commit_claim",
                "a commit is required when git metadata is readable",
                "null",
            )
        )
    elif not re.fullmatch(r"[0-9a-f]{40}", str(commit)):
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "commit_claim",
                "the commit is not a 40-character lower-case SHA-1",
                str(commit),
            )
        )
    if control.get("git_available") and control.get("git_dirty"):
        paths = control.get("git_dirty_paths") or []
        if len(paths) != control.get("git_dirty_count"):
            hits.append(
                Hit(
                    "REVIEW_ARCHIVE_PROVENANCE.json",
                    "dirty_claim",
                    "git_dirty_count disagrees with the length of git_dirty_paths",
                    f"{control.get('git_dirty_count')} vs {len(paths)}",  # type: ignore[arg-type]
                )
            )
        if not paths or "uncommitted" not in str(control.get("note", "")):
            hits.append(
                Hit(
                    "REVIEW_ARCHIVE_PROVENANCE.json",
                    "dirty_claim",
                    "a --allow-dirty build must list the uncommitted paths and disclaim HEAD",
                    "-",
                )
            )
    anonymisation = provenance["anonymisation"]
    decision = anonymisation.get("licence_decision") or {}
    if decision.get("resolved") is not True:
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "unresolved_licence_mode",
                "the licence mode is recorded without a resolved decision",
                "-",
            )
        )
    if decision.get("decided_mode") != RESOLVED_LICENSE_MODE:
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "unresolved_licence_mode",
                f"the recorded decision must be the resolved mode {RESOLVED_LICENSE_MODE!r}",
                str(decision.get("decided_mode")),
            )
        )
    if decision.get("this_build") != anonymisation.get("licence_mode"):
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "licence_mode_mismatch",
                "licence_decision.this_build disagrees with anonymisation.licence_mode",
                "-",
            )
        )
    if anonymisation.get("licence_mode") != license_mode:
        hits.append(
            Hit(
                "REVIEW_ARCHIVE_PROVENANCE.json",
                "licence_mode_mismatch",
                "the sidecar records a different licence mode from the one being checked",
                f"{anonymisation.get('licence_mode')} vs {license_mode}",
            )
        )
    recomputed = language_mix(entries)
    declared_mix = provenance.get("language_mix") or {}
    for key in ("documents_with_chinese", "total_chinese_characters"):
        if declared_mix.get(key) != recomputed[key]:  # type: ignore[union-attr]
            hits.append(
                Hit(
                    "REVIEW_ARCHIVE_PROVENANCE.json",
                    "language_mix",
                    f"{key} does not match the shipped bytes",
                    f"{declared_mix.get(key)} vs {recomputed[key]}",  # type: ignore[union-attr]
                )
            )
    return hits


def _allowed_by_mode(hit: Hit, license_mode: str) -> bool:
    """In 'preserved' mode the copyright line in LICENSE is a declared exception."""
    return license_mode == "preserved" and hit.where == "LICENSE"


LICENSE_WARNING = (
    "WARNING: --license preserved ships LICENSE verbatim, so the archive does name the "
    "copyright holders (Zichao Zeng and Yinzhao Wang, per CM-21). The guard accepts exactly "
    "that one hit, in exactly that one file, and records it as a declared exception; "
    "tests/test_release_hygiene.py, which is scrubbed like every other text file, then fails "
    "its licence-attribution assertion inside the archive. Say so in the cover letter."
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--out",
        default="dist/review_archive",
        help="output directory (default: dist/review_archive, already git-ignored)",
    )
    ap.add_argument("--format", choices=("zip", "tar.gz"), default="zip")
    ap.add_argument(
        "--check-only", metavar="ARCHIVE", help="verify an existing archive instead of building one"
    )
    ap.add_argument(
        "--license",
        choices=("placeholder", "preserved"),
        default="placeholder",
        help="'placeholder' (default) withholds the copyright holder's name; "
        "'preserved' ships LICENSE verbatim as a declared exception",
    )
    ap.add_argument(
        "--allow-dirty",
        action="store_true",
        help="build even though the working tree carries uncommitted edits; the sidecar then "
        "lists them in source_control.git_dirty_paths and disclaims HEAD as a full description "
        f"(default: refuse, {DIRTY_MESSAGE.splitlines()[0]})",
    )
    ap.add_argument(
        "--allow-no-git",
        action="store_true",
        help="build even though no git metadata is readable; the sidecar then records "
        "commit = null with git_unavailable = true, i.e. no commit claim at all "
        "(default: refuse rather than write a placeholder hash)",
    )
    ap.add_argument("--verbose", action="store_true", help="list scrubbed and excluded files")
    ap.add_argument("--root", default=str(ROOT), help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    out_dir = (root / args.out).resolve() if not Path(args.out).is_absolute() else Path(args.out)
    if args.license == "preserved":
        print(LICENSE_WARNING, file=sys.stderr)
    if args.license != RESOLVED_LICENSE_MODE:
        print(
            f"note: the recorded licence decision is {RESOLVED_LICENSE_MODE!r}; this build "
            f"deviates from it and the sidecar labels the deviation.",
            file=sys.stderr,
        )

    if args.check_only:
        path = Path(args.check_only)
        if not path.is_absolute():
            path = root / path
        if not path.is_file():
            print(f"no such archive: {path}", file=sys.stderr)
            return 2
        hits = check(path, license_mode=args.license, verbose=args.verbose)
        if hits:
            report_hits(hits, label=path.name)
            return 1
        print(
            f"identity guard PASSED: {path} ({len(path.read_bytes()):,} bytes, "
            f"{len(read_archive(path))} entries, licence mode {args.license})"
        )
        return 0

    state = worktree_state()
    enforce_provenance_policy(state, allow_dirty=args.allow_dirty, allow_no_git=args.allow_no_git)
    archive = build(root, out_dir, license_mode=args.license, verbose=args.verbose, state=state)
    path = finalise(out_dir, args.format, archive, args.verbose)

    # The guard runs on the bytes that were actually written, not on the
    # in-memory model: what a reviewer receives is what is checked.
    hits = check(path, license_mode=args.license, verbose=args.verbose)
    if hits:
        report_hits(hits, label=path.name)
        return 1
    summary(path, archive, verbose=args.verbose, state=state)
    print(
        f"\nidentity guard PASSED: no author identity, local path or cache artefact survives, "
        f"no unresolved action marker survives, and no malformed ORCID iD survives "
        f"(licence mode: {args.license}, the recorded decision being "
        f"{RESOLVED_LICENSE_MODE} = BSD-3-Clause with the holder names tokenised)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
