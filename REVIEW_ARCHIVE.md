# The anonymised peer-review archive

Chinese version: [REVIEW_ARCHIVE_CN.md](REVIEW_ARCHIVE_CN.md).

Ecography's Software Notes rules let a submission be reviewed from a
supplementary archive instead of an anonymous repository mirror: *"It can be
made available to reviewers either via a service such as Anonymous Github … or
by submitting the full source-code package and data as a supplementary zip
file."*  This document describes the packager that produces that archive, what
it removes, what it deliberately keeps, and how to verify the result.

The claim this backs up in the submission checklist is narrow and testable:
**CI builds the review archive through the packager, and the job fails if any
personal path or author identity survives.**  That claim is enforced by
`tools/make_review_archive.py` (which re-reads and scans its own output) and by
`tests/test_review_archive.py` (which unpacks the archive and runs it), and it
runs as a named step in `.github/workflows/tests.yml`.

## Build it

```bash
python tools/make_review_archive.py                      # -> dist/review_archive/PhyloSlicer_review_<version>.zip
python tools/make_review_archive.py --format tar.gz
python tools/make_review_archive.py --out /tmp/review --verbose
python tools/make_review_archive.py --allow-dirty        # only if the sidecar may advertise uncommitted paths
python tools/make_review_archive.py --allow-no-git       # only where no .git exists; records commit = null
```

Stdlib only, no install step and no `PYTHONPATH`: the version is read out of
`src/phyloslicer/__init__.py` rather than imported, so the tool carries no
runtime dependency.  The output directory is git-ignored (`dist/review_archive/`
and `review_archive/` are both listed in `.gitignore`).

Exit code is `0` only when the guard passed; `1` when identity or an artefact
survives, and `1` with a `PROVENANCE GATE` message when the archive could not make
an honest provenance claim (no git metadata, or a dirty worktree) and the matching
`--allow-*` flag was not given; `2` when the named archive does not exist.

## Verify it

```bash
# re-check an archive you already have, without rebuilding anything
python tools/make_review_archive.py --check-only dist/review_archive/PhyloSlicer_review_0.1.0.zip

# the whole contract, including the negative cases
python -m pytest -q tests/test_review_archive.py
```

`--check-only` re-reads every entry, scans it for the identity patterns, checks
the SHA-256 of each entry against the provenance sidecar, and rejects any build
or VCS artefact.  Two of the tests plant an author name and a `__pycache__` file
into a copy of a good archive and assert that the guard *fails* — a check that
cannot fail is worth nothing.

Reproducibility: entries are sorted, every timestamp is the commit date of HEAD
(override with `SOURCE_DATE_EPOCH`), archive metadata is normalised (mode 0644,
uid/gid 0, no owner names), and gzip is written with `mtime=0`.  Two consecutive
builds are byte-identical; `test_two_builds_are_byte_identical` asserts it, and
the recorded hashes are:

| archive | SHA-256 |
| --- | --- |
| `PhyloSlicer_review_0.1.0.zip` | `d3750aa2e2abd553cf19fc3343a9440d270b953d01dea7dfa8e22a61b81d1e79` |
| `PhyloSlicer_review_0.1.0.tar.gz` | `0af47241b70ed8069f895fc774870495d3806a2181a6fd1004f3b7e1198dd59a` |

(Built from commit `0cba60f` with the worktree clean, three times with CPython 3.11;
all three builds are byte-identical, and the archive holds 263 entries. Rebuild and
compare if you need to - the point of printing the digests is that they must match.  `SOURCE_DATE_EPOCH` is honoured and is the only knob that moves the
timestamps.)

Unpack it and it is a working checkout:

```bash
unzip PhyloSlicer_review_0.1.0.zip && cd PhyloSlicer_review_0.1.0
python -m pytest -q
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python -m pip install -e .
```

## What is in it, and what is kept out

Included: `pyproject.toml`, `MANIFEST.in`, `.gitignore`, `LICENSE`,
`CITATION.cff`, `.zenodo.json`, `src/`, `tests/`, `examples/`, `docs/`,
`benchmarks/`, `validation/` (reference CSVs, `manifest.txt` and the shipped
`equivalence_report.csv` included), `testing_suite/`, and every root document —
`README`, `CHANGELOG`, `CONTRIBUTING`, `RELEASE`, `CODE_OF_CONDUCT`, `SECURITY`,
`NOTICE`.

Two rules drive that list rather than taste:

* a reviewer must be able to install, run and re-verify the package;
* the package's own hygiene tests must pass inside the archive, so anything they
  read has to be there — `test_version_licence_and_release_date_are_stated`
  reads `CITATION.cff`, `.zenodo.json` and all ten manual/changelog files, and
  `test_gitignore_excludes_artifacts_that_carry_absolute_paths` reads
  `.gitignore` itself.

Kept out, with the reason recorded in the sidecar and printed by `--verbose`:
`.git/`, `.github/` (its `CODEOWNERS` declares `@ZengZichao` as owner of every
file), `.pre-commit-config.yaml`, `__pycache__/`, `.pytest_cache/`,
`.ruff_cache/`, `.mypy_cache/`, `*.egg-info/`, `build/`, `dist/`, `.DS_Store`,
`*.pyc`, the output directory itself, and — importantly — `tools/` and
`REVIEW_ARCHIVE.md`/`REVIEW_ARCHIVE_CN.md`.  The packager's scrub table and this
document necessarily contain the literals they remove; shipping them would undo
the anonymisation, and neither is needed to run the package.  The archive
instead carries `PEER_REVIEW_NOTICE.md` / `PEER_REVIEW_NOTICE_CN.md`, generated
by the packager, which tell the reviewer what the placeholder tokens mean.

Every archived document keeps its Chinese counterpart, because the project's own
hygiene rule requires it: the generated notice ships as a pair, and excluding
`.github/` removes a whole directory rather than half of a pair.

## What was scrubbed, and how

The audit's starting list was checked against the tree before anything was
written.  Hits found, with file and line, and the handling:

| where | what | handling |
| --- | --- | --- |
| `.github/CODEOWNERS:3`, `.github/ISSUE_TEMPLATE/config.yml:4` | owner handle, repo URL | whole `.github/` excluded |
| `LICENSE:3` | `Copyright (c) 2026, Zichao Zeng and Yinzhao Wang` | both holder names tokenised as a unit, licence text intact (see below) |
| `pyproject.toml:12,56-60` | author name + e-mail, five repository URLs | tokenised |
| `CITATION.cff:14-27,31,35` | both authors: names, e-mails, ORCIDs, affiliation, Zenodo DOI, repo URL | tokenised |
| `.zenodo.json:11-19,31,35` | same, as `creators`/`related_identifiers` | tokenised |
| `README.md:3-6,12,128,143-145` / `README_CN.md:3-6,12,117,132-134` | CI and Zenodo badges, repo URL, DOI, author block with affiliation and ORCID | tokenised |
| `SECURITY.md:19` / `SECURITY_CN.md:18` | reporter e-mail | tokenised |
| `CODE_OF_CONDUCT.md:65,137` / `CODE_OF_CONDUCT_CN.md:50` | reporter e-mail | tokenised |
| `CONTRIBUTING.md:13` / `CONTRIBUTING_CN.md:14`, `RELEASE.md:85` / `RELEASE_CN.md:80` | clone/install URL, Zenodo GitHub settings URL | owner segment tokenised |
| `tests/test_release_hygiene.py:191,194,202` | the hygiene test *asserts* the real name, e-mail and URL | tokenised by the same map (see below) |
| `docs/USER_GUIDE.md`, `docs/USER_GUIDE_CN.md`, `validation/*.py`, `examples/*.py`, `benchmarks/*`, `testing_suite/*` | checked: **no** author string, repo URL or absolute local path found | shipped unchanged |

Nothing matched a home-directory path (the macOS or Windows user-directory
prefixes, or this checkout's own absolute path), and nothing matched the grant
numbers `42422209`, `42272354` or `2023YFC3108600`: those identifiers live in the
manuscript, not in the repository.  They are still in the guard's reject list, so
if one is ever committed the build fails instead of leaking it.  (This document
names no absolute path precisely because `test_release_hygiene.py` scans it, too:
it is a scanned `.md` file in the source tree, and quoting the pattern it forbids
would fail the suite that the archive is meant to protect.)

Scrubbing is a **uniform literal-to-token substitution applied to every archived
text file at once** — documents, code, metadata *and tests*.  That is the
load-bearing choice: `tests/test_release_hygiene.py` requires `pyproject.toml`,
`CITATION.cff`, `.zenodo.json`, `README*.md` and `LICENSE` to name the same
author and repository, so scrubbing only some of them would hand the reviewer a
package whose own consistency test fails.  Substituting everywhere keeps the
literals the test reads and the literals the metadata carries identical, and the
suite passes unmodified inside the archive.  Tokens used:
`Author One Withheld`, `AuthorTwo`, `Withheld`, `AUTHOR_WITHHELD_FOR_PEER_REVIEW`,
`author.one.withheld@example.invalid` / `author.two.withheld@example.invalid`
(RFC 2606 `.invalid`, so no address can be contacted by accident),
`WITHHELD_FOR_PEER_REVIEW` in the ORCID fields, `Institution Withheld For Peer Review`,
`School Withheld For Peer Review`, `10.5281/zenodo.WITHHELD_FOR_PEER_REVIEW`,
`GRANT_WITHHELD_FOR_PEER_REVIEW`.  Two extra rules exist because of the
substitution's own blind spots: the home directory and this checkout's absolute
path (derived at build time from `Path.home()` and the account name, not
hard-coded), and the fact that `README.md:143` wraps `Shanghai⏎Jiao Tong
University` across a line break, which the full institution string cannot match —
without the `Jiao Tong University` / `Shanghai` rules the affiliation survived in
plain sight.

The source files are never touched: the scrubbed content is written into the
archive only, and `test_source_files_were_not_modified` pins that the working
tree still names the real authors.

**ORCID fields carry a token, never an identifier.**  An earlier build replaced the
two authors' iDs with the synthetic strings `0000-0000-0000-0001` and `-0002`, which
look like registered identifiers and resolve nowhere — a metadata harvester reads
them as a defect in the deposit rather than as anonymisation.  (The report called
them checksum-invalid; strictly they are not, because the all-zero stem happens to
make `1` the ISO 7064 check character, which is precisely why a checksum alone is
not a safe test.)  The scrub now writes `WITHHELD_FOR_PEER_REVIEW` into
`CITATION.cff` and `.zenodo.json`, keeping the key so the metadata stays parseable,
and the guard enforces the policy in two ways: the synthetic form is a scan-only
literal, and `orcid_violations()` re-checks **every** ORCID-shaped run in every
archived entry against the ISO 7064 Mod 11-2 rule and fails the build on any that
does not satisfy it.  A genuine third-party iD therefore still passes, and a planted
fake one cannot.

### Names that are kept on purpose

`NOTICE.md`, `NOTICE_CN.md` and the `references:` entry in `CITATION.cff` name
the treesliceR authors (Araujo, Ferreira, Nakamura, Coelho, Rangel) and their
ORCIDs.  They are third parties, not the authors of this submission; stripping
them would remove a required attribution and would not anonymise anything.  The
guard therefore rejects the two submission authors' ORCIDs by literal, not any
`orcid.org/…` string, and does not treat `Araujo` as identity.  This is recorded
in the sidecar under `third_party_names_kept`.

## `LICENSE`: the one requirement that cannot be satisfied twice

BSD-3-Clause is granted *by a named holder*; condition 1 requires the copyright
notice to be retained, and the project's own hygiene test asserts
`assert "Zichao Zeng" in license_text, "the licence must name the copyright
holder"`.  So either the archive names the author, or it does not ship a
licence.  There is no third option that is both honest and legal, and a
licence-stripped archive is not acceptable — that would misrepresent the terms
under which the code is offered.

**The decision is made and recorded: the submitted archive is the `placeholder`
build of the BSD-3-Clause licence — the notice ships in full and only the holders'
names are tokenised**, giving `Copyright (c) 2026, Author One Withheld and Author
Two Withheld`.  All three
conditions and the disclaimer are verbatim, so nothing about the grant changes
for the reviewer; the notice is restored to the real holders in the public
release, which is the copy that actually has to carry it.  The alternative mode is
implemented rather than merely described:

```bash
python tools/make_review_archive.py --license preserved
```

`preserved` ships `LICENSE` verbatim with the real names, records it in the
sidecar under `declared_exceptions` with the reason, and the guard then permits
that one hit in that one file and still fails on any other.  Choose `preserved`
if you would rather the licence stay untouched and rely on the archive going only
to the editorial office; choose the default if you want the reviewer-facing copy
to contain no name at all.  Note that `preserved` mode makes the archive
non-anonymous in one line, so the accompanying cover letter must say so, and that
it costs the archive its own suite: measured on a `preserved` build,
`tests/test_release_hygiene.py` inside the archive gives 10 passed / 1 failed,
the failure being `assert "Author One Withheld" in license_text` — the scrubbed
test compared against the verbatim licence.  That is the collision, stated
plainly, rather than a check weakened to hide it.

Nothing about this is left as an open question for the author: the mode is pinned
in the packager as `RESOLVED_LICENSE_MODE = "placeholder"`, every build writes it
into the sidecar as `anonymisation.licence_decision` (`resolved: true`, the
decided mode, the mode this build used, the BSD-3-Clause statement above, and a
`deviation` note when a build deliberately departs from it), the guard's closing
line prints the resolved value, and `check()` rejects any archive whose sidecar
leaves the decision unresolved or whose shipped files carry an unfinished-decision
marker (`AUTHOR ACTION`, `confirm which mode`, `TODO:`).

## Provenance sidecar

`REVIEW_ARCHIVE_PROVENANCE.json` is written into the archive and records: the
package version read from `src/phyloslicer/__init__.py`; the commit of the
working tree, its branch, `git describe`, the released tag `v0.1.0` and its
commit, and how HEAD relates to the release-preparation commit `b151683`
(**HEAD is a post-release maintenance commit, not the released tag** — the
deposited Zenodo record was cut from `b151683`, so an editor comparing the
archive with the DOI lands on a different commit and needs to be told); the
Python version and platform that built it; the single fixed timestamp and where
it came from; the licence mode **and the resolved licence decision with this build's
mode and any deviation from it**; the ORCID policy; the language census
(`language_mix`: how many archived documents carry Chinese characters, how many of
those are the `*_CN.md` counterparts, the character totals and the twelve largest —
the evidence behind the note's "English user guide, a Chinese user manual" claim);
the scrub rules that fired, by name, and the
per-file list of scrubbed paths with their occurrence counts; every excluded
path with its reason; and a SHA-256 plus byte count for every file in the
archive.

The sidecar deliberately never quotes a *removed* literal — only the rule name
("first author e-mail") and its replacement.  A provenance record that listed
`Zichao Zeng → Author One Withheld` would leak the very thing it documents, so
that mapping lives here and in the packager, both of which stay out of the
archive.  `test_provenance_never_quotes_a_removed_literal` enforces it.  The
sidecar hashes every entry except itself (a file cannot contain its own digest);
`--check-only` knows that one exception and rejects anything else that is
undeclared, missing, or whose bytes differ.

One consequence of reading the working tree at build time is that the archive can
only claim what the tree can prove, and the packager now refuses to overclaim.  Two
states stop the build:

* **no git metadata** (no `.git`, or no `git` on PATH).  This used to write the
  literal `"unknown"` into `source_control.commit`, which then failed the shipped
  suite as `assert 7 == 40` — a reviewer whose tooling drops `.git` met a red test
  and a hash-shaped lie.  The build is refused with a `PROVENANCE GATE` message and
  an explanation; `--allow-no-git` makes it proceed and write `commit: null` with
  `git_unavailable: true`, i.e. no commit claim at all, which
  `tests/test_review_archive.py` reads and skips the commit assertions on instead
  of failing.
* **uncommitted edits in the working tree.**  Naming HEAD while dozens of paths
  differ from it is a false provenance claim, so the build is refused unless
  `--allow-dirty` is passed; when it is, the sidecar records `git_dirty: true`,
  `git_dirty_count` and the full `git_dirty_paths` list, and `source_control.note`
  disclaims HEAD in as many words.  `check()` rejects a dirty build that did not
  say so, and rejects any commit that is not a 40-character SHA-1 when git *was*
  readable.

So the honest workflow is: let the orchestrator commit, build from the settled
checkout (no flag needed), and let the CI step — which runs after the suite — be
the one that produces the file you send.  A `--allow-dirty` archive is legitimate
only because it advertises itself as one.

## What a reviewer is expected to be able to run

Everything in the list below was executed inside an unpacked archive, from the
archive alone, with no install step:

| command | result (measured, CPython 3.11) |
| --- | --- |
| `python -m pytest -q` | the full suite passes; 12 s inside the test harness, 27 s on a cold first extraction, release-hygiene tests included |
| `python examples/quickstart.py --sites 24 --species 60 --trees 3` | runs in ~0.6 s and prints the fitted table |
| `python validation/compare.py` | `verdict: PASS` in ~2.5 s, with no `PYTHONPATH` and no install |
| `python benchmarks/recompute_ratios.py` | re-derives every quoted speed-up range |

Every one of those runs works as written, with no install step and no
`PYTHONPATH`.  `examples/quickstart.py` and the benchmark runners put `src/` on
`sys.path` themselves, pytest gets it from `pyproject.toml`'s
`pythonpath = ["src"]`, and `validation/compare.py` now does the same at import
(`sys.path.insert(0, <repo>/src)` above its `phyloslicer` imports), so the
README's bare `python validation/compare.py` is literally correct.  Measured on
8 October 2026 with the variable unset: `python validation/compare.py` exits 0
and prints `verdict: PASS`.  An earlier revision of this file claimed the
comparator was the one script that still needed `PYTHONPATH=src` or an install;
that was read from a stale snapshot and the claim is withdrawn.

`validation/compare.py` is the one that matters for the manuscript's
hash-registered-reference claim: it reads `validation/reference/manifest.txt`,
whose entries name bare filenames (the comparator resolves them relative to the
`reference/` directory, not to any recorded path), verifies each SHA-256 before
comparing, and refuses unregistered references.  Nothing in that chain refers to
the author's machine, so it survives anonymisation and re-runs green in the
archive — which is what the "reviewers can check the hashes themselves" sentence
in the paper requires.

The in-archive execution tests are the slowest part of
`tests/test_review_archive.py`: measured, the module runs 19 tests in 26 s, of
which roughly 15 s is the unpacked suite plus the example and the comparator.
That is well inside the 3-minute budget, so nothing is skipped by default; on a
machine where even that is unwelcome, set `PHYLOSLICER_SKIP_ARCHIVE_RUN=1` and
the three run tests skip while the build, scrub, provenance, reproducibility and
guard tests still execute.  They are *not* skipped in CI.

## For the editors

Author identity is withheld **only** for the duration of review, and only in the
archive.  It is carried by the manuscript title page and by the version of record
in the repository: `CITATION.cff` (names, e-mails, ORCIDs, affiliation), `LICENSE`
(`Copyright (c) 2026, Zichao Zeng and Yinzhao Wang`), `.zenodo.json`, and the README author block
of the public release at the tagged commit `v0.1.0`.  Nothing in the anonymised
archive is a permanent change: the packager writes scrubbed bytes into the zip
and never edits a source file, so the deposit that Zenodo and GitHub carry on
acceptance still names its authors, and the archive's own metadata fields read
`Withheld` rather than a placeholder for unfinished work.  The provenance sidecar
names the exact commit the archive was built from so the two can be reconciled,
and the archive's `PEER_REVIEW_NOTICE.md` / `PEER_REVIEW_NOTICE_CN.md` explain the
tokens to reviewers without explaining who they are.
