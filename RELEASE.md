# Release procedure

PhyloSlicer version 0.1.0 · BSD-3-Clause

A release here is the version that the documentation, the citation metadata and
the archived validation references all describe at the same time. Cut it in this
order.

## 1. Fix the version of record

`src/phyloslicer/__init__.py` is the only place a version literal is allowed;
`pyproject.toml` reads it dynamically. Update every file that advertises it:

- `README.md` and `README_CN.md` banner
- `docs/USER_GUIDE.md` and `docs/USER_GUIDE_CN.md` banner
- `CONTRIBUTING.md` and `CONTRIBUTING_CN.md` banner
- `CHANGELOG.md` and `CHANGELOG_CN.md` - add a dated `## X.Y.Z — YYYY-MM-DD` entry
- `CITATION.cff` - `version:` and `date-released:`
- `.zenodo.json` - `version` and `publication_date`

`python -m pytest tests/test_release_hygiene.py` fails if any of them drifts.

## 2. Run every gate that CI runs

```bash
ruff check src tests benchmarks examples testing_suite validation
ruff format --check src tests benchmarks examples testing_suite validation
python -m pytest -q
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python - <<'PY'
import hashlib, re
from pathlib import Path
man = Path("validation/reference/manifest.txt").read_text(encoding="utf-8")
entries = re.findall(r"\[file ([^\]]+)\]\nsha256 = ([0-9a-f]{64})", man)
bad = [n for n, h in entries
       if hashlib.sha256((Path("validation/reference") / n).read_bytes()).hexdigest() != h]
print(f"{len(entries)} registered references, {len(bad)} hash mismatches")
assert not bad and entries
PY
python benchmarks/recompute_ratios.py
```

`validation/compare.py` rewrites `validation/reference/equivalence_report.csv`.
Commit that file only if the report is meant to move with this release; the
pinned references themselves must never change without re-running the R
generator in `validation/r/`.

## 3. Tag and publish

```bash
git tag -a "v0.1.0" -m "PhyloSlicer 0.1.0"
git push origin main --follow-tags
```

Then create the GitHub release for the tag, using the matching `CHANGELOG.md`
section as the notes.

## 4. Mint and record the DOI

Enable the Zenodo GitHub integration once, at
<https://zenodo.org/account/settings/github/>, and turn on this repository.
Zenodo then picks up the tagged release: `.zenodo.json` pre-fills the draft,
and its `version` and `publication_date` are overwritten from the tag metadata,
which is expected. Do not put a `conceptdoi` field in `.zenodo.json` - Zenodo
mints the DOIs, and writing one there beforehand produces a conflict.

Two DOIs come back: a **version DOI** (one per release) and a
**version-agnostic concept DOI** (stable, cites every release). Record the
concept DOI:

- add an `identifiers:` entry to `CITATION.cff` holding the concept DOI, and
  cite it beside the repository URL in `README.md` / `README_CN.md`
- never commit a placeholder DOI: a reserved value that resolves to nothing is
  worse than no DOI, because someone will cite it in good faith
- use the version DOI for the release badge, if you want one:
  `https://zenodo.org/badge/DOI/<version-doi>.svg`
- do not publish with that placeholder still in place; it exists so the
  metadata can be complete before the first tag does

## 5. Confirm from a clean environment

```bash
python -m venv /tmp/phyloslicer-check && source /tmp/phyloslicer-check/bin/activate
pip install "phyloslicer @ git+https://github.com/ZengZichao/PhyloSlicer@v0.1.0"
python -c "import phyloslicer; print(phyloslicer.__version__)"
phyloslicer --version
```

The release is done when the installed version, the tag, the changelog entry,
the citation metadata and the DOI all name the same build.
