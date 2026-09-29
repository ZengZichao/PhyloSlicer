"""Release-hygiene invariants for the published package.

Each test here pins a property that a public release must keep true:

* the tree carries no build artefact that embeds a personal absolute path, and
  ``.gitignore`` keeps such artefacts out of the repository;
* no source file may cite an internal review item, plan section or document
  that is not part of the release;
* the version, licence and release date are stated in the package and agree
  across README, manuals, changelog and citation metadata;
* the metadata names one public repository, consistently, with no placeholder
  placeholder left behind;
* the documented out-of-the-box example actually runs;
* the origin-time range flag fires on the bundled dataset it was added for.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT_SUFFIXES = {".py", ".md", ".toml", ".yml", ".yaml", ".R", ".csv", ".txt", ".drawio"}

SELF = Path(__file__).name  # this file quotes the patterns it forbids

# every file that advertises a version banner or a version field, and must
# therefore name the version of record literally
VERSION_BEARING = (
    "README.md",
    "README_CN.md",
    "docs/USER_GUIDE.md",
    "docs/USER_GUIDE_CN.md",
    "CHANGELOG.md",
    "CHANGELOG_CN.md",
    "CONTRIBUTING.md",
    "CONTRIBUTING_CN.md",
    "RELEASE.md",
    "RELEASE_CN.md",
    "CITATION.cff",
    ".zenodo.json",
)
# pyproject.toml is scanned for drift but must NOT carry a literal: it reads the
# version dynamically from the package (see test_defect_regressions.TestMetadata)
VERSION_SCANNED = VERSION_BEARING + ("pyproject.toml",)


def _tracked_text_files():
    for p in sorted(ROOT.rglob("*")):
        if p.name == SELF:
            continue
        if not p.is_file() or p.suffix not in TEXT_SUFFIXES:
            continue
        if any(part in {".git", "__pycache__", ".ruff_cache"} for part in p.parts):
            continue
        yield p


def test_gitignore_excludes_artifacts_that_carry_absolute_paths():
    """Caches and build output embed the builder's home path; keep them out."""
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    for entry in (
        "__pycache__/",
        ".ruff_cache/",
        ".mypy_cache/",
        ".pytest_cache/",
        "build/",
        "dist/",
    ):
        assert entry in gitignore, f"{entry} must stay out of the repository"


def test_no_personal_paths_or_names_anywhere_in_the_tree():
    """A public repository is no place for a local filesystem layout."""
    home = str(Path.home())
    # assembled from parts so this file is not itself a hit for the check
    users_dir = "/" + "Users" + "/"
    hits = []
    for p in _tracked_text_files():
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if home in text or users_dir in text or "C:" + "\\" + "Users" in text:
            hits.append(str(p.relative_to(ROOT)))
    assert not hits, f"personal absolute paths found in: {hits}"


def test_no_dangling_internal_references():
    """No comment may cite a review item / plan section outside this release."""
    pattern = re.compile(r"review(?:\s+items?)?\s+A\d|item A\d|plan section \d|§\d+ asks|审阅报告")
    hits = []
    for p in _tracked_text_files():
        for n, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            if pattern.search(line):
                hits.append(f"{p.relative_to(ROOT)}:{n}")
    assert not hits, f"references to documents not shipped with the release: {hits}"


def test_version_licence_and_release_date_are_stated():
    import phyloslicer

    version = phyloslicer.__version__
    assert re.fullmatch(r"\d+\.\d+\.\d+", version), version
    doc = sys.modules["phyloslicer"].__doc__ or ""
    assert "BSD-3-Clause" in doc, "the licence must be stated inside the package"

    changelog = (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    assert re.search(rf"^## {re.escape(version)} — \d{{4}}-\d{{2}}-\d{{2}}", changelog, re.M), (
        "the version of record needs a dated changelog entry"
    )

    for name in VERSION_BEARING:
        text = (ROOT / name).read_text(encoding="utf-8")
        assert version in text, f"{name} does not mention version {version}"


def test_no_document_advertises_a_stale_version():
    """Every version banner and version field equals ``phyloslicer.__version__``."""
    import phyloslicer

    version = phyloslicer.__version__
    stale = []
    for name in VERSION_SCANNED:
        for n, line in enumerate((ROOT / name).read_text(encoding="utf-8").splitlines(), 1):
            # a banner line ("Version 0.1.0 · ..." / "版本 0.1.0 · ...") or a
            # metadata field ("version: 0.1.0", '"version": "0.1.0"')
            if re.match(r"^\s*(?:\*\*)?(?:Version|版本)\s", line) or re.match(
                r"^\s*\"?version\"?\s*[:=]", line
            ):
                for found in re.findall(r"\d+\.\d+\.\d+", line):
                    if found != version:
                        stale.append(f"{name}:{n} advertises {found}")
    assert not stale, f"version drift: {stale}"


def test_every_english_document_ships_a_chinese_counterpart():
    """``X.md`` is the version of record; ``X_CN.md`` must exist beside it."""
    missing = []
    for path in sorted(ROOT.rglob("*.md")):
        if any(part in {".git", "__pycache__", ".venv", "node_modules"} for part in path.parts):
            continue
        if path.name.endswith("_CN.md"):
            continue
        counterpart = path.with_name(f"{path.stem}_CN.md")
        if not counterpart.exists():
            missing.append(str(path.relative_to(ROOT)))
    assert not missing, f"documents without a Chinese counterpart: {missing}"


def test_metadata_names_the_public_repository_consistently():
    """One repository URL, used everywhere, with real author attribution."""
    repo = "https://github.com/ZengZichao/PhyloSlicer"
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert f'Repository = "{repo}"' in pyproject
    assert 'name = "Zichao Zeng"' in pyproject and "zengzichao@sjtu.edu.cn" in pyproject

    assert f'repository-code: "{repo}"' in (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert repo in (ROOT / ".zenodo.json").read_text(encoding="utf-8")

    for name in ("README.md", "README_CN.md"):
        assert repo in (ROOT / name).read_text(encoding="utf-8"), name
    license_text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert "Zichao Zeng" in license_text, "the licence must name the copyright holder"

    # no placeholder or unattributed authorship may survive
    for name in ("pyproject.toml", "CITATION.cff", ".zenodo.json", "LICENSE"):
        assert "PhyloSlicer developers" not in (ROOT / name).read_text(encoding="utf-8"), name


def test_quickstart_example_runs_out_of_the_box():
    """The documented worked example must run with one command."""
    script = ROOT / "examples" / "quickstart.py"
    assert script.exists(), "an out-of-the-box example is part of the release"
    out = subprocess.run(
        [sys.executable, str(script), "--sites", "24", "--species", "40", "--trees", "3"],
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert out.returncode == 0, out.stderr[-800:]
    assert "converged" in out.stdout and "rate_median" in out.stdout


def test_origin_range_flag_fires_on_the_bundled_dataset():
    """The shipped example contains origins the fitted model puts before the root."""
    from phyloslicer import cpd_rate
    from phyloslicer.datasets import sim_passerines

    d = sim_passerines(n_sites=208, n_species=308, n_trees=12, seed=7)
    res = cpd_rate(d["trees"][0], d["mat"], n_slices=50)
    assert "origin_within_tree" in res.columns
    flagged = res[~res["origin_within_tree"]]
    assert len(flagged) == 12, "the bundled example's extrapolated origins are a known count"
    assert (flagged["pDO"] > d["trees"][0].root_age).all()
    assert "origin extrapolated beyond root" in " ".join(flagged["reason"])


def test_posterior_outputs_carry_the_documented_columns():
    from phyloslicer import posterior_origin, posterior_rate
    from phyloslicer.datasets import sim_passerines

    d = sim_passerines(n_sites=24, n_species=60, n_trees=3, seed=11)
    pr = posterior_rate(d["trees"], d["mat"], rate="cpd", n_slices=10)
    po = posterior_origin(d["trees"], d["mat"], rate="cpd", n_slices=10)
    for frame in (pr, po):
        assert {"n_ultrametric_draws", "origin_within_tree"} <= set(frame.columns)
    assert (pr["n_ultrametric_draws"] == pr["n_trees"]).all()


def test_rate_profile_and_sensitivity_report_ultrametricity():
    from phyloslicer import rate_profile
    from phyloslicer.datasets import sim_passerines

    d = sim_passerines(n_sites=24, n_species=60, n_trees=2, seed=5)
    out = rate_profile(d["trees"][0], d["mat"], n_slices=10, index="PD", normalize=False)
    assert out["ultrametric"] is True
