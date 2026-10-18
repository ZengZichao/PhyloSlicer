# Continuous-integration workflows (reviewer mirror)

Ecography's Software Notes instructions "strongly recommend using continuous
integration services that trigger a test run every time the source code is
modified". The live definitions live in `.github/workflows/`, which the review
archive cannot ship as-is: `.github/CODEOWNERS` in that directory names the
repository owner and would defeat double-blind review, so the whole directory is
excluded by the packager. This file mirrors the workflow definitions themselves,
which carry no author identity, so that a reviewer can read exactly what CI runs.

Mirrored from the working tree on 2026-10-18. `tests/test_ci_mirror.py` fails if a workflow changes without refreshing this file.

## `.github/workflows/build.yml`

```yaml
name: build

# 构建门禁：sdist 与 wheel 必须能从干净检出构建，且元数据通过
# twine check。打 v* 标签时把产物直接挂到对应 GitHub Release，
# 发布不再依赖本地手动 build。

on:
  push:
    branches: [main]
    tags: ["v*"]
  pull_request:

permissions:
  contents: read

concurrency:
  group: build-${{ github.ref }}
  cancel-in-progress: true

jobs:
  sdist-wheel:
    runs-on: ubuntu-latest
    permissions:
      # 挂载产物到 GitHub Release 需要写权限；顶层保持 read，仅本作业放开
      contents: write
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - name: Build sdist and wheel
        run: |
          python -m pip install --upgrade pip build twine
          python -m build
      - name: Check distribution metadata
        run: python -m twine check dist/*
      - name: Upload build artifacts
        uses: actions/upload-artifact@v4
        with:
          name: dist
          path: dist/
      - name: Attach distributions to the release
        if: startsWith(github.ref, 'refs/tags/v')
        env:
          GH_TOKEN: ${{ github.token }}
        run: gh release upload "$GITHUB_REF_NAME" dist/* --clobber

```

## `.github/workflows/codeql.yml`

```yaml
name: codeql

# GitHub 代码扫描：对公开仓库免费，覆盖 Python 的安全类缺陷
# （注入、不安全反序列化、路径穿越等）。结果出现在仓库的
# Security → Code scanning 下，不阻塞 tests 门禁。

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]
  schedule:
    - cron: "23 4 * * 1"

permissions:
  contents: read
  security-events: write

concurrency:
  group: codeql-${{ github.ref }}
  cancel-in-progress: true

jobs:
  analyze:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: github/codeql-action/init@v3
        with:
          languages: python
      # 纯解释型语言无需 autobuild；依赖由 init 的 extract 按需解析
      - uses: github/codeql-action/analyze@v3
        with:
          category: "/language:python"
```

## `.github/workflows/tests.yml`

```yaml
name: tests

on:
  # 只在落上 main 时跑 push 事件：Dependabot/功能分支由 pull_request 覆盖，
  # 避免同一个 PR 的分支被 push 和 pull_request 各跑一遍全矩阵。
  push:
    branches: [main]
  pull_request:
  schedule:
    - cron: "17 3 * * 1"

permissions:
  contents: read

concurrency:
  group: tests-${{ github.ref }}
  # 同一分支连续推送时取消已被超越的旧跑，省 16 个矩阵 job 的排队时间
  cancel-in-progress: true

jobs:
  suite:
    runs-on: ${{ matrix.os }}
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, macos-latest]
        # every Python version pyproject.toml claims to support (3.10 - 3.14)
        python-version: ["3.10", "3.11", "3.12", "3.13", "3.14"]
        # the accelerated and the pure kernel must agree, so test both
        accel: ["with-accel", "no-accel"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install the project with its optional extras
        # installing the *project* (not a hand-written dependency list) is what
        # makes the [io] adapters and the [spatial] geopandas path, and their
        # importorskip-guarded tests, actually run here
        run: |
          python -m pip install --upgrade pip
          if [ "${{ matrix.accel }}" = "with-accel" ]; then
            python -m pip install -e ".[dev,io,spatial,accel]"
          else
            python -m pip install -e ".[dev,io,spatial]"
          fi
      - name: Check the environment matches the matrix
        run: |
          python - <<'PY'
          import importlib
          import importlib.util
          import sys

          want_numba = "${{ matrix.accel }}" == "with-accel"
          required = ["numpy", "scipy", "pandas", "matplotlib", "pytest",
                      "dendropy", "Bio", "geopandas"]
          missing = []
          for module in required:
              try:
                  importlib.import_module(module)
              except ImportError:
                  missing.append(module)
          if missing:
              sys.exit(f"extras are not installed, so their tests would be "
                       f"silently skipped: missing {missing}")
          from phyloslicer._accel import HAS_NUMBA

          has_numba = importlib.util.find_spec("numba") is not None
          if has_numba != want_numba or HAS_NUMBA != want_numba:
              sys.exit(f"accel axis says {want_numba}, but numba is importable="
                       f"{has_numba} and phyloslicer._accel.HAS_NUMBA={HAS_NUMBA}")
          PY
      - name: Lint and format
        # every directory that ships in the release is gated, not just src
        run: |
          ruff check src tests benchmarks examples testing_suite validation
          ruff format --check src tests benchmarks examples testing_suite validation
      - name: Test
        run: python -m pytest -q --durations=10
      - name: Build the anonymised review archive and check no personal path or author identity survives
        # the submission checklist tells the Ecography editors that CI packs the
        # double-blind review archive through the packager and fails the job if any
        # personal path or author identity survives.  This step is that guard: the
        # packager scrubs, then re-reads its own bytes and exits non-zero on any
        # surviving hit, and the test re-checks the archive and unpacks and runs it.
        # No dependency setup is needed beyond the checkout above: the packager is
        # stdlib-only (build-time tooling must not add a runtime dependency) and the
        # test drives the interpreter that the matrix already installed.
        run: |
          python tools/make_review_archive.py --verbose
          python tools/make_review_archive.py --format tar.gz
          python -m ruff check tools tests
          python -m ruff format --check tools tests
          python -m pytest -q --durations=10 tests/test_review_archive.py
      - name: Run the end-to-end test harness
        # every case the test plan defines, including the CLI and performance
        # budgets; the harness exits nonzero if any of them fails
        run: python testing_suite/scripts/run_all_tests.py
      - name: Run the documented out-of-the-box example
        # a reader must be able to run the example without acquiring anything
        # from outside the repository; this step fails if the example and the
        # docs drift apart
        run: python examples/quickstart.py --sites 24 --species 60 --trees 3
      - name: Re-run the cross-validation and check its provenance
        run: |
          python validation/compare.py
          python - <<'PY'
          import hashlib, re, sys
          from pathlib import Path
          man = Path("validation/reference/manifest.txt").read_text(encoding="utf-8")
          entries = re.findall(r"\[file ([^\]]+)\]\nsha256 = ([0-9a-f]{64})", man)
          bad = [n for n, h in entries
                 if hashlib.sha256((Path("validation/reference") / n).read_bytes()).hexdigest() != h]
          print(f"manifest: {len(entries)} registered references, {len(bad)} hash mismatches")
          sys.exit(1 if bad or not entries else 0)
          PY
      - name: Recompute every documented speed-up ratio from the archived tables
        # the documentation quotes ranges, not raw timings; this fails if a
        # quoted range no longer follows from the CSVs it is said to come from
        run: python benchmarks/recompute_ratios.py
      - name: Check the release metadata is consistent
        run: |
          python - <<'PY'
          import re, sys
          from pathlib import Path
          import phyloslicer

          version = phyloslicer.__version__
          for name in ("CITATION.cff", ".zenodo.json", "README.md", "README_CN.md"):
              text = Path(name).read_text(encoding="utf-8")
              assert version in text, f"{name} does not mention {version}"
          assert re.search(rf'^version: {re.escape(version)}$',
                           Path("CITATION.cff").read_text(encoding="utf-8"), re.M)
          print(f"release metadata consistent at {version}")
          PY
```


## `.github/workflows/publish.yml`

```yaml
name: publish

# PyPI 发布门禁：打 v* 标签时从干净检出构建 sdist 与 wheel，校验元数据后
# 通过 Trusted Publishing（OIDC）上传 PyPI——仓库里不存任何 API token。
# 首次使用前需在 PyPI（Account settings → Publishing）登记 pending
# publisher：owner=ZengZichao、repository=PhyloSlicer、workflow=publish.yml、
# environment=pypi；登记后重新运行本工作流即可。

on:
  push:
    tags: ["v*"]
  workflow_dispatch:

permissions:
  contents: read

jobs:
  pypi-publish:
    name: Build and publish to PyPI
    runs-on: ubuntu-latest
    environment: pypi
    permissions:
      # OIDC claim for Trusted Publishing; PyPI rejects the upload without it
      id-token: write
      contents: read
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.12"
      - name: Verify tag matches the package version
        if: startsWith(github.ref, 'refs/tags/v')
        run: |
          TAG_VERSION="${GITHUB_REF_NAME#v}"
          PKG_VERSION="$(python - <<'EOF'
          import pathlib
          import re

          text = pathlib.Path("src/phyloslicer/__init__.py").read_text(encoding="utf-8")
          print(re.search(r'__version__ = "([^"]+)"', text).group(1))
          EOF
          )"
          if [ "$TAG_VERSION" != "$PKG_VERSION" ]; then
            echo "tag v$TAG_VERSION does not match package version $PKG_VERSION" >&2
            exit 1
          fi
      - name: Build sdist and wheel
        run: |
          python -m pip install --upgrade pip build twine
          python -m build
      - name: Check distribution metadata
        run: python -m twine check dist/*
      - name: Publish to PyPI
        uses: pypa/gh-action-pypi-publish@release/v1
```
