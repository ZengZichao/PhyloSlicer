# Contributing

Thank you for working on PhyloSlicer. This project is a research package, so the
bar for a change is that it stays *checkable*: every documented number must
follow from a file in the repository, and every claim about agreement with the
reference implementation must be reproducible.

Version 0.1.0 · BSD-3-Clause · Python ≥ 3.10

## Setting up a development environment

```bash
git clone https://github.com/ZengZichao/PhyloSlicer.git
cd PhyloSlicer
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,io,spatial,accel]"
pre-commit install
```

`[accel]` (numba) and `[spatial]` (geopandas) are optional backends. Their tests
skip when the backend is missing, so a full local run should install both; CI
tests the accelerated and the pure kernel against each other on purpose.

## The gates a change must pass

```bash
ruff check src tests benchmarks examples testing_suite validation
ruff format --check src tests benchmarks examples testing_suite validation
python -m pytest -q
python testing_suite/scripts/run_all_tests.py
python examples/quickstart.py --sites 24 --species 60 --trees 3
python validation/compare.py
python benchmarks/recompute_ratios.py
```

`.pre-commit-config.yaml` runs the two ruff gates on commit; CI runs all seven on
Python 3.10-3.13 on Ubuntu and macOS.

## Conventions

- **One source of truth for the version.** `src/phyloslicer/__init__.py` holds
  `__version__`; `pyproject.toml` reads it dynamically. Do not add a second
  literal. `tests/test_release_hygiene.py` fails if a banner or a citation
  field disagrees with it.
- **Numerical changes need a reference.** A new or changed index/rate path is
  scored against treesliceR output registered in
  `validation/reference/manifest.txt` (SHA-256, shape, column pattern,
  orientation). An unregistered reference is refused rather than trusted.
- **Never quote a benchmark number you have not regenerated.** The four
  protocols in `benchmarks/` are not interchangeable; read
  `benchmarks/README.md`, keep the `source` column honest, and report medians
  (not minimums) for cross-language comparisons.
- **Diagnostics are part of the API.** A fit that fails returns a row with
  `converged=False` and a human-readable `reason`; do not raise away a
  degeneracy that a reader has to be able to see.
- **Docs ship in pairs.** The English file (no suffix) is the version of record;
  every `X.md` has an `X_CN.md` counterpart with the same section structure. If
  you change the English text, change the Chinese one in the same pull request.
- **No personal paths.** Absolute local paths anywhere in the tree fail
  `tests/test_release_hygiene.py::test_no_personal_paths_or_names_anywhere_in_the_tree`.

## Proposing a change

Open an issue first for anything that changes public behaviour, an index
definition, or a tolerance in `validation/`. Pull requests should state which
gate output you ran, link the issue, and add a regression test that fails without
the fix (see `tests/test_defect_regressions.py` for the house style: one pinned
test per defect, named for the symptom).

## Releasing

See [RELEASE.md](RELEASE.md) · Chinese: [RELEASE_CN.md](RELEASE_CN.md).
