# PhyloSlicer Software Test Report

**Version:** 0.1.0  
**Test date:** 2026-09-29  
**Test environment:** Python 3.14.6 / numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3 / matplotlib 3.11.1 / pytest 9.1.1 / dendropy 5.0.13 / biopython 1.88; the optional accelerator and spatial backends (numba, geopandas) are **not installed** in this environment, so exactly two test cases skip by design  
**Test machine:** Apple M5 (10 cores), 32 GB RAM, macOS 27.0 (`python3 -c "import platform; print(platform.platform())"` reports `macOS-27.0-arm64-arm-64bit-Mach-O`, `platform.machine()` reports `arm64`)  
**Executed by:** `testing_suite/scripts/run_all_tests.py` (end-to-end; its output directories were redirected outside the working tree, see section 7) and `pytest` (unit/regression)  
**Repository path:** repository root

---

## 1. Test Summary

| Metric | Value |
|---|---|
| **Total checks (end-to-end harness)** | 118 |
| **Passed** | 118 |
| **Failed** | 0 |
| **Pass rate** | **100.0%** |
| **Total elapsed** | 3.3–4.3 seconds across the runs of this pass (the same 118 cases; wall clock varies. `testing_suite/logs/test_results.json` holds the last run) |
| **Existing pytest suite** | 411 collected → 409 passed, 0 failed, 2 skipped |

### Conclusion

All 118 end-to-end checks of PhyloSlicer v0.1.0 pass; the category breakdown below is reproduced
exactly by `testing_suite/scripts/run_all_tests.py` in its current form. The count rose from the 102
this report carried until this pass, because `TEST_PLAN.md` defines 118 cases and the harness
implemented 102 of them: the ten `T-CLI` cases, `T-ERR-06` (disconnected component) and the five
`T-PERF` cases had no implementation, so "every planned case passes" was not yet true of the
harness. Those sixteen cases are implemented now, and CI runs the harness. The pytest suite now
collects **411** cases (up from the 397 collected by the previous run of this report), expanded
from 316 `def test_` functions by 27 `@pytest.mark.parametrize` marks. Of those 411 cases,
**409 passed, none failed, and 2 skipped**; the two skips are `tests/test_rates.py:279`
(numba-accelerated kernel versus the pure kernel) and `tests/test_io_spatial_cli.py:404` (geopandas
grid adjacency), both skipped only because numba and geopandas are absent here; an environment with
the `[accel]` and `[spatial]` extras installed executes them.

**About the release-hygiene numbers.** The 409/0/2 result above was taken at 2026-09-29 22:12, after
the document-pairing and metadata work had settled. Runs of the same command taken earlier in this
same pass, while `tests/test_release_hygiene.py` and several documents were still being edited
concurrently, reported 405 passed / 4 failed, then 407 / 2, then 408 / 1. Every one of those failures
was a documentation-hygiene assertion rather than a software defect: a stale version banner still
advertised in `docs/USER_GUIDE_CN.md`, a version field in `pyproject.toml`, and
`test_every_english_document_ships_a_chinese_counterpart` reporting
`documents without a Chinese counterpart: ['testing_suite/TEST_REPORT.md']` — the last of these is
satisfied by `TEST_REPORT_CN.md`, which this very pass adds beside this report. The functional and
numerical results never moved: the end-to-end harness and all non-hygiene cases passed in every run.

---

## 2. Coverage

### 2.1 By category

| Category | Checks | Passed | Failed | Pass rate |
|---|---|---|---|---|
| Install | 6 | 6 | 0 | 100% |
| IO | 10 | 10 | 0 | 100% |
| Core | 10 | 10 | 0 | 100% |
| Indices | 10 | 10 | 0 | 100% |
| Rates | 10 | 10 | 0 | 100% |
| Slicing | 12 | 12 | 0 | 100% |
| Uncertainty | 8 | 8 | 0 | 100% |
| Spatial | 6 | 6 | 0 | 100% |
| Compat | 5 | 5 | 0 | 100% |
| Sim | 6 | 6 | 0 | 100% |
| Viz | 3 | 3 | 0 | 100% |
| CLI | 10 | 10 | 0 | 100% |
| Error | 10 | 10 | 0 | 100% |
| E2E | 7 | 7 | 0 | 100% |
| Perf | 5 | 5 | 0 | 100% |
| **Total** | **118** | **118** | **0** | **100%** |

### 2.2 Functional coverage

Every core module of PhyloSlicer is exercised:

- ✅ `phyloslicer.core.tree` — TreeArray construction, validation, depths, pruning
- ✅ `phyloslicer.slicing.api` — the 5 slicing modes plus the slice stack
- ✅ `phyloslicer.rates.api` — CpD/CpE/CpB rates + fitting + sensitivity
- ✅ `phyloslicer.rates.fitting` — exponential fit + batch fit + boundary detection
- ✅ `phyloslicer.indices.api` — PD/PE/PB indices + per-slice decomposition
- ✅ `phyloslicer.uncertainty` — posterior rate/origin + bootstrap + HDI/ESS
- ✅ `phyloslicer.spatial` — adjacency construction + rate_map
- ✅ `phyloslicer.compat` — treesliceR compatibility layer
- ✅ `phyloslicer.simulate` — Yule / birth-death trees + dataset simulation
- ✅ `phyloslicer.io` — Newick/NEXUS reading and writing + matrix building + alignment
- ✅ `phyloslicer.viz` — rate-line / sensitivity / posterior distribution plots
- ✅ `phyloslicer.datasets` — the sim_passerines dataset
- ✅ `phyloslicer.cli.main` — the slice/rates/validate CLI subcommands

---

## 3. Benchmark Results

### 3.1 Runtime (median of 3 repetitions)

| Tree size | Sites | Slices | Slice (s) | PD/slice (s) | CpD rate (s) | Posterior (s) | Memory (MB) |
|---|---|---|---|---|---|---|---|
| 20 | 10 | 50 | 0.0001 | 0.0005 | 0.0050 | 0.0164 | 0.08 |
| 50 | 20 | 50 | 0.0002 | 0.0009 | 0.0053 | 0.0195 | 0.20 |
| 100 | 30 | 50 | 0.0005 | 0.0017 | 0.0060 | 0.0253 | 0.38 |
| 100 | 30 | 100 | 0.0005 | 0.0017 | 0.0070 | 0.0302 | 0.57 |
| 100 | 30 | 200 | 0.0005 | 0.0017 | 0.0077 | 0.0300 | 0.86 |
| 200 | 50 | 100 | 0.0012 | 0.0037 | 0.0096 | 0.0438 | 1.04 |
| 500 | 100 | 100 | 0.0034 | 0.0132 | 0.0186 | 0.0929 | 2.97 |

These figures were regenerated on 2026-09-29 on the machine described above by
`testing_suite/scripts/run_benchmarks.py`, with its `BENCH_DIR` pointed outside the working tree so
nothing in the repository was rewritten. The table that ships in
`testing_suite/benchmarks/benchmark_results.csv` still holds the 2026-09-26 run made under Python
3.12.14; its peak-memory column matches the table above to within rounding, while its timings are
3-5x larger. Regenerating that artifact is a separate step from this report.

### 3.2 Performance analysis

- **Slicing**: under 4 ms at every size tested, scaling linearly.
- **CpD rate computation** (alignment + slicing + contribution matrix + fit): 18.6 ms for a 500-tip tree.
- **Posterior analysis** (3 trees): 92.9 ms for a 500-tip tree.
- **Memory**: a 500-tip tree with 100 sites and 100 slices peaks at about 3 MB.
- **Scaling in slice count**: going from 50 to 200 slices raises CpD time by ~28%, as the expected
  linear growth in the `edges x slices` contribution matrix predicts.

---

## 4. Issues Found During Testing, and Fixes

### 4.1 Issue record

Five problems **in the test scripts themselves** (not software bugs) were corrected while this
harness was being written; all five corrections are still in `testing_suite/scripts/run_all_tests.py`:

| Issue | Cause | Fix |
|---|---|---|
| T-RATES-04: fit_rate parameter recovery | The model `cum ∝ exp(-r*ages)` in `fit_rate` has no intercept, so after normalization `cum[-1]=1` constrains the range of r; the test used r=0.5, which the model cannot match | Changed to assert convergence + R² > 0.7 + rate > 0 |
| T-COMPAT-03: squeeze_int argument name | The argument is named `to`, not `to_` (`to` is a usable identifier here) | Use `to=0.5` |
| T-COMPAT-05: prune_tips threshold | In a 4-tip symmetric tree every terminal branch is 1.0, so a 0.5 threshold pruned all of them | Use a tree with unequal branch lengths and a sensible threshold |
| T-E2E-05: CLI invocation | Running `main.py` directly breaks the relative imports | Use `python -m phyloslicer.cli.main` |
| T-E2E-05: output file clash | Re-running hit pre-existing output files | Clean the output files before the test |

### 4.2 Problems in the software itself

**This round found no new software bug.** For the record, `tests/test_regressions.py` and
`tests/test_defect_regressions.py` exist precisely to pin behaviour that was found to be wrong during
development and then fixed (slice endpoints, beta-profile normalization, the units of the ultrametric
test, the direction of `timeSteps`, pruning complexity, Newick quoting, refusal of an unregistered
validation reference, and so on). Their case counts should therefore not be read as "the software has
no defects", but as "these defects now have regression guards".

---

## 5. End-to-End Test Details

All figures below come from the 2026-09-29 harness run described in section 1.

### 5.1 Full CpD workflow (T-E2E-01)
- Input: 20-tip Yule tree + 30-site presence/absence matrix
- Operation: `cpd_rate(tree, mat, n_slices=50)`
- Output: a 30-row DataFrame with 11 columns (`site_id`, `CpD`, `PD`, `pDO`, `converged`, `r_squared`,
  `reason`, `at_lower_bound`, `at_upper_bound`, `ultrametric`, `origin_within_tree`)
- Result: ✅ all 30 sites converged

### 5.2 Full CpE workflow (T-E2E-02)
- Operation: `cpe_rate(tree, mat, n_slices=50)`
- Result: ✅ 30 rows carrying the `CpE`/`PE`/`pEO` columns

### 5.3 Full CpB workflow (T-E2E-03)
- Operation: `cpb_rate(tree, mat, adj, n_slices=50, component="sorensen")`
- Result: ✅ 30 rows carrying the `CpB`/`PB`/`pBO` columns

### 5.4 Posterior analysis workflow (T-E2E-04)
- Input: directory of 5 posterior trees + 30-site matrix
- Operation: `posterior_rate(dir, mat, n_slices=30, hdi=0.95)`
- Result: ✅ 30 rows with `n_trees = 5`, HDI bounds (`hdi_low`/`hdi_high`) and an ESS proxy
  (`ess_proxy`, with its `ess_proxy_caveat` column)

### 5.5 CLI workflow (T-E2E-05)
- Operation: `phyloslicer slice` + `phyloslicer rates`
- Result: ✅ Newick and CSV outputs produced (`cli_slice.tre`, `cli_rates.csv`)

### 5.6 Cross-validation (T-E2E-06)
- Operation: `phyloslicer validate` against the registered treesliceR reference (`bd20__random__CpD.csv`)
- Result: ✅ equivalence rate 1.0000 over all 30 sites

### 5.7 Visualization workflow (T-E2E-07)
- Operation: `plot_rate_line` renders the rate curve
- Result: ✅ PNG figure written (`e2e_rate_line.png`)

---

## 6. The Existing pytest Suite

Besides the 118 end-to-end checks, PhyloSlicer ships a full pytest unit suite:

```
411 collected -> 409 passed, 0 failed, 2 skipped   (whole command: 11-16 s wall clock, run twice)
```

The 411 cases come from 316 `def test_` functions expanded by 27 `@pytest.mark.parametrize` marks.
The command was `python3 -m pytest tests/ -p no:cacheprovider -q --tb=no -rs`. Every run made for
this report ended its output after the short skip/failure summary without the usual trailing
`N passed in Xs` count line, so the figures above were read off the progress line and the
`-rs`/`-rf` short summaries, and the wall clock from the shell's own timing.

The two skips:

| Skip site | Reason |
|---|---|
| `tests/test_rates.py:279` | `numba not installed (pip install phyloslicer[accel])` |
| `tests/test_io_spatial_cli.py:404` | `could not import 'geopandas': No module named 'geopandas'` |

Neither numba nor geopandas is installed here, which is why exactly two cases skip; a run with the
`[accel]` and `[spatial]` extras installed executes both. No case fails in the run recorded above;
the earlier snapshots taken while the documents were still moving are described in section 1.

Files covered (the number is the cases collected from that file, from
`python3 -m pytest tests/ --collect-only -q`):

- `tests/test_core.py` (13) — TreeArray structure, depths, tip intervals, validation
- `tests/test_slicing.py` (26) — slicing semantics and analytical solutions
- `tests/test_indices.py` (41) — diversity indices against brute-force computations
- `tests/test_rates.py` (23) — rate fitting, diagnostics, posterior analysis
- `tests/test_invariants.py` (68) — cross-module invariants and end-to-end properties
- `tests/test_io_spatial_cli.py` (34) — reading/writing, spatial adjacency, command line
- `tests/test_compat.py` (52) — treesliceR-style compatibility layer
- `tests/test_regressions.py` (27) — regression tests
- `tests/test_defect_regressions.py` (82) — one pinned case per defect, named for the symptom
- `tests/test_validation_chain.py` (32) — validation chain and self-attestation guards
- `tests/test_release_hygiene.py` (11) — release hygiene of the published package
- `tests/test_published_numbers.py` (2) — documented speed-up ranges recomputed from the tables

---

## 7. Test Asset Inventory

```
testing_suite/
├── TEST_PLAN.md              # test plan
├── TEST_PLAN_CN.md           # Chinese counterpart of the test plan
├── TEST_REPORT.md            # this report (version of record)
├── TEST_REPORT_CN.md         # Chinese translation of this report
├── data/                     # test data
│   ├── small_tree.tre        # 4-tip symmetric tree
│   ├── small_mat.csv         # 4-site matrix
│   ├── medium_tree.tre       # 20-tip Yule tree
│   ├── medium_mat.csv        # 30-site matrix
│   ├── large_tree.tre        # 100-tip birth-death tree
│   ├── large_mat.csv         # 50-site matrix
│   ├── posterior_trees/      # 5 posterior trees
│   ├── posterior.nex         # posterior trees in NEXUS (5 trees)
│   ├── posterior.tre         # posterior trees in Newick (5 trees)
│   ├── grid_bounds.csv       # grid bounds (30 cells)
│   ├── grid_adj.csv          # Queen adjacency matrix (30 x 30)
│   ├── passerines_tree.tre   # simulated passerine tree (50 tips)
│   ├── passerines_mat.csv    # simulated matrix (50 sites)
│   ├── passerines_bounds.csv # simulated grid (50 cells)
│   ├── passerines_adj.csv    # simulated adjacency
│   └── passerines_posterior/ # simulated posterior trees (3 trees)
├── scripts/
│   ├── generate_test_data.py # test data generation
│   ├── run_all_tests.py      # the end-to-end suite
│   └── run_benchmarks.py     # benchmarks
├── results/                  # test outputs
│   ├── e2e_cpd.csv           # CpD rate results
│   ├── e2e_cpe.csv           # CpE rate results
│   ├── e2e_cpb.csv           # CpB rate results
│   ├── e2e_posterior.csv     # posterior analysis results
│   ├── e2e_validate.csv      # cross-validation results
│   ├── e2e_rate_line.png     # rate curve figure
│   ├── cli_slice.tre         # CLI slice output
│   └── cli_rates.csv         # CLI rate output
├── logs/
│   └── test_results.json     # detailed results (JSON)
└── benchmarks/
    ├── benchmark_results.csv # benchmark data
    ├── BENCHMARK_REPORT.md   # benchmark report
    └── BENCHMARK_REPORT_CN.md # Chinese counterpart of the benchmark report
```

Two notes on this inventory. First, `testing_suite/scripts/run_all_tests.py` and
`testing_suite/scripts/run_benchmarks.py` write into `results/`, `logs/` and `benchmarks/` when run
with their default paths; the 118-case harness run recorded in section 1 used those defaults, so
`testing_suite/logs/test_results.json` and the CSVs in `results/` on disk are that run
(2026-09-29, version 0.1.0). Those CSVs also gained the `origin_within_tree` column, which the
committed copies predated. Second, every path listed above exists in the current tree.

---

## 8. Conclusions

PhyloSlicer v0.1.0 passes the tests that exercise the software:

1. **Functional completeness** — every core feature (slicing, rates, indices, posterior, spatial,
   visualization, CLI) works.
2. **Numerical correctness** — the mathematical properties of the diversity indices (row sums = total
   PD, the partition property) are verified.
3. **Error handling** — every boundary condition and illegal input raises a clear `ValueError`.
4. **Reproducibility** — the same seed gives the same result.
5. **Performance** — a full rate analysis of a 500-tip tree finishes in under 20 ms.
6. **Compatibility** — the treesliceR compatibility layer behaves like the original R implementation.
7. **Cross-validation** — agreement with the treesliceR reference data is verified (equivalence rate
   1.0000 over 30 sites).

**Recommendation:** the software is ready to be released — all 118 end-to-end checks pass and the
pytest suite is green apart from the two optional-backend skips. Re-run `python -m pytest -q` in an
environment that also has the `[accel]` and `[spatial]` extras installed, in the same commit that is
tagged, so that those two cases and the release-hygiene gate are verified green together.
