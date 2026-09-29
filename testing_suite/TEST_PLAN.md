# PhyloSlicer Software Test Plan

**Version:** 0.1.0  
**Date:** 2026-09-26  
**Environment:** conda `phyloslicer` (Python 3.12, numpy/scipy/pandas/matplotlib/numba)  
**Working path:** `testing_suite/` (relative to the repository root)

---

## 1. Test objectives

A complete, systematic end-to-end test of PhyloSlicer, covering:

1. **Installation and dependencies** — confirm the package installs and imports, and that every core/optional dependency is available.
2. **Command-line interface (CLI)** — end-to-end tests of the four subcommands `phyloslicer slice | rates | posterior | validate`.
3. **Input parsing** — Newick/NEXUS reading, presence/absence matrix construction, tree-matrix alignment.
4. **Core algorithms** — contribution matrix, slicing windows, PD/PE/PB indices, rate fitting.
5. **Output results** — CSV/Newick file correctness, field completeness, numerical consistency.
6. **Exception handling** — error messages for invalid input, boundary conditions and degenerate scenarios.
7. **Boundary conditions** — 2-tip trees, empty sites, degenerate neighbourhoods, NaN handling.
8. **Full workflow on realistic data** — use the simulated datasets to run the complete workflow from tree reading through rate calculation, posterior analysis and visualization.
9. **Performance/benchmark tests** — runtime and memory on datasets of different sizes.
10. **Reproducibility verification** — consistency of results under the same random seed, numerical stability.

---

## 2. Test environment

| Item | Version |
|---|---|
| Python | 3.12 |
| numpy | >= 1.24 |
| scipy | >= 1.11 |
| pandas | >= 2.0 |
| matplotlib | >= 3.7 |
| numba | >= 0.59 (optional acceleration) |
| dendropy | >= 5.0 (I/O adapter) |
| pytest | >= 8.0 |

---

## 3. Test case categories

### 3.1 Installation and dependency tests (T-INSTALL)

| ID | Description | Expected result |
|---|---|---|
| T-INSTALL-01 | `import phyloslicer` succeeds | no ImportError |
| T-INSTALL-02 | `phyloslicer.__version__` returns a string | non-empty string |
| T-INSTALL-03 | core dependencies numpy/scipy/pandas/matplotlib are importable | no exception |
| T-INSTALL-04 | optional dependency numba is importable | no exception |
| T-INSTALL-05 | optional dependency dendropy is importable | no exception |
| T-INSTALL-06 | CLI entry point `phyloslicer --version` prints normally | exit 0 |

### 3.2 Command-line interface tests (T-CLI)

| ID | Description | Expected result |
|---|---|---|
| T-CLI-01 | `phyloslicer slice tree.tre --rootward 5.0 --out out.tre` | exit 0, Newick output file |
| T-CLI-02 | `phyloslicer slice tree.tre --tipward 5.0 --out out.tre` | exit 0 |
| T-CLI-03 | `phyloslicer slice tree.tre --n-slices 10 --out dir/` | exit 0, 10 .tre files |
| T-CLI-04 | `phyloslicer rates tree.tre sites.csv --rate cpd --out r.csv` | exit 0, CSV output |
| T-CLI-05 | `phyloslicer rates tree.tre sites.csv --rate cpe --out r.csv` | exit 0 |
| T-CLI-06 | `phyloslicer posterior trees/ sites.csv --out p.csv` | exit 0 |
| T-CLI-07 | `phyloslicer validate tree.tre sites.csv --reference ref.csv --out v.csv` | exit 0 (equivalence rate = 1) |
| T-CLI-08 | `phyloslicer slice` without a required argument | exit 2 |
| T-CLI-09 | `phyloslicer slice tree.tre --rootward 5.0 --tipward 3.0 --out o.tre` | exit 2 (mutually exclusive arguments) |
| T-CLI-10 | `--log` writes a JSON run log | log file contains version/elapsed_s |

### 3.3 Input parsing tests (T-IO)

| ID | Description | Expected result |
|---|---|---|
| T-IO-01 | Newick string parsing | correct tip_labels/edges/lengths |
| T-IO-02 | Newick file reading | consistent with string parsing |
| T-IO-03 | NEXUS multi-tree file reading | returns all trees |
| T-IO-04 | Newick with square-bracket annotations | annotations stripped correctly |
| T-IO-05 | Newick with quoted labels | quoted labels parsed correctly |
| T-IO-06 | Newick round trip (parse → write → parse) | structure identical |
| T-IO-07 | long-format table to presence/absence matrix | correct sites×species DataFrame |
| T-IO-08 | wide-format table to presence/absence matrix | correct conversion |
| T-IO-09 | align_tree_matrix tree-matrix alignment | correct pruning/dropping + warning |
| T-IO-10 | tree and matrix share no species | ValueError |

### 3.4 Core algorithm tests (T-CORE)

| ID | Description | Expected result |
|---|---|---|
| T-CORE-01 | TreeArray.from_newick builds correctly | n_tips/n_edges/root_age correct |
| T-CORE-02 | node depth computation | each node's depth = parent depth + branch length |
| T-CORE-03 | tip_intervals correctness | each edge's [lo, hi) covers the correct tip range |
| T-CORE-04 | contribution matrix C[e,j] | equals max(0, min(d_e, t_j) - max(b_e, t_{j-1})) |
| T-CORE-05 | slicing window (equal-time) | windows are equally wide and cover [0, T] |
| T-CORE-06 | slicing window (equal-PD) | every slice has the same PD |
| T-CORE-07 | subarray pruning | keeps the MRCA subtree |
| T-CORE-08 | ultrametric detection | returns True for an ultrametric tree |
| T-CORE-09 | non-ultrametric tree detection | returns False + warning/error |
| T-CORE-10 | total_pd computation | equals the sum of all branch lengths |

### 3.5 Diversity index tests (T-INDICES)

| ID | Description | Expected result |
|---|---|---|
| T-INDICES-01 | pd() total PD computed correctly | matches a hand calculation |
| T-INDICES-02 | pe() phylogenetic endemism computed correctly | matches a hand calculation |
| T-INDICES-03 | pb() whole-tree beta diversity | value lies in [0, 1] |
| T-INDICES-04 | pd_per_slice PD per slice | row sums equal total PD |
| T-INDICES-05 | pe_per_slice PE per slice | row sums equal total PE |
| T-INDICES-06 | pb_per_slice row sums equal 1 (partition) | np.allclose(row_sums, 1.0) |
| T-INDICES-07 | pb_per_slice sorensen/turnover/nestedness | every component correct |
| T-INDICES-08 | pb_per_slice multisite/pairwise | both approaches correct |
| T-INDICES-09 | multisite_domain mixed vs paired | values as expected |
| T-INDICES-10 | degenerate neighbourhood handling | converged=False + reason |

### 3.6 Cumulative rate tests (T-RATES)

| ID | Description | Expected result |
|---|---|---|
| T-RATES-01 | cpd_rate returns the correct DataFrame | contains site_id/CpD/PD/pDO/converged columns, among others |
| T-RATES-02 | cpe_rate returns the correct DataFrame | contains CpE/PE/pEO |
| T-RATES-03 | cpb_rate returns the correct DataFrame | contains CpB/PB/pBO |
| T-RATES-04 | fit_rate recovers the parameters of a known index curve | r error < 1% |
| T-RATES-05 | origin_time computation | pXO = -ln(1-p)/r |
| T-RATES-06 | empty-site handling | converged=False, non-empty reason |
| T-RATES-07 | boundary-fit flags | at_lower_bound/at_upper_bound |
| T-RATES-08 | sensitivity slice-count sensitivity | returns suggested_slices |
| T-RATES-09 | rate_profile (normalize=True) | row sums equal 1 |
| T-RATES-10 | rate_profile (normalize=False) | row sums equal total PD |

### 3.7 Slicing tests (T-SLICING)

| ID | Description | Expected result |
|---|---|---|
| T-SLICING-01 | slice_rootward(time=0) → empty slice | PD=0 |
| T-SLICING-02 | slice_rootward(time=T) → whole tree | PD=total_pd |
| T-SLICING-03 | slice_tipward(time=0) → whole tree | PD=total_pd |
| T-SLICING-04 | slice_tipward(time=T) → empty slice | PD=0 |
| T-SLICING-05 | slice_interval(start, stop) | correct interval |
| T-SLICING-06 | slice_interval invert=True | returns an (older, younger) tuple |
| T-SLICING-07 | slice_pieces(n=10) → 10 slices | n_slices=10 |
| T-SLICING-08 | slice_pieces(width=0.5) | correct width |
| T-SLICING-09 | slice_pieces criterion="pd" | equal-PD slices |
| T-SLICING-10 | prune_tips prunes correctly | correct tips kept/dropped |
| T-SLICING-11 | time beyond root_age → ValueError | clear error message |
| T-SLICING-12 | negative time → ValueError | clear error message |

### 3.8 Posterior uncertainty tests (T-UNCERTAINTY)

| ID | Description | Expected result |
|---|---|---|
| T-UNC-01 | posterior_rate with a directory input | correct per-site summary |
| T-UNC-02 | posterior_rate with a Newick string input | same as above |
| T-UNC-03 | posterior_origin computation | contains origin_median/hdi_low/high |
| T-UNC-04 | bootstrap_sites | returns n_boot rows |
| T-UNC-05 | iter_trees with multi-format input | iterates over all trees correctly |
| T-UNC-06 | bad tree skipped + warning | no interruption, n_trees decreases |
| T-UNC-07 | HDI interval computation | hdi_low <= rate_median <= hdi_high |
| T-UNC-08 | ESS proxy computation | ess_proxy <= n_trees |

### 3.9 Spatial utility tests (T-SPATIAL)

| ID | Description | Expected result |
|---|---|---|
| T-SPATIAL-01 | adjacency_from_bounds queen | diagonal = 1, corner contact = 1 |
| T-SPATIAL-02 | adjacency_from_bounds rook | corner contact = 0 |
| T-SPATIAL-03 | adjacency_from_coords knn | symmetric, diagonal = 1 |
| T-SPATIAL-04 | adjacency_from_coords radius | symmetric |
| T-SPATIAL-05 | to/from_adjacency_dict round trip | consistent |
| T-SPATIAL-06 | rate_map produces a figure | returns a matplotlib Axes |

### 3.10 Compatibility layer tests (T-COMPAT)

| ID | Description | Expected result |
|---|---|---|
| T-COMPAT-01 | squeeze_root == slice_rootward | identical results |
| T-COMPAT-02 | squeeze_tips == slice_tipward | identical results |
| T-COMPAT-03 | squeeze_int invert return order | (younger, older), R style |
| T-COMPAT-04 | phylo_pieces timeSteps | monotonically increasing depth vector |
| T-COMPAT-05 | prune_tips method=1/2 | side after/before |

### 3.11 Simulation and dataset tests (T-SIM)

| ID | Description | Expected result |
|---|---|---|
| T-SIM-01 | yule_tree generates correctly | n_tips correct, ultrametric |
| T-SIM-02 | birth_death_tree generates correctly | n_tips correct, ultrametric |
| T-SIM-03 | make_dataset random | matrix with the correct dimensions |
| T-SIM-04 | make_dataset clustered | spatial autocorrelation |
| T-SIM-05 | sim_passerines dataset | contains trees/mat/coords/bounds/adj |
| T-SIM-06 | seed reproducibility | same seed, same result |

### 3.12 Visualization tests (T-VIZ)

| ID | Description | Expected result |
|---|---|---|
| T-VIZ-01 | plot_rate_line | returns an Axes |
| T-VIZ-02 | plot_sensitivity | returns an Axes |
| T-VIZ-03 | plot_posterior | returns an Axes |

### 3.13 End-to-end tests on realistic data (T-E2E)

| ID | Description | Expected result |
|---|---|---|
| T-E2E-01 | complete CpD workflow (API) | meaningful rate results |
| T-E2E-02 | complete CpE workflow (API) | meaningful rate results |
| T-E2E-03 | complete CpB workflow (API) | meaningful rate results |
| T-E2E-04 | complete posterior analysis workflow (API) | meaningful HDI intervals |
| T-E2E-05 | complete CLI workflow | all subcommands run correctly |
| T-E2E-06 | cross-validation (validate) | equivalence rate = 1.0 |
| T-E2E-07 | full visualization workflow | figure files generated |

### 3.14 Exception and boundary tests (T-ERROR)

| ID | Description | Expected result |
|---|---|---|
| T-ERR-01 | 1-tip tree → ValueError | explicit error message |
| T-ERR-02 | NaN branch length → ValueError | explicit error message |
| T-ERR-03 | negative branch length → ValueError | explicit error message |
| T-ERR-04 | duplicate tip labels → ValueError | explicit error message |
| T-ERR-05 | cycle → ValueError | explicit error message |
| T-ERR-06 | disconnected tree → ValueError | explicit error message |
| T-ERR-07 | n_slices < 1 → ValueError | explicit error message |
| T-ERR-08 | percentile outside (0,1) → ValueError | explicit error message |
| T-ERR-09 | hdi outside (0,1] → ValueError | explicit error message |
| T-ERR-10 | invalid criterion → ValueError | explicit error message |

### 3.15 Performance and benchmark tests (T-PERF)

| ID | Description | Expected result |
|---|---|---|
| T-PERF-01 | CpD rates on a 100-tip tree | < 5s |
| T-PERF-02 | CpD rates on a 1000-tip tree | < 60s |
| T-PERF-03 | slice-count scaling (10→200) | linear time |
| T-PERF-04 | reasonable memory usage | no memory leak |
| T-PERF-05 | comparison against the bundled benchmarks (`benchmarks/`) | results reproducible |

---

## 4. Test data

### 4.1 Simulated data

- **Small tree**: 4-tip symmetric tree `((A:1,B:1):1,(C:1,D:1):1);`
- **Medium tree**: 20-tip Yule tree (seed=42, age=10.0)
- **Large tree**: 100-tip birth-death tree (seed=42, age=10.0)
- **Posterior tree set**: 5 fifty-tip Yule trees (different seeds)
- **Presence/absence matrix**: 30 sites × 20 species (seed=7)
- **Adjacency matrix**: queen adjacency (6×6 grid, first 30 cells)

### 4.2 Existing validation reference data

Cross-validation uses the treesliceR reference files in `validation/reference/`.

---

## 5. Test execution flow

```
1. Generate test data (testing_suite/scripts/generate_test_data.py)
2. Run the unit tests (pytest, tests/)
3. Run the end-to-end tests (testing_suite/scripts/run_all_tests.py)
4. Run the benchmark tests (testing_suite/scripts/run_benchmarks.py)
5. Aggregate the results → test report (testing_suite/logs/test_results.json → testing_suite/TEST_REPORT.md)
```

---

## 6. Pass/fail criteria

- **Pass**: every assertion holds, no exception escapes, and the output files are correct.
- **Fail**: any assertion fails, an exception is uncaught, or output is missing/incorrect.
- **Overall**: a 100% pass rate is required before the work moves to the repository tidy-up stage.
