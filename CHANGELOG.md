# Changelog

All notable changes to PhyloSlicer are documented here. The version of record
is `phyloslicer.__version__`; the date below is the date that entry was tagged.

## 0.1.0 — 2026-09-29

Initial implementation.

### Added

- Slicing as a single `(E, k)` interval-intersection contribution matrix, so no
  tree is ever copied per slice (`slicing/`).
- Sparse membership products for per-site, per-slice PD, PE and PB(+RW)
  (`indices/`).
- Certified cumulative-rate fitting: bounded multi-start least squares with
  `converged`, `r_squared`, `residuals` and an explicit `reason` per row
  (`rates/`).
- `origin_within_tree`, which flags a fitted diversity origin older than the
  root, with `reason` recording `origin extrapolated beyond root` (`rates/`).
- Posterior-tree uncertainty: `posterior_rate`, `posterior_origin`,
  `bootstrap_sites` (`uncertainty/`).
- Spatial workflow helpers: rook/queen/knn/radius adjacency builders and rate
  maps accepting GeoDataFrames (`spatial/`, `viz/`).
- `treesliceR`-style compatibility layer (`compat/`) and bundled simulated
  datasets (`datasets/`).
- CLI: `phyloslicer slice | rates | posterior | validate` with JSON run
  summaries (`cli/`).
- `examples/quickstart.py`, a single-command out-of-the-box run on the bundled
  simulated dataset; executed by CI so the documented example cannot rot.
- Cross-validation against treesliceR v1.1.0 on 20/100-tip birth-death and
  coalescent scenarios, with a SHA-256 provenance manifest
  (`validation/`).
- The benchmark ladder: in-process naive baseline, cross-language slicing kernel
  at fixed `k`, scaling across tree sizes, and end-to-end CpD at 1000/2000 tips,
  plus `benchmarks/benchmark_memory.py` for Python-side peak memory and
  `benchmarks/recompute_ratios.py` to recompute every quoted speed-up range from
  the archived tables.
- Test suite (`tests/`) and end-to-end test plan/report (`testing_suite/`).
- `.pre-commit-config.yaml` (ruff check, ruff format) matching the `dev` extra.
