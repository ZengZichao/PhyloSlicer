# Changelog

All notable changes to PhyloSlicer are documented here. The version of record
is `phyloslicer.__version__`; the date below is the date that version was made
public, and an entry records an earlier source-finalisation date in its own body.

## 0.1.1 — 2026-10-18

First release on PyPI.  The 0.1.0 lineage was hardened per the 2026-10-08 review
report: the one real library defect was fixed, the verification suite was widened
to pin the corrected behaviour, and the release pipeline gained PyPI publication.
The source was finalised and made public on 2026-10-18.

### Added

- PyPI publication through GitHub Actions Trusted Publishing
  (`.github/workflows/publish.yml`): pushing a `v*` tag builds the sdist and
  wheel, checks the metadata and uploads both with OIDC, so the repository
  holds no API token.
- Python 3.14 in the CI matrix and the classifiers.
- Regression coverage that pins the corrected behaviour: the posterior audit
  columns against mixed hypometric pools, the published numbers, the output
  schema, a CI mirror, and the review archive; the β-partition invariants
  cover all 24 slice-and-habitat combinations instead of 6.
- Review-archive guards: DOI shape and AST-fragment checks, an explicit
  `--allow-no-git` / `--allow-dirty` with dirty paths recorded, and identity
  literals split out to `tools/identity_literals.py`.

### Fixed

- `posterior_rate` / `posterior_origin` computed their audit columns on the
  wrong axis (`np.vstack` flattened the draws before indexing by site); they
  now stack to `(n_draws, n_sites, 2)` and read the columns from there.
- `validation/compare.py` refuses to overwrite the version-controlled
  `equivalence_report.csv` by default (add `--in-place`), and the report gains
  an `n_finite` column so the quoted "35 950 / 36 000" is re-derivable from the
  delivered table.
- CI maintenance hardening: double-run guard, CodeQL, build verification and
  protection of the pinned ruff version; Dependabot covers the GitHub Actions
  ecosystem with a major-version ignore rule.
- `LICENSE` credits the second author; `NOTICE` and the READMEs state that the
  English documents are authoritative and the Chinese ones are structural
  translations of the same structure.

## 0.1.0 — 2026-09-30

Initial implementation.  The 0.1.0 source was finalised on 2026-09-29; the public
GitHub release and the Zenodo deposit were both made on 2026-09-30, which is the
release date of the version and the date the Software Notes rules ask for.

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
- A release-documentation set beyond the manuals: `docs/CI_WORKFLOWS*.md` (a
  verbatim mirror of the CI workflow definitions, since `.github/` and `tools/`
  cannot ship in a review archive), `docs/WORD_COUNT*.md` (the counting rule that
  makes the manuscript's length statistics re-derivable) and
  `docs/ZENODO_DEPOSIT_CHECKLIST*.md` (the field-by-field correction of the
  deposited record). Every English document ships a Chinese counterpart, which
  `tests/test_release_hygiene.py` enforces.
