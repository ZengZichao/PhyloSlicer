# PhyloSlicer

**Temporal slicing of phylogenies and cumulative-diversity rate analysis for macroecology in Python.**

Version 0.1.0 · BSD-3-Clause · Python ≥ 3.10

Repository: <https://github.com/ZengZichao/PhyloSlicer>

---

## Overview

PhyloSlicer is an independent, vectorised Python implementation of the
"slice phylogenies to infer phylogenetic patterns over evolutionary time"
paradigm introduced by **treesliceR** (Araujo et al. 2025, *Ecography* 2025,
e07364, <https://doi.org/10.1111/ecog.07364>), extended with posterior-tree
uncertainty quantification, robust fitting diagnostics, spatial workflow
helpers and a command-line interface.

Speed is a structural consequence, not a marketing word: slicing is a single
`E × k` interval-intersection contribution matrix instead of the `k` complete
tree copies the reference design materialises. The one array still grows with `k`,
but the `k` whole-tree copies do not. Python-side peaks are measured
(`benchmarks/benchmark_memory.py` → `benchmarks/results_memory_py.csv`): at 2000 tips the
single `(E, k)` array peaks at 1.93 / 3.49 / 6.69 MB for k = 25 / 50 / 100, against
2.81 / 5.62 / 11.23 MB for `k` array-level tree copies — a 1.45× to 1.68× saving that
widens with `k`. Because this package's own matrix also grows with `k`, the honest
claim is about **object count and copy structure**, not about a memory ratio that
diverges; the R-side peak of treesliceR is **not** measured, and the documentation
keeps to the same qualification.

### Key features

1. **Array-based slicing** — slicing is expressed as an interval-intersection
   contribution matrix (`C[e, j] = max(0, min(d_e, t_j) − max(b_e, t_{j−1}))`),
   so no tree is ever copied.
2. **Sparse-matrix indices** — per-site, per-slice PD/PE/PB(+RW) via a single
   sparse product `A @ C`.
3. **Robust rate fitting** — bounded multi-start optimisation with `converged`,
   `r_squared`, `residuals` and an explicit `reason` for every failure.
4. **Posterior-tree uncertainty** — `posterior_rate` / `posterior_origin`
   credible intervals and `bootstrap_sites` stability analysis.
5. **Spatial workflow** — rook/queen/knn/radius adjacency builders and rate
   maps that accept GeoDataFrames (geopandas optional).
6. **CLI** — `phyloslicer slice | rates | posterior | validate` with JSON run
   summaries.

## Documentation

- **English user guide (primary):** [docs/USER_GUIDE.md](docs/USER_GUIDE.md)
- **Chinese user guide:** [docs/USER_GUIDE_CN.md](docs/USER_GUIDE_CN.md)
- **Changelog:** [CHANGELOG.md](CHANGELOG.md) · Chinese: [CHANGELOG_CN.md](CHANGELOG_CN.md)
- **Contributing:** [CONTRIBUTING.md](CONTRIBUTING.md) · Chinese: [CONTRIBUTING_CN.md](CONTRIBUTING_CN.md)
- **Release procedure:** [RELEASE.md](RELEASE.md) · Chinese: [RELEASE_CN.md](RELEASE_CN.md)

## Installation

```bash
pip install -e .            # core (numpy, scipy, pandas, matplotlib)
pip install -e ".[io]"      # + dendropy/biopython adapters
pip install -e ".[spatial]" # + geopandas
pip install -e ".[accel]"   # + numba batched kernels
pip install -e ".[dev]"     # + pytest, ruff, mypy, pre-commit
```

## Quick start

Run the bundled example first — it needs no data, no download and no setup, and
it is the workflow the user guide walks through:

```bash
python examples/quickstart.py
```

The snippet below is the same workflow on **your own** files:

```python
import phyloslicer as ps

tree = ps.TreeArray.from_newick("birds.tre")
mat = ps.io.presence_matrix(long_table, site_col="grid", species_col="species")
adj = ps.spatial.adjacency_from_coords(x, y, method="knn", k=8)

rates = ps.cpd_rate(tree, mat, n_slices=100)      # CpD, PD, pDO + diagnostics
post = ps.posterior_rate("posterior_trees/", mat, n_slices=100)
beta = ps.cpb_rate(tree, mat, adj, component="sorensen")
```

Command line:

```bash
phyloslicer rates input.tre sites.csv --rate cpd --n-slices 100 --out rates.csv
phyloslicer slice input.tre --rootward 5.0 --out crown.tre
```

## Testing and validation

```bash
python -m pytest tests/                     # unit + regression tests
python validation/compare.py                # Python vs treesliceR equivalence
```

The comparison set is `{PD, PE, CpD, CpE}` plus the beta families (`CpB`,
`CpB_RW`, `r_phylo(index = "PB")`) — 11 of the 14 beta component-by-approach
families are cross-checked numerically against treesliceR 1.1.0. See
[`validation/NOTES.md`](validation/NOTES.md) for the full gap list and
tolerance rationale.

## Benchmarks

`benchmarks/` holds four independent protocols and the archived tables they
produced. Read [`benchmarks/README.md`](benchmarks/README.md) before quoting a
number from them: two column names mean different things in different files, and
only the fixed-k scaling tables are cross-language. `python
benchmarks/recompute_ratios.py` recomputes every documented speed-up range from
the CSVs and is run by CI.

## Citation

If you use PhyloSlicer, please cite both this package and the methodological
source. The package citation is in [`CITATION.cff`](CITATION.cff); the
version-agnostic Zenodo concept DOI for all releases is
[10.5281/zenodo.23057620](https://doi.org/10.5281/zenodo.23057620).

> Araujo, M.L., Ferreira, L.G.S.S., Nakamura, G., Coelho, M.T.P., Rangel, T.F.
> (2025) 'treesliceR': a package for slicing phylogenies and inferring
> phylogenetic patterns over evolutionary time. *Ecography* 2025, e07364.
> https://doi.org/10.1111/ecog.07364

## License

Code: BSD-3-Clause (see [LICENSE](LICENSE)). Documentation: CC-BY 4.0.
Third-party attribution, including the origin of the `treesliceR` reference
values: [NOTICE.md](NOTICE.md) · Chinese: [NOTICE_CN.md](NOTICE_CN.md).

## Author

Zichao Zeng (曾子超) — School of Life Sciences and Biotechnology, Shanghai
Jiao Tong University · <zengzichao@sjtu.edu.cn> ·
ORCID [0000-0001-6553-970X](https://orcid.org/0000-0001-6553-970X)
