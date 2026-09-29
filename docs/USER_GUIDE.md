# PhyloSlicer — User Guide

**Temporal slicing of phylogenies and cumulative-diversity rate analysis for macroecology in Python.**

Speed here is a structural consequence, not a marketing word: slicing is a single
`E x k` interval-intersection contribution matrix instead of the `k` complete tree
copies the reference design materialises, so neither runtime nor memory is multiplied
by the number of slices. That is the claim this sentence supports; the measured
evidence lives in `benchmarks/`, where the cross-language speed-ups are read off
the fixed-k scaling tables (`results_scaling_py.csv` vs `results_scaling_r.csv`,
slicing kernel only, byte-identical input trees, median of 15 timed repetitions
per point) while `results.csv` compares this package against *its own* naive
implementation and is therefore *not* a cross-language benchmark. No timing figure
is quoted here: the R columns are the archived benchmark outputs, and R timings
on a laptop vary between sessions - the reference-package notes warn about a
factor of two. The three R scripts were re-run here under R 4.5.3 (aarch64-apple-darwin20.0.0) with ape 5.8.1 and treesliceR 1.1.0 into
scratch files and reproduced the archived medians to within 6 % (median ratio
1.06 on the validation scenarios, 1.02-1.08 on the 1000/2000-tip ones), but the
committed CSVs were left in place so that the numbers the documentation quotes
are the ones they were measured against.

Version 0.1.0  ·  BSD-3-Clause  ·  Python ≥ 3.10

PhyloSlicer is an independent, vectorised Python implementation of the
"slice phylogenies to infer phylogenetic patterns over evolutionary time"
paradigm introduced by **treesliceR** (Araujo et al. 2025, *Ecography* 2025,
e07364, doi:10.1111/ecog.07364; first published online 28 October 2024),
extended with posterior-tree uncertainty quantification, robust fitting
diagnostics, spatial workflow helpers and a command-line interface.
Chinese documentation: [用户手册（中文版）](USER_GUIDE_CN.md).

---

## Table of contents

1. [Installation](#1-installation)
2. [Five-minute quick start](#2-five-minute-quick-start)
3. [Core concepts and conventions](#3-core-concepts-and-conventions)
4. [Input/output: trees and matrices](#4-inputoutput-trees-and-matrices)
5. [Slicing](#5-slicing)
6. [Diversity indices](#6-diversity-indices-pd-pe-pb)
7. [Cumulative rates](#7-cumulative-rates-cpd-cpe-cpb)
8. [Uncertainty across posterior trees](#8-uncertainty-across-posterior-trees)
9. [Spatial helpers](#9-spatial-helpers)
10. [Visualisation](#10-visualisation)
11. [Simulation and bundled data](#11-simulation-and-bundled-data)
12. [Command-line interface](#12-command-line-interface)
13. [Moving from treesliceR (R)](#13-moving-from-treeslicerr)
14. [Diagnostics, statuses and error messages](#14-diagnostics-statuses-and-error-messages)
15. [Numerical validation](#15-numerical-validation)
16. [FAQ and troubleshooting](#16-faq-and-troubleshooting)
17. [API quick reference](#17-api-quick-reference)

---

## 1. Installation

```bash
# core package (numpy, scipy, pandas, matplotlib)
pip install -e .

# optional extras
pip install -e ".[io]"       # + DendroPy / Biopython adapters
pip install -e ".[spatial]"  # + GeoPandas
pip install -e ".[accel]"    # + numba (batched fitting kernels)
pip install -e ".[dev]"      # + pytest, ruff, mypy, ...
```

Everything works without the extras; missing optional dependencies produce
clear `ImportError` messages (or, for spatial grids, an automatic fallback
to rectangle bounds).

## 2. Five-minute quick start

```python
import phyloslicer as ps

# 1. read a dated (ultrametric) Newick tree
tree = ps.TreeArray.from_newick("birds.tre")

# 2. sites x species presence/absence matrix (DataFrame; see section 4)
mat = ps.io.presence_matrix(long_table, site_col="grid_cell", species_col="species")

# 3. slice the phylogeny into 100 equal-time slices and compute CpD rates
rates = ps.cpd_rate(tree, mat, n_slices=100)
print(rates.head())
#    site_id       CpD        PD       pDO  converged  r_squared  reason
# 0      s_1  1.234...  45.6...   2.42...       True     0.98...

# 4. phylogenetic beta-diversity rates over queen-adjacent neighbours
adj = ps.spatial.adjacency_from_coords(x, y, method="knn", k=8)
beta = ps.cpb_rate(tree, mat, adj, n_slices=100, component="sorensen")

# 5. uncertainty across a folder of posterior trees
post = ps.posterior_rate("posterior_trees/", mat, n_slices=100)
```

Command line equivalents:

```bash
phyloslicer rates   birds.tre sites.csv --rate cpd --n-slices 100 --out rates.csv
phyloslicer slice   birds.tre --rootward 5.0 --out crown.tre
phyloslicer posterior posterior_trees/ sites.csv --hdi 0.95 --out posterior.csv
phyloslicer validate birds.tre sites.csv --reference r_reference.csv --out report.csv
```

## 3. Core concepts and conventions

### 3.1 Time, depth and ages

PhyloSlicer assumes a **dated, ultrametric** tree. Two coordinate systems
are used internally:

* **depth from the root** (0 at the root, `root_age` at the tips) — the
  natural coordinate for interval arithmetic;
* **age before present** (`root_age − depth`) — the coordinate users see
  (Ma, or any branch-length time unit your tree uses).

Slice windows are always reported as ages, ordered root → tip, and the
`ages` vector used by the rate model is the **older boundary** of each
slice, expressed in **Ma before present** (`root_age − depth`), not in
depth-from-root. Non-ultrametric trees are refused by default:
`slice_pieces()` raises `ValueError` unless you pass
`allow_non_ultrametric=True`, which downgrades the check to a warning and
records `ultrametric = False` in the result table; no rescaling of the tree
is performed, and `root_age` is then the deepest tip depth.

> Why the default is an error: a non-ultrametric tree still produces numbers,
> but they no longer mean "diversity before time *t*", and a silent warning is
> easy to lose in log output. Opt in explicitly with
> `allow_non_ultrametric=True` when you want the audit columns instead.

### 3.2 The contribution matrix

Every edge `e` spans the interval `[b_e, d_e]` in depth coordinates, and a
slice is a window `[t_{j−1}, t_j]`. The branch length edge `e` contributes
to slice `j` is

```
C[e, j] = max(0, min(d_e, t_j) − max(b_e, t_{j−1}))
```

computed for all edges × slices in one vectorised operation. **No tree is
ever copied**: a stack of k slices is one `E × k` array, not k tree objects
(this is the structural difference from treesliceR's `phylo_pieces`, which
materialises one modified copy of the whole phylogeny per slice).

### 3.3 Presence/absence matrices

All index and rate functions take a **sites × species** `pandas.DataFrame`:
rows are assemblages (the `index` becomes `site_id` in every output),
columns are species names. Values are thresholded as `> 0`. Species names
are matched against tree tip labels; `align_tree_matrix()` (applied
automatically by the rate functions) prunes tree tips absent from the
matrix, drops matrix species absent from the tree and drops never-observed
species — **with a warning for every adjustment**.

### 3.4 Rate model

Per assemblage, per-slice index values are accumulated root → tip,
normalised to sum to 1, and fitted to `cumulative ∝ exp(−r · age)` by
bounded least squares (bounds `[1e-6, 50]`). Every fit returns:

| field | meaning |
|---|---|
| `rate` | fitted r (NaN when the fit failed) |
| `converged` | True only when the optimiser certified the optimum |
| `r_squared` | 1 − SSE/SST of the fitted cumulative curve |
| `residuals` | per-slice residuals (empty when `diagnostics=False`) |
| `start_used`, `n_iter` | optimiser provenance |
| `reason` | empty when converged; a plain-language explanation otherwise |

Origin times follow `pXO = −ln(1 − percentile) / r`; the default
`percentile=0.95` matches treesliceR's default `pDO = 5`: both return
`-ln(0.05)/r`, the age at which the fitted curve is 0.05 — the reference names the
remaining 5 % directly, this package names the accumulated 95 % complement.

## 4. Input/output: trees and matrices

### 4.1 Trees

```python
from phyloslicer import TreeArray
from phyloslicer.io import read_newick, read_nexus, parse_newick, write_newick, write_nexus

tree  = TreeArray.from_newick("tree.tre")     # file path or Newick string
tree  = read_newick("(A:1.0,B:1.0):1.0;")     # string directly
trees = TreeArray.from_nexus("posterior.nex") # list[TreeArray], multi-tree
tree  = parse_newick("((A:1,B:1)95:1,(C:1,D:1)100:1);")  # internal labels read but not stored

tree.to_newick()                # Newick string (root has no stem length)
write_nexus(trees)              # NEXUS string with all trees
tree.to_dendropy()              # requires phyloslicer[io]
tree.to_biopython()
TreeArray.from_dendropy(dnd_tree)
TreeArray.from_biopython(bp_tree)
```

Notes and guarantees:

* Comments in square brackets (including NHX annotations) are stripped
  before parsing; quoted names (`'sp one'`) are supported.
* A `TREE` statement inside a NEXUS `TREES` block **may span several
  lines**.
* Tip labels are renumbered in depth-first order at construction; every
  clade's tips form a contiguous index range, which is what makes the
  interval kernels possible.
* Internal (support) labels are **not preserved** on round-trip —
  PhyloSlicer's data model stores only tips, edges and lengths. Keep the
  original file if you need support values.
* Use `read_newick` for the first tree of a multi-tree Newick file; use
  `iter_trees()` (section 8) to stream all of them.

Key `TreeArray` members:

```python
tree.n_tips, tree.n_edges, tree.root, tree.tip_labels
tree.root_age          # deepest tip depth (== tip depth for ultrametric trees)
tree.tip_ages()        # root_age - depth, per tip
tree.is_ultrametric    # tolerance-based check
tree.check_ultrametric()  # warns when violated
tree.total_pd()        # sum of branch lengths
tree.node_depths()     # (n_nodes,) depth from root
tree.tip_intervals()   # (E, 2) [lo, hi) tip-index ranges covered by each edge
tree.subarray(["A", "B"])  # prune to a subset (MRCA-rooted)
tree.with_lengths(v)   # SlicedTree sharing the topology with new lengths
```

### 4.2 Presence/absence matrices

```python
from phyloslicer.io import presence_matrix, align_tree_matrix, as_matrix

# wide table: sites as rows, species as columns (light validation only)
wide = presence_matrix(wide_df)

# long table: one row per (site, species[, value]) observation
wide = presence_matrix(long_df,
                       site_col="grid", species_col="species",
                       value_col="count")       # value > 0  -> presence
wide = presence_matrix(long_df, site_index_col=True)  # site ids from the index

# manual alignment with explicit reporting (rate functions do this for you)
tree2, mat2 = align_tree_matrix(tree, mat)   # both realigned + warnings
mat_df = as_matrix(numpy_array_or_sparse)    # coerce to DataFrame
```

`align_tree_matrix` performs, with a warning each: removal of all-zero
species columns, removal of matrix species absent from the tree, pruning of
tree tips absent from the matrix, and reordering of columns to the (pruned)
tree's tip order.

## 5. Slicing

### 5.1 Single slices

```python
from phyloslicer.slicing import slice_rootward, slice_tipward, slice_interval

crown = slice_rootward(tree, 5.0)               # keep the youngest 5 Ma
back  = slice_tipward(tree, 5.0)                # drop the youngest 5 Ma
mid   = slice_interval(tree, start=8.0, stop=5.0)   # the [5, 8] Ma band
older, younger = slice_interval(tree, 8.0, 5.0, invert=True)

pd_slice = slice_rootward(tree, 2.0, criterion="pd")  # remove root-side PD of 2
```

Semantics (the **endpoint behaviour** mirrors treesliceR; the placement of
equal-PD windows diverges on purpose — see §13.1):

| call | keeps | `criterion="time"` | `criterion="pd"` |
|---|---|---|---|
| `slice_rootward(tree, x)` | crown slice | the youngest `x` Ma | tree minus root-side PD `x` |
| `slice_tipward(tree, x)` | trunk | everything except the youngest `x` Ma | root-side PD `x` |
| `slice_interval(tree, a, b)` | the `[b, a]` band (needs `a > b` in Ma) | band between ages | band between root-side PD depths (`a < b`) |

Edge cases are validated: `time < 0` raises; `time > root_age` raises;
`slice_rootward(tree, 0)` returns an **empty** slice (PD = 0) and
`slice_tipward(tree, 0)` the full tree — this endpoint behaviour matches the
reference. Thresholds measured in PD run through the cumulative-PD curve of the
tree (`depth_for_pd`), which counts the branches **actually active** at each
depth. On a strictly binary tree with distinct node depths that reproduces the
reference cut points; on a multifurcating tree it deliberately does not, and the
difference is unvalidated (§13.1).

The result is a `SlicedTree`: same topology as the parent, truncated branch
lengths, zero-length edges retained. Use

```python
sliced.pd              # PD stored in the slice
sliced.tree            # a TreeArray copy
sliced.collapse()      # TreeArray with zero-length internal edges removed
```

All slicing functions accept `collapse=True` to return the collapsed
`TreeArray` directly.

### 5.2 Many slices — `SliceStack`

```python
from phyloslicer.slicing import slice_pieces

stack = slice_pieces(tree, n=100)                     # 100 equal-time slices
stack = slice_pieces(tree, width=0.5)                 # fixed 0.5-Ma width (rounded
                                                      # treesliceR-style to fit)
stack = slice_pieces(tree, n=50, criterion="pd")      # equal-PD slices
```

`n` and `width` are mutually exclusive (this replaces R's easily-confused
`n` + `method` pair). The `SliceStack` gives you:

```python
stack.n_slices
stack.windows        # (k, 2) ages [older, younger], root -> tip
stack.ages           # (k,) older boundary of each slice  <- rate-model ages
stack.widths         # (k,) slice widths
stack.contribution()             # (k,) branch length per slice (all tips)
stack.contribution(["sp_1"])     # idem restricted to a tip subset
stack.contribution_matrix()      # (E, k) full edge x slice matrix
stack.to_sliced_trees()          # materialise slices as SlicedTree objects
```

### 5.3 Pruning tips

```python
from phyloslicer.slicing import prune_tips

kept = prune_tips(tree, 0.4)                    # keep tips with terminal branch >= 0.4
kept = prune_tips(tree, 0.4, side="before")     # drop them instead (R method=2)
kept = prune_tips(tree, 0.75, quantiles=True)   # threshold = 75th percentile
lists = prune_tips(tree, np.array([0.2, 0.4]))  # one tree per threshold
```

A threshold that would leave fewer than two tips raises `ValueError` naming the
threshold and the workable range. treesliceR returns `NULL` there, which in a
script surfaces later as an attribute error on `None`; one element of a vector
call cannot now quietly yield a missing list entry.

## 6. Diversity indices (PD, PE, PB)

> **Naming note:** `from phyloslicer.indices import pd` binds the *PD
> function* to the name `pd`, shadowing the conventional `import pandas as
> pd` alias. In scripts that use pandas, prefer
> `from phyloslicer import indices` and call `indices.pd(...)`, or import
> the function under a different alias.

```python
from phyloslicer.indices import pd, pe, pb, pd_per_slice, pe_per_slice, pb_per_slice

faith  = pd(tree, mat)        # (n_sites,) total PD per site
endo   = pe(tree, mat)        # (n_sites,) phylogenetic endemism (range-weighted)
beta   = pb(tree, mat, adj)   # DataFrame: site_id, PB   (full-tree beta)

stack = slice_pieces(tree, n=100)
prof_pd = pd_per_slice(tree, mat, stack)   # (n_sites, k)
prof_pe = pe_per_slice(tree, mat, stack)   # (n_sites, k)
prof_pb = pb_per_slice(tree, mat, adj, stack,
                       component="sorensen",     # "sorensen"|"turnover"|"nestedness"
                       approach="multisite",     # "multisite"|"pairwise"
                       weighted=False,           # True -> PB_RW (Laffan et al. 2018)
                       multisite_domain="mixed") # "paired" removes the drift, §6.1
# prof_pb["values"]  (n_focal, k) per-slice values, summing to the total
# prof_pb["total"]   (n_focal,) full-tree beta of the neighbourhood
# prof_pb["status"]  per-focal status string (section 14)
```

Implementation notes:

* Membership of sites in edges is built from the tip-interval table with
  one binary search per site — no `(sites × tips)` intermediate.
* All per-site, per-slice quantities are one sparse matrix product `A @ C`.
* Range weights for endemism are `1 / (number of sites covering the edge)`,
  matching treesliceR and Laffan et al. (2018).
* Every beta series **partitions the neighbourhood's total beta**: for every
  finite row, the per-slice values sum to 1. That is implemented behaviour, not
  an assumption - `pb_per_slice` divides each per-slice series by the sum of its
  *own* shares, and the returned `normalisation` column records which route a row
  took: `"full_tree"` when the sliced shares already add up to the whole-tree
  index (sorensen, whose `b + c` numerator is additive over slices, and any
  neighbourhood where the two denominators happen to coincide), `"profile"` when
  the denominator had to be rebuilt from `sum_j series_j`.
  The property is *newly* true for turnover and nestedness: `min(b, c)` is not
  additive across slices (`sum_j min(b_j, c_j) <= min(sum_j b_j, sum_j c_j)`), so
  dividing by a whole-tree denominator - what the code did before the A2 fix, and
  what treesliceR still does - leaves turnover rows summing to 0.96-0.98 and
  nestedness rows to 1.05-1.40 on the workload below.
  The assertion behind this bullet is
  `tests/test_regressions.py::TestBetaNormalisation::test_all_beta_profiles_partition_the_neighbourhood_beta`:
  the `bd20_random` fixture (the committed cross-validation inputs
  `bd20__tree.tre` + `bd20__random__mat.csv`, 20 tips / 30 sites, queen-grid
  adjacency) sliced with `slice_pieces(tree, n=30, criterion="time")`, all six
  component x approach combinations, checked with
  `np.allclose(row_sums, 1.0, atol=1e-8)` and a stricter
  `max |row_sum - 1| < 1e-12`. It replaced a 6-site / 5-slice fixture that passed
  at `atol=0.05`. Two companions keep the claim honest:
  `test_profile_and_total_share_the_same_quantities` pins `total` against `pb()`
  to 1e-12 (the profile is normalised, the reported index is not rescaled), and
  `test_r_compatible_profile_keeps_the_unnormalised_curve` pins
  `indices.r_compatible_profile` to treesliceR's *non*-partitioning curve.


### 6.1 The multisite index depends on neighbourhood size (`multisite_domain`)

treesliceR's multisite β mixes two summation domains (this is inherited from
the reference, not a porting error):

```
a_tot  = Σ_(m sites) PD_i − PD_union        grows ~ linearly in m
bc_tot = Σ_(C(m,2) pairs) (b + c)_pair      grows ~ quadratically in m
β      = bc_tot / (2·a_tot + bc_tot)
```

The shared term `a` is summed over the *m* sites while the turnover terms
`b + c` (and `min(b, c)`) are summed over the *C(m, 2)* pairs, so the
fraction is **not an invariant of the composition**: it rises monotonically
with neighbourhood size even when every individual pair is unchanged. On the
reference construction — one shared skeleton branch plus one private terminal
branch per site, all of length 1, neighbourhood = all *m* sites — the Sørensen
value is analytically `m / (m + 2)`, which the default `multisite_domain="mixed"`
reproduces exactly:

| m | `mixed` (default, = treesliceR) | analytic `m/(m+2)` | every pairwise comparison | `paired` |
|---|---|---|---|---|
| 2 | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| 3 | 0.6000 | 0.6000 | 0.5000 | 0.5000 |
| 5 | 0.7143 | 0.7143 | 0.5000 | 0.5000 |
| 10 | 0.8333 | 0.8333 | 0.5000 | 0.5000 |
| 20 | 0.9091 | 0.9091 | 0.5000 | 0.5000 |

**Therefore: never compare multisite values computed on neighbourhoods of
different size** (different `k` in a k-NN adjacency, or a grid whose focal
cells have different neighbour counts). Read them as a *pair-summed γ
decomposition*, not as a size-independent β.

* `multisite_domain="mixed"` (default; `"reference"` is an accepted synonym,
  see `indices.MULTISITE_DOMAIN_ALIASES`) reproduces the reference arithmetic.
  It stays the default because switching it would silently redefine a published
  metric mid-release.
  Note the scope of that caveat: it is about the **index** (`PB`), not about
  every downstream number. `pb_per_slice()` divides each row by its own sum, so
  for `sorensen` and `turnover` the multisite profile — and therefore the
  `CpB` fitted by `cpb_rate()` — is *identical* under both domains; only `PB`
  moves. `nestedness` is the one component whose rate changes with the domain,
  because its two terms carry different denominators.
* `multisite_domain="paired"` puts the shared term in the same pair-summed
  domain as `b + c` and `min(b, c)`, which removes the drift: the value becomes
  0.5000 at every *m* above and coincides with the `approach="pairwise"`
  aggregation of the same neighbourhood, while leaving `m = 2` untouched (the
  two domains coincide there). This is a **deliberate departure from
  treesliceR**: with `"paired"` the multisite output is no longer
  cross-comparable with R.

Note that following the reference arithmetic is not the same as satisfying a
published definition. The method paper behind treesliceR (Araujo et al. 2025,
*Ecography*, e07364, doi:10.1111/ecog.07364, as cited by this package) is the
only such entry this project carries, but the reference package's own
vignette still says "Araujo et al. (in review)" and does not print the
multisite formulas, so we could not verify `a`-by-site / `b+c`-by-pair against
a published equation. The internal inconsistency above is stated as such and
left opt-in until that check is done. If a published multiple-site index
(Baselga 2010 / Carrasco et al. 2013 family) is required, use
`approach="pairwise"` and aggregate explicitly.

## 7. Cumulative rates (CpD, CpE, CpB)

```python
import phyloslicer as ps

res = ps.cpd_rate(tree, mat, n_slices=100)                 # CpD + PD + pDO
res = ps.cpe_rate(tree, mat, n_slices=100)                 # CpE + PE + pEO
res = ps.cpb_rate(tree, mat, adj,
                  component="sorensen",                    # or turnover/nestedness
                  approach="multisite",                    # or pairwise
                  weighted=False,                          # True -> CpB_RW
                  n_slices=100)                            # CpB + PB + pBO
```

Each returns a DataFrame, one row per assemblage:

```
site_id, CpD, PD, pDO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
site_id, CpE, PE, pEO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
site_id, CpB (CpB_RW), PB (PB_RW), pBO, converged, r_squared, reason, at_lower_bound, at_upper_bound, ultrametric, origin_within_tree
```

The last four columns are the fit and slice audit trail: `at_lower_bound` /
`at_upper_bound` say the reported rate is a boundary optimum rather than an
estimate, `ultrametric` records whether the sliced tree passed the
ultrametricity check (`False` is only reachable via
`allow_non_ultrametric=True`), and `origin_within_tree` says whether the fitted
origin time lies inside the tree. `pXO = -ln(1 - percentile) / rate` is an
extrapolation of the fitted exponential, so a slow rate can place the origin
**before the root**, where the phylogeny contains no corresponding node; the
value is reported unchanged and the flag plus `reason` say so. `posterior_rate`
and `posterior_origin` carry `n_ultrametric_draws` and `origin_within_tree`
across the whole pool.

Common parameters:

* `criterion="time"` (default) or `"pd"` — how slices are cut.
* `percentile=0.95` — the origin time is the moment the fitted curve stands at
  `1 - percentile`, i.e. when only 5 % of the diversity had accumulated (a near-root
  moment); it is **not** the time by which 95 % had accumulated
  (`pXO = −ln(1 − percentile)/r`; `0.95` ≡ R's `pDO = 5`).
* `diagnostics=True` — keep residuals in the `RateFit` objects (set to
  `False` to save memory for very large runs).
* `multistart=8` — starts for the robust per-row refit fallback.
* parallelism is internal to the vectorised kernel; the single sparse product
  `A @ C` already solves all assemblages simultaneously.

Without fitting:

```python
prof = ps.rate_profile(tree, mat, n_slices=100, index="PD")   # "PD"|"PE"|"PB"|"PB_RW"
prof["values"]   # (n_sites, k) relative per-slice values (PB: beta shares)
prof["ages"]     # (k,)
prof["status"]   # per-site status ("ok", "empty_site", ... )
```

Slice-number sensitivity (mirrors `CpR_sensitivity`, plus a heuristic
suggestion):

```python
sens = ps.sensitivity(tree, mat, slice_grid=[10, 25, 50, 100, 200],
                      rate="cpd",              # "cpe"|"cpb"|"cpb_rw" (need adj)
                      sample_sites=0,          # 0 = use every site (R requires > 0)
                      seed=42)                 # deterministic site subsampling
# columns: site_id, n_slices, rate, converged, at_lower_bound,
#          at_upper_bound, suggested_slices
```

`suggested_slices` is the smallest grid count after which a site's rate
changes by less than 5% between consecutive grid points.

### 7.1 Fitting functions directly

```python
from phyloslicer.rates import fit_rate, origin_time

fit = fit_rate(profile, stack.ages)            # profile: (k,) root -> tip
fit.rate, fit.converged, fit.r_squared
fit.residuals, fit.start_used, fit.reason
origin_time(fit.rate, percentile=0.95)         # pXO in Ma
```

`fit_rate` runs a log-linear coarse estimate, a geometric multi-start grid
and a bounded least-squares polish; one failing start never aborts the
call. `method="mle"` switches to a Gaussian maximum-likelihood objective
with Nelder–Mead. The batched vectorised kernel `fit_rates_batch(profiles,
ages)` fits an `(n_sites, k)` matrix at once and is what `cpd_rate` /
`cpe_rate` use internally.

## 8. Uncertainty across posterior trees

```python
from phyloslicer.uncertainty import posterior_rate, posterior_origin, bootstrap_sites

# sources: a list of TreeArray, a directory of .tre files, or a
# multi-tree Newick/NEXUS file — trees are streamed, not all loaded
post = posterior_rate("posterior_trees/", mat, rate="cpd", n_slices=100, hdi=0.95)
```

`posterior_rate` returns per-site summaries:

```
site_id, rate_mean, rate_median, hdi_low, hdi_high, ess_proxy,
ess_proxy_caveat, n_trees, n_converged
```

* Trees are processed one at a time. A draw that cannot be *read* is skipped
  with a warning naming the file and the record (`<file>: skipping record j of
  N: ...`), and a draw that reads but cannot be analysed warns as
  `skipping tree i: ...`; if any draw was dropped at all, a closing warning
  states how many, because the summaries then rest on a smaller pool than the
  one on disk. `iter_trees` also accepts a multi-tree Newick **string**, a file
  or a directory.
* `hdi_low/high` are the empirical highest-density interval of the
  per-site rate across trees.
* `ess_proxy` is an effective-sample-size proxy from lagged
  autocorrelation; values far below `n_trees` indicate ordered (possibly
  autocorrelated) draws.
* Each tree is aligned to the matrix independently, so posterior pools
  with slightly different tip sets are handled correctly.

`posterior_origin` returns `site_id, origin_median, hdi_low, hdi_high` for
the origin time pXO.

```python
boot = bootstrap_sites(tree, mat, n_boot=1000, rate="cpd", n_slices=100, seed=0)
# DataFrame: replicate, rate   (rate refit on the pooled assemblage of each
#                              bootstrap resample of sites)
```

## 9. Spatial helpers

```python
from phyloslicer.spatial import (adjacency_from_grid, adjacency_from_bounds,
                                 adjacency_from_coords, to_adjacency_dict,
                                 from_adjacency_dict, rate_map)

adj = adjacency_from_grid(geodataframe, method="queen")  # or "rook"
adj = adjacency_from_bounds(bounds_df, method="rook")    # (n, 4) minx/miny/maxx/maxy
adj = adjacency_from_coords(x, y, method="knn", k=8)     # symmetrised k-nearest
adj = adjacency_from_coords(x, y, method="radius", r=2.0)
d   = to_adjacency_dict(adj)          # {"0": ["1", "5", ...], ...} for interop
adj = from_adjacency_dict(d)
```

Conventions:

* The **diagonal is 1**: a focal assemblage always includes itself,
  matching treesliceR.
* `queen` counts any boundary contact (edges or corners); `rook` requires
  a shared boundary segment of positive length. The GeoPandas path
  implements the same semantics on true polygons; without GeoPandas,
  provide `minx/miny/maxx/maxy` columns and the rectangle arithmetic is
  used.
* `from_adjacency_dict`/`to_adjacency_dict` interoperate with the R
  `AU_adj`-style sparse representation.

Mapping:

```python
ax = rate_map(rates, grid_df, rate="CpD")            # bounds columns ...
ax = rate_map(post, grid_gdf, rate="rate_median")    # ... or GeoDataFrame
ax = rate_map(post, grid_gdf, rate="hdi_high", quantiles=True)
# NaN rates render in light grey; ax is a matplotlib Axes
```

## 10. Visualisation

```python
from phyloslicer.viz import plot_rate_line, plot_sensitivity, plot_posterior

ax = plot_rate_line(profile, fit=fit, ages=stack.ages)   # cumulative + fit + pXO line
ax = plot_sensitivity(sens_df)                           # rate vs n_slices per site
ax = plot_posterior(post)                                # ECDF with HDI shading
```

> `ps.viz.plot_rate_line(..., direction="from_root")` is the default: the x-axis runs from the root towards the present, the same orientation the worked examples in this guide use. Pass `direction="from_present"` for the treesliceR view (x is age in Ma before present).

All functions return a matplotlib `Axes`, use static publication styling
by default and never call `plt.show()`.

## 11. Simulation and bundled data

```python
from phyloslicer.simulate import yule_tree, birth_death_tree, make_dataset

tree = yule_tree(100, age=1.0, seed=1)             # ultrametric, exactly n tips
tree = birth_death_tree(100, age=1.0, birth_rate=1.0, death_rate=0.2, seed=1)
# grows until 100 lineages are extant, then rescales the shape onto age=1.0
mat, coords = make_dataset(n_sites=100, n_species=50,
                           structure="clustered",  # disk-shaped ranges
                           seed=1)                 # or "random"
```

Bundled example (simulation mirroring the scale of the treesliceR
Australian-passerine case study — grid, 308 species, posterior pool):

```python
from phyloslicer.datasets import sim_passerines

data = sim_passerines(n_sites=208, n_species=308, n_trees=5, seed=42)
# data["trees"] (list[TreeArray], labels sp_1..), data["mat"] (sites x species,
# columns == tip labels), data["coords"], data["bounds"], data["adj"]
```

The columns of `data["mat"]` carry the trees' tip labels, so tree and
matrix can be used together directly:

```python
res = ps.cpd_rate(data["trees"][0], data["mat"], n_slices=50)
```

## 12. Command-line interface

```
phyloslicer [--version]
phyloslicer slice     TREE (--rootward X | --tipward X | --n-slices N)
                      [--criterion {time,pd}] --out PATH
phyloslicer rates     TREE SITES [--rate {cpd,cpe}] [--n-slices N]
                      [--criterion {time,pd}] [--percentile P] --out CSV
phyloslicer posterior TREES SITES [--rate {cpd,cpe}] [--n-slices N]
                      [--hdi H] --out CSV
phyloslicer validate  TREE SITES --reference REF_CSV [--manifest PATH]
                      [--rate-column NAME] [--n-slices N] [--tolerance T]
                      [--no-strict] --out CSV
```

* `slice` writes one Newick file, or a directory of per-slice files
  (`slice_0001.tre`, ...) for `--n-slices`.
* Exactly one of `--rootward` / `--tipward` / `--n-slices` is accepted.
* An existing `--out` is refused unless `--force` is given. That matters most
  for `--n-slices`: without the check, re-running `slice` into the same
  directory leaves the previous generation's `slice_*.tre` files behind, and
  `posterior` then pools two runs' slices while `slice` reports only the count
  it just wrote.
* `--log FILE` writes a JSON run summary (version, command, `elapsed_s`,
  `started_at`, key statistics). There is no `--seed`: the CLI performs no
  random sampling, and `sensitivity(..., seed=)` is seedable through the API.
* `validate` compares a cumulative-rate run against a treesliceR reference CSV.
  The compared column is the one the **provenance manifest** registers as
  `rate_column` (so both `CpD` and `CpE` references work, and the report names
  the columns after the index actually compared); a `--reference` with no
  manifest entry is refused, exit 2. Provide the manifest with `--manifest`,
  `$PHYLOSLICER_REFERENCE_MANIFEST`, a `manifest.txt` next to the reference, or
  `validation/reference/manifest.txt`. Beta references (`CpB`, `CpB_RW`) need an
  adjacency matrix, which this subcommand does not take — score those with
  `python validation/compare.py`.
* `--n-slices` must equal the `n_slices` recorded in the manifest header,
  because two different window counts are not the same quantity.
* Agreement rules: finite pairs compare by relative error ≤ tolerance
  (default `1e-4`, matching the optimiser agreement the shipped references
  actually reach, about 2e-6); sites NA on **both** sides count as agreement; a
  one-sided NA counts as disagreement. The run prints the observed maximum
  relative error next to the pass rate.
* Exit codes: `0` success, `1` numeric disagreement (suppress with
  `--no-strict`, for exploratory use), `2` a refused or malformed request. So
  `validate` can gate CI: an equivalence rate below 1 is a nonzero exit.

## 13. Moving from treesliceR (R)

| treesliceR | PhyloSlicer | notes |
|---|---|---|
| `squeeze_root(t, x, criterion)` | `ps.slice_rootward(t, x, criterion)` | `"my"` → `"time"`, `"PD"` → `"pd"` |
| `squeeze_tips` | `ps.slice_tipward` | idem |
| `squeeze_int` | `ps.slice_interval` | `invert=TRUE` → `invert=True`, but the main API returns `(older, younger)` where R returns `(younger, older)` — use `ps.compat.squeeze_int` for R's order |
| `phylo_pieces(n, method)` | `ps.slice_pieces(n=` or `width=)` | `method=1` → `n=`, `method=2` → `width=`; `SliceStack.ages` is decreasing Ma-before-present where R's `timeSteps` is ascending cumulative depth — use `ps.compat.phylo_pieces` for R's vector |
| equal-`"PD"` slicing | equal-`"pd"` slicing | **not numerically comparable**: see §13.1 |
| `prune_tips(method=1/2)` | `ps.prune_tips(side="after"/"before")` | `qtl=TRUE` → `quantiles=True` |
| `r_phylo(index=)` | `ps.rate_profile(index=, normalize=False)` | the R function returns raw per-slice values for PD/PE, so use `normalize=False` to reproduce it; the PhyloSlicer default returns profiles normalised to sum to 1 |
| `CpD / CpE` | `ps.cpd_rate / ps.cpe_rate` | + `converged`, `r_squared`, `reason` |
| `CpB / CpB_RW` | `ps.cpb_rate(weighted=, multisite_domain=)` | **cross-validated for the 11 beta families the reference can produce** (§15); the 3 it does not return have no reference, and the turnover/nestedness rates are a documented normalisation divergence rather than a disagreement. Under the default `multisite_domain="mixed"` the multisite series follows the reference arithmetic, which makes it drift with neighbourhood size (§6.1); `multisite_domain="paired"` is an **opt-in pair-summed domain** (a parameter value, not a separate `approach`/"paired method") that removes the drift at the cost of no longer being comparable with R. The `approach="pairwise"` aggregation uses the pairwise mean because R's `CpB` pairwise branch is dimensionally inconsistent - measured: that branch aborts inside its own `nls` call |
| `CpR_sensitivity(_plot)` | `ps.sensitivity` / `ps.viz.plot_sensitivity` | `sample_sites=0` = all sites; + suggestion |
| `CpR_graph` | `ps.spatial.rate_map` / `ps.viz.plot_rate_line` | matplotlib backend |
| — | `ps.posterior_rate / posterior_origin` | new: posterior uncertainty |
| — | `ps.bootstrap_sites` | new: site resampling |
| — | adjacency builders, CLI, GeoPandas interop | new |

A drop-in compatibility layer keeps the R-style signatures:

```python
from phyloslicer.compat import squeeze_root, squeeze_tips, squeeze_int, phylo_pieces, prune_tips

pieces, time_steps = phylo_pieces(tree, n=100, criterion="my", method=1, timeSteps=True)
# time_steps is R's ascending cumulative depth from the root, cumsum(rep(n, nslices)):
# time_steps[j] is the depth reached BY pieces[j], and time_steps[-1] == tree depth.
# (The main API's SliceStack.ages is the opposite convention; ask for it with
#  phylo_pieces(..., timeSteps_as_depth=False).)

complement = squeeze_int(tree, from_=7.0, to=1.0, invert=True)
# complement[0] == squeeze_root(tree, 1.0)   -> the younger/crown piece (R's tree1)
# complement[1] == squeeze_tips(tree, 7.0)   -> the older piece        (R's tree2)
```

Both orderings are asserted by `tests/test_compat.py`.

### 13.1 Behavioural differences that are *not* just renames

* Silent `NA` rows are replaced by `converged=False` plus a `reason`.
* Every tree/matrix adjustment (dropped species, pruned tips) warns.
* The CLI, seeds and JSON run summaries make runs auditable.
* `n_slices` corresponds to R's `n` (number of slices), not a width.
* **Equal-PD slicing (`criterion="PD"`) is a deliberate, unvalidated
  deviation.** R hard-codes the branch counts as `2, 3, 4, …`
  (`squeeze_root.R:98`, `phylo_pieces.R:159`:
  `nBranch = 2:(length(unique(nodes$YearBegin)) + 1)`), which assumes a strictly
  binary tree with all-distinct node depths. On a multifurcating tree the
  arithmetic does not even conserve total PD: for
  `(((A:2,B:2):1,C:3,D:3,E:3):1,(F:4,G:4):0);` (T = 4, total PD = 23) R's cut
  points accumulate 13 PD where the tree has 23. PhyloSlicer counts the
  *actually active* branches per depth interval, so its `"pd"` windows cumulate
  to the full tree PD by construction. That is a correction, not a port, and it
  is **not cross-validated** — every reference file in `validation/reference/`
  was generated with `criterion = "my"` (see `validation/NOTES.md`).
* `squeeze_root`/`squeeze_tips`/`squeeze_int`/`phylo_pieces` reject a
  `criterion` outside `{"my", "million years", "pd", "time"}` (case-insensitive)
  and `phylo_pieces`/`prune_tips` reject a `method` outside `{1, 2}`, instead of
  R's silent "no cut" / `NULL`.
* The multisite β domain is opt-in changeable (§6.1); the default keeps R's
  arithmetic.

## 14. Diagnostics, statuses and error messages

### 14.1 `pb_per_slice` status strings

| status | meaning |
|---|---|
| `ok` | neighbourhood analysed normally |
| `fewer_than_two_sites` | focal row of `adj` has < 2 sites (incl. self) |
| `degenerate_neighbourhood` | all sites empty or all sites share the full composition (beta undefined) |
| `zero_denominator` | a beta denominator vanished (numerically degenerate) |

Skipped focal sites come back as `converged=False` rows with the status in
`reason` — never silent `c(NA, NA, NA)`.

### 14.2 `RateFit.reason` values

`"fewer than two slices"`, `"non-finite input values"`,
`"profile sums to zero (empty or degenerate assemblage)"`,
`"no start converged"`, `"optimiser did not report convergence"`,
`"vectorised Newton did not satisfy criteria"`.

### 14.3 Common validation errors

* `slice_rootward(tree, t)` with `t > tree.root_age` → *"the threshold
  inputted ... is bigger than the available by the phylogeny"*.
* `slice_interval(tree, a, b)` with `a <= b` (time criterion) → *"the
  thresholds set in arguments [start] and [stop] are incompatible"*.
* `slice_pieces(tree)` without `n`/`width`, or with both → *"pass either n
  or width"*.
* `sensitivity(..., rate="cpb")` without `adj` → *"rate=cpb requires an
  adjacency matrix"*.
* `cpb_rate(...)` with an `adj` whose row count differs from the matrix →
  *"adj must have one row per site of mat"*.

## 15. Numerical validation

The comparison set is `{PD, PE, CpD, CpE}`, plus the beta-diversity families the
reference can return on a given input, run against treesliceR
1.1.0 strictly as a black box (R script in `validation/r/`, comparison in
`validation/compare.py`, provenance in `validation/reference/manifest.txt`,
tolerance rationale and gap list in `validation/NOTES.md`). Within that set,
only the `CpD` and `CpE` files are genuinely validated against R today:

* CpD and CpE rates: maximum relative error ≈ 2 × 10⁻⁶ across all
  assemblages (tolerance 10⁻⁴; the residual gap is optimiser tolerance); the PD
  and PE *totals* in the same files agree to ≈ 4 × 10⁻¹³. These are the only
  genuinely external comparisons the project can currently make.
* per-slice PD and PE profiles (`*__rPD.csv` / `*__rPE.csv`): these match to
  ≈ 3.4 × 10⁻¹¹ over 18000 comparisons, and they are a genuine Python-vs-R
  comparison. They used to be the exception: the archived copies were stored
  slice-by-site with 30 slice columns, a layout the committed driver (50 slices,
  site-by-slice) cannot emit, so the manifest declared them `status = unverified`
  and `validation/compare.py` skipped them rather than silently transposing them.
  Running `validation/r/generate_reference.R` under R 4.5.3 (aarch64-apple-darwin20.0.0) with ape 5.8.1 and treesliceR 1.1.0 settled it -
  every other reference file came out byte-identical, and the profiles written by
  that run are registered now. The superseded copies are not shipped; the
  re-run itself is the evidence.

What this comparison does **not** cover (all of it is spelled out in
`validation/NOTES.md`; do not read the four numbers above as a general
equivalence proof):

* **3 beta families** (`CpB` Sorensen/pairwise, and `r_phylo(index = "PB")` for
  Sorensen- and turnover-pairwise) — treesliceR 1.1.0 returns nothing comparable for them
  on these inputs, so there is no reference to score. The other 11 beta families
  **are** compared: 24 per-slice cumulative-β profile rows to 1.1e-9,
  1260 whole-tree `PB`/`PB_RW` totals to 7.2e-13, and 540
  verified rate fits to 2.9e-6 at the optimiser tolerance.
* **the turnover and nestedness beta *rates*** — compared, and diverging by up to
  196% in relative terms:
  this package renormalises the profile to sum to 1, the reference fits
  `cumsum(profile) / PB_whole_tree`. Each such row is scored on both curves, and the
  reference-convention refit reproduces R, so the gap is the pinned normalisation choice
  (see validation/NOTES.md), not the arithmetic; the whole-tree `PB` of the same file agrees to
  7.2e-13.
* **equal-`"PD"` slicing** — every scenario runs `criterion = "my"`. The `"PD"`
  criterion deviates from R on purpose (§13.1) and has never been compared.
* **the tip-pruning path** — the generator guarantees every species occurs at
  least once, so `align_tree_matrix` never prunes and that branch is not
  exercised.
* **degenerate-assemblage NA patterns** — `compare.py` now compares the NA
  *patterns* element-wise rather than counting non-finite values, but the
  reference generator gives every species at least one occurrence, so there are
  no NA rows on either side and the check has nothing to bite on within the
  archived set. Any earlier claim that "NA patterns coincide on every
  degenerate assemblage" is unsupported and should be read as removed.
* The reference trees are the ones the shipped generator emits (20/100-tip
  birth-death and coalescent trees, unit age, `n_slices = 50`, 30 sites);
  nothing here validates a tree of a different shape or depth.

The test suite additionally pins analytic slice solutions, brute-force index
cross-checks, per-slice partition invariants, compat-layer return orders, and
I/O round-trips. Its size is deliberately not quoted here: a number in a manual
goes stale the moment anyone adds a test. Ask the tree instead -

```bash
grep -c "^[[:space:]]*def test_" tests/test_*.py                               # functions
python -m pytest --collect-only -q -o addopts="" 2>/dev/null | grep -c "::"     # cases
```

- the first counts `def test_` functions, the second counts what pytest actually
*collects*, and the two differ because parametrised functions expand into one
case per parameter set. `-o addopts=""` is required: `pyproject.toml` sets
`addopts = "-q"`, so a bare `pytest --collect-only -q` runs at `-qq` and prints one
count per file rather than node IDs (the `grep` then returns 0).
Run everything with:

```bash
python -m pytest tests/
python validation/compare.py validation/reference
```

## 16. FAQ and troubleshooting

**Q: My matrix columns don't match the tree tips.**
A: `align_tree_matrix` (called automatically by the rate functions) fixes
the overlap and warns about every adjustment. To make the alignment
explicit — or to silence avoidable warnings — subset your matrix first.

**Q: Can I use a non-ultrametric tree?**
A: Not by default: `slice_pieces()` raises `ValueError`, because slicing
semantics assume all tips end at the same age. Pass
`allow_non_ultrametric=True` for the warn-and-go behaviour, which records
`ultrametric = False` in the result table but still leaves the slice ages
without their "diversity before time *t*" meaning. Rescale or re-ultrametrise
the tree first.

**Q: Some sites return `converged=False` with "profile sums to zero".**
A: The assemblage is empty. That is reported, not fatal. Fill or drop the
site depending on your question.

**Q: The CpD of a site is huge / the pDO is absurd.**
A: Check `r_squared` and `plot_rate_line(profile, fit, ages)`. Nearly-flat
or single-slice-dominated profiles push the rate to the bound (50); such
boundary fits are flagged through the diagnostics rather than hidden.

**Q: How many slices should I use?**
A: Run `sensitivity` over a grid; use `suggested_slices` per site, and
report the grid in your methods.

**Q: Where did the tree's support values go?**
A: `TreeArray` stores tips, edges and lengths only; internal labels are
dropped on read. Keep your original Newick/NEXUS for support values.

**Q: Is `percentile=0.95` the same as `pDO = 5` in R?**
A: Yes — `pXO = −ln(1 − 0.95)/r = −ln(0.05)/r`, which is R's
`−log(pDO/100)/r` with `pDO = 5`.

**Q: Is there built-in parallelism?**
A: Parallelism is internal to the vectorised kernel — the single sparse product
`A @ C` already computes all sites simultaneously; use a workflow engine to
parallelise across trees or chunks if needed.

## 17. API quick reference

```python
# core
ps.TreeArray.from_newick(x) / .from_nexus(p) / .from_dendropy(o) / .from_biopython(o)
tree.root_age / .n_tips / .tip_labels / .total_pd() / .is_ultrametric
tree.subarray(labels) / .to_newick() / .to_dendropy() / .to_biopython()

# slicing
ps.slice_rootward(tree, time, criterion="time", collapse=False)
ps.slice_tipward(tree, time, criterion="time", collapse=False)
ps.slice_interval(tree, start, stop, invert=False, criterion="time", collapse=False)
ps.slice_pieces(tree, n=None, width=None, criterion="time")  -> SliceStack
ps.prune_tips(tree, threshold, quantiles=False, side="after")
SliceStack: .windows .ages .widths .contribution(subset) .contribution_matrix()
            .to_sliced_trees(collapse=False)

# indices
ps.indices.pd(tree, mat) / .pe(tree, mat) / .pb(tree, mat, adj, component, approach, weighted)
ps.indices.pd_per_slice / .pe_per_slice / .pb_per_slice
ps.indices.site_edge_membership / .edge_range_sizes

# rates
ps.fit_rate(profile, ages, method="lsq", multistart=8, bounds=(1e-6, 50)) -> RateFit
ps.origin_time(rate, percentile=0.95)
ps.cpd_rate(tree, mat, n_slices=100, criterion="time", percentile=0.95, ...)
ps.cpe_rate(...);  ps.cpb_rate(tree, mat, adj, component=..., approach=..., weighted=..., multisite_domain=...)
ps.rate_profile(tree, mat, adj=None, n_slices=100, index="PD")
ps.sensitivity(tree, mat, adj=None, slice_grid=(...), rate="cpd", sample_sites=0, seed=None)

# uncertainty
ps.posterior_rate(trees, mat, rate="cpd", n_slices=100, hdi=0.95)
ps.posterior_origin(trees, mat, rate="cpd", percentile=0.95, hdi=0.95)
ps.uncertainty.bootstrap_sites(tree, mat, n_boot=1000, rate="cpd", seed=None)
ps.uncertainty.iter_trees(source)

# spatial
ps.spatial.adjacency_from_grid(gdf, method="queen")
ps.spatial.adjacency_from_bounds(bounds, method="rook")
ps.spatial.adjacency_from_coords(x, y, method="knn", k=8, r=None)
ps.spatial.to_adjacency_dict(adj) / .from_adjacency_dict(d)
ps.spatial.rate_map(rates_df, grid, rate="CpD", quantiles=False)

# viz / simulate / datasets / compat
ps.viz.plot_rate_line / .plot_sensitivity / .plot_posterior
ps.simulate.yule_tree / .birth_death_tree / .make_dataset
ps.datasets.sim_passerines()
ps.compat.squeeze_root / .squeeze_tips / .squeeze_int / .phylo_pieces / .prune_tips

# io
ps.io.read_newick / .parse_newick / .read_nexus / .write_newick / .write_nexus
ps.io.presence_matrix / .align_tree_matrix / .as_matrix
```

## Citation

If you use PhyloSlicer, please cite both this package and the
methodological source:

> Araujo, M.L., Ferreira, L.G.S.S., Nakamura, G., Coelho, M.T.P., Rangel, T.F. (2025)
> 'treesliceR': a package for slicing phylogenies and inferring phylogenetic patterns
> over evolutionary time. *Ecography* 2025, e07364. https://doi.org/10.1111/ecog.07364
>
> **On the year.** The article was first published online on 28 October 2024 and forms part of the 2025 volume of *Ecography* as article e07364; it is cited here by the volume year. Verified against the Crossref record for doi:10.1111/ecog.07364 on 26 September 2026 (journal-article, Wiley, volume 2025, article-number e07364, issued 2024-10-28, print 2025-01); the five authors it lists, in order, are the ones in the DESCRIPTION of the treesliceR 1.1.0 build that generated the validation references. The package's own vignette still carries an older "in review" placeholder, which the journal record supersedes. If your citation style follows the online-first date, cite 2024 instead - the DOI is unchanged.

## License

Code: BSD-3-Clause. Documentation: CC-BY 4.0.
