# Cross-validation notes (step 4 of the validation protocol)

## Tolerance rationale

| quantity | tolerance | why |
|---|---|---|
| per-slice PD / PE profiles (`r_phylo` raw slices) | rtol 1e-8 | fully deterministic arithmetic; residual differences come only from floating-point summation order (sparse product vs per-edge sums) and the 12-digit Newick round-trip of the tree |
| CpD / CpE rates | rtol 1e-4 | the *rate* is the solution of a non-linear optimisation: treesliceR fits `stats::nls` (single start r = 0.2, Gauss-Newton) while phyloslicer uses bounded multi-start least squares. Both converge to the same optimum of the same objective, but to optimiser-tolerance accuracy (~1e-5 relative), so 1e-4 is the appropriate equivalence band |
| pDO / pEO, pBO | follows the rate | derived as `-ln(1 - percentile) / rate` |
| per-slice PB profiles (`r_phylo(index = "PB")`) | rtol 1e-8 | the same deterministic arithmetic as PD/PE, driven by the queen adjacency the generator exports as an input rather than rebuilds from the matrix |
| CpB / CpB_RW rates | rtol 1e-4 | the same optimiser difference as CpD/CpE. Where the two implementations normalise the fitted curve differently (see *Beta-diversity divergences* below) the row is scored twice: once on the shipped profile and once on the reference's own curve, and only the second is allowed to close the gap |
| PB / PB_RW whole-tree totals | rtol 1e-8 | the reference reports them next to its own fit, so they are a deterministic quantity, not an optimiser output |

## Degenerate cases

treesliceR returns `c(NA, NA, NA)` for empty assemblages and lets `nls`
errors abort the whole call. phyloslicer returns `converged=False` with an
explicit `reason` and `NaN` values. The comparator treats an NA on both
sides as agreement and reports NA-pattern differences in the `note` column.

## Beta-diversity divergences (documented design decisions)

treesliceR's `CpB` pairwise branch builds a k x P matrix of per-pair relative
values and then feeds `cumsum()` of it into `nls` together with a length-k
`age` vector - dimensionally inconsistent.
`CpB_RW`'s pairwise branch instead averages the per-pair ratios within each
slice. PhyloSlicer adopts the pairwise mean for every component/approach
combination (single code path, documented in the docstring of
`indices.api._pb_series`).

That reading of the source is now a measurement: in the 2026-09-21 run under R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 the Sorensen
pairwise branch fails inside its own optimiser with `task 1 failed - "'qr' and
'y' must have the same number of rows"` (the per-pair vector at `CpB.R:479`
recycled against the slice vector at `CpB.R:525`), and the pairwise
`r_phylo(index = "PB")` output for Sorensen and turnover comes back as a
flattened k x P block whose elements are 300 / 750 / 1800 long rather than one
50-slice curve per assemblage.  `write_beta_rates()` / `write_beta_profile()`
in the generator refuse those shapes instead of forcing them into a reference,
so the three families are *absent* from `validation/reference/` and the report
carries them as `no reference` rows - a named gap, never a silent pass.

The second divergence is a convention, and it is the reason a beta rate row can
read `documented divergence` while everything around it verifies: the reference
fits `cumsum(profile) / PB_whole_tree`, whereas PhyloSlicer renormalises the
profile to sum to 1 before fitting (a pinned design decision).  On turnover and
nestedness, in both approaches, that costs up to a factor of 1.96 in the fitted
rate.  The comparator therefore refits each of those rows on the reference's own
curve (`fit_reference_convention`); the refit reproduces the reference rate, and
the whole-tree `PB` / `PB_RW` column of the very same reference file agrees with
ours at 7.2e-13.  What that row documents is the normalisation choice, not a
difference in the arithmetic.

> **What is still provisional.**  The multisite beta *index* grows with
> neighbourhood size, because the reference convention sums the shared-term
> component over sites while the turnover components are summed over pairs.
> PhyloSlicer follows it deliberately; values from neighbourhoods of different
> sizes must not be compared with one another.

## Verdict

Run `python validation/compare.py` after `Rscript validation/r/generate_reference.R`
and append the printed summary here, together with the scenario list used.

### Current status

- **PD/PE per-slice profiles** (rPD/rPE): **validated.** Ran on 2026-09-21 under `micromamba env r-4.5.3`: R 4.5.3 (aarch64-apple-darwin20.0.0), ape 5.8.1, treesliceR 1.1.0. The
  profiles that run writes are registered in
  `validation/reference/manifest.txt` (`site-by-slice`, 30 rows x 50 slice
  columns, sha256 pinned) and `compare.py` scores them: max relative error
  **3.396e-11 over 18000 comparisons** against the 1e-8 deterministic tolerance,
  100% pass rate, with **157** of those cells bit-for-bit identical - the
  `identical_cells` column of the report, counted on the same finite cells the
  maximum is taken over. Until that run there was no way
  to establish the archived files' provenance (item 4) and these two families
  were declared `status = unverified` and skipped; the history is kept in the
  manifest's `limitations` block rather than edited out.
- **CpD / CpE rates**: genuinely compared, max relative error ~2e-6 (tolerance
  1e-4); the PD and PE *totals* in the same files agree to ~4e-13. These are the
  only two *alpha* rate families; the beta rates are `CpB` and `CpB_RW`, and the
  full family list is in the list below.
- **CpB, CpB_RW and `r_phylo(index = "PB")`**: **cross-validated where the
  reference can produce a counterpart.**  Eleven of the fourteen beta families
  have a registered reference and are scored: the 24 per-slice cumulative-beta
  profile rows agree to **1.08e-09 over 36000 values** (tolerance 1e-8; 10244 of
  those cells bit-for-bit identical, most of them the numerically-zero oldest slices),
  the 42
  whole-tree `PB`/`PB_RW` rows to **7.2e-13 over 1260**, and the 18 verified rate
  rows to **2.88e-06 over 540 assemblage fits** (tolerance 1e-4).  Twenty-four
  further rate rows - turnover and nestedness, both approaches - are recorded as
  `documented divergence` for the normalisation reason above, and eighteen rows
  (three families x six scenarios) as `no reference` because treesliceR 1.1.0
  returns nothing comparable.  None of the three is counted as a pass, and
  `compare.py --strict` fails on all of them.
- **`r_phylo` raw PD/PE/PB profiles**: these *are* the `:PD` / `:PE` rows of
  the equivalence report, and after the 2026-09-21 run under R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0 they are
  genuine Python-vs-R comparisons (3.396e-11 max relative error), not
  Python-vs-a-file-of-unknown origin - item 4 records what changed.
  `r_phylo(index = "PB")` is compared too, wherever the reference returns one
  curve per assemblage (see the beta bullet above and item 2).
- **`criterion="PD"` (equal-PD slicing)**: **still not cross-validated**, and
  now for a demonstrated reason rather than an absence of tooling. The
  generator's guarded `criterion = "PD"` block *was* run (the 2026-09-21 run under R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0) and
  failed on all six scenarios (18 warnings, no `*__PD.csv` written).  The failure
  is `Error: object 'cutted_tree' not found`, not the
  `length of 'dimnames' [2] not equal to array extent` an earlier revision of
  these notes recorded: `phylo_pieces.R:42` and `:152` branch on `criterion ==
  "pd"` in lower case while `CpD`/`CpE`/`r_phylo` document and forward `"PD"`, so
  no slice is ever computed.  The documented spelling therefore reaches no
  comparison at all.  (Whether the undocumented lower-case spelling *should* be
  used is a separate decision, and the structural divergence in item 3 applies to
  it too.)  Nothing may be claimed about equal-PD slicing numerically in either
  direction from these references; the structural divergence described in item 3
  stands as a reading of the source, not as a measurement.
- **Pruning path** (matrix missing entire clades): not covered by the
  current validation set; the reference generator guarantees every species
  occurs at least once, so `align_tree_matrix` never triggers pruning.
- **Comparison set is `{PD, PE, CpD, CpE}` plus the beta set** (`CpB`,
  `CpB_RW`, `r_phylo(index = "PB")`) **for the families treesliceR can produce.**
  Anything outside that, equal-PD slicing included, has no external reference,
  whatever this file's headline numbers suggest.


## Not covered by this cross-validation (read before quoting any pass rate)

1. **Three beta families have no external reference**: `CpB`
   `sorensen`/`pairwise`, and `r_phylo(index = "PB")` for `sorensen`/`pairwise` and
   `turnover`/`pairwise`. They are not gaps in the harness - treesliceR 1.1.0 does not
   return a comparable object for them (measured, see *Beta-diversity divergences*), and
   the generator refuses to write a shape it would have to invent. Every statement about
   those three remains an internal-consistency result.
2. **Cumulative-β references now exist, for the families the reference
   produces.** `generate_reference.R` calls `CpB`, `CpB_RW` and
   `r_phylo(index = "PB")` for every component-by-approach combination, exports the
   queen adjacency it used as an input (`*__adj.csv`, so the Python side is scored on
   exactly the neighbourhood sets the reference consumed), and registers what came
   back: 84 of the 126 beta rows in the report are verified. The three families item 1
   names produce no file at all, and the four rate families of the A2 normalisation
   convention are scored as divergences rather than passes. The `r_phylo` PD/PE
   outputs *are* written out and back the report's `:PD` / `:PE` rows with a genuine
   Python-vs-R comparison - see item 4 for how that was settled on 2026-09-21.
3. **`criterion="PD"` (equal-PD slicing) was never compared.** All six archived
   scenarios run `criterion="my"`. The deviation is deliberate and documented:
   PhyloSlicer counts real active branches, so total PD is conserved even on
   multifurcating trees, whereas the reference hard-codes branch counts
   `2, 3, 4, ...` per depth (`squeeze_root.R:98`, `phylo_pieces.R:159`,
   `nBranch = 2:(length(unique(nodes$YearBegin)) + 1)`) and therefore does not
   even conserve total PD on such a tree — e.g. on
   `(((A:2,B:2):1,C:3,D:3,E:3):1,(F:4,G:4):0);` (T = 4, total PD = 23) the
   reference cut points accumulate 13 while PhyloSlicer accumulates 23.
   `validation/r/generate_reference.R` *attempts* `CpD`/`CpE`/`r_phylo`
   references under `criterion = "PD"` (files suffixed `__PD.csv`). **They could
   not be produced**: the block was run (the 2026-09-21 run under R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0) and every call failed
   inside treesliceR with `object 'cutted_tree' not found` (the lower-case
   `criterion == "pd"` branch problem described under *Verdict*), so no
   `*__PD.csv` exists in `validation/reference/` and none may be added to the
   pass-rate table. Even if a future treesliceR
   version produces them, the warning block in the script still applies: the
   archived scenario trees are binary, so agreement there would not validate
   the multifurcating case where the two implementations diverge.
4. **Provenance of `*__rPD.csv` and `*__rPE.csv`: RESOLVED on 2026-09-21 by
   running the generator.** Ran on 2026-09-21 under `micromamba env r-4.5.3`: R 4.5.3 (aarch64-apple-darwin20.0.0), ape 5.8.1, treesliceR 1.1.0. `Rscript
   validation/r/generate_reference.R <dir>` reproduced **21 of the 21** time-criterion files it
   shares with `validation/reference/` byte-for-byte - the three Newick trees,
   the six presence matrices and all twelve `*__CpD.csv` / `*__CpE.csv`
   references - which pins the seed, the RNG stream and the package versions
   behind the archived set. The only outputs that differed were the twelve
   profile files: the generator writes **30 rows × 51 columns** (`site_id` +
   `slice_1..slice_50`), i.e. **site-by-slice with 50 slice columns**, exactly
   as its source says, while the committed files were **50 rows × 31 columns**,
   **slice-by-site with 30 columns** and a `site_id` that runs 1..30 then 1..20
   (the R recycling signature of a `seq_len(nrow(mat)) = 30` vector glued onto
   a 50-row block). So the committed profiles came from an older revision of
   this script and were never the companion of the registered inputs; the
   regenerated ones are, and `manifest.txt` now registers them as
   `status = registered`, which is what turns the report's `:PD` / `:PE` rows
   from "Python agrees with these CSV files" into "Python agrees with
   treesliceR". What the manifest still cannot do: its hashes detect
   after-the-fact edits, they do not by themselves re-derive that a CSV is what
   the committed script writes - only the run above does that, and it is
   recorded in each profile's `note`. The superseded files are not part of the
   repository: they exist only in the authors' working history, so the
   contradiction is settled by the re-run above rather than by exhibiting the old
   bytes.
   The original measurement, kept for the record: each committed file is
   **50 rows × 31 columns**
   (`site_id` + `slice_1..slice_30`), i.e. stored **slice-by-site with 30 slice
   columns**, because the archived scenario set uses 30 sites. The committed
   `generate_reference.R`, however, sets `n_slices <- 50L` and writes
   `matrix(unlist(rpd_ref), nrow = length(rpd_ref), byrow = TRUE)`, which for a
   30-site matrix emits **30 rows × 51 columns** — **site-by-slice with 50 slice
   columns**. The two orientations and the two column counts are different, so
   either a different driver produced the archived files, or they were transposed
   from PhyloSlicer's own output - and the 2026-09-21 run says the former. An
   earlier revision of `validation/compare.py` transposed an out-of-shape
   reference automatically ("legacy column-bound layout"), which is what let the
   mismatch pass unnoticed for as long as it did. `manifest.txt` pins sha256,
   rows, cols, orientation and rate/total/origin columns for every shipped
   file.
5. **The pruning path is unit-tested but never exercised by the R comparison**, because
   every fixture matrix contains every tree tip.
6. **NA-pattern agreement on degenerate assemblages** is now a real comparison of the
   patterns (not of the count of non-finite values) in `validation/compare.py`
   (`na_pattern_report`, element-wise), but the six reference scenarios are constructed
   so that no assemblage is empty — every archived `*__CpD.csv` / `*__CpE.csv` is fully
   finite (measured: 0 NaN in `CpD`/`PD`/`pDO` and `CpE`/`PE`/`pEO` throughout) — so
   this check has nothing to bite on within the archived reference set.
7. **The R and package versions are now recorded, but not checked.**
   `validation/reference/manifest.txt` carries `written_in`, `r_available` and
   `generator_sha256`, and each regenerated profile's `note` names
   `R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0`; the vendored copy's `DESCRIPTION`
   says `Version: 1.1.0`. Nothing *asserts* those versions at run time, so
   "validated against treesliceR 1.1.0" is a documented property of the
   reference set rather than an enforced precondition of the comparator - a
   future run under a different treesliceR could silently disagree with the
   pinned bytes. The sha256 pins catch that (a different output fails
   `check_registration`), which is the weaker but automatic guarantee.

A pass rate of 100% in this directory therefore means: *phyloslicer reproduces the numbers
in the archived `criterion = "my"` reference CSVs to the stated tolerance, within the four
alpha families `{PD, PE, CpD, CpE}` and the eleven beta families treesliceR can produce* -
the per-slice PD/PE profiles included, since the generator was run and the manifest
registers them (item 4). It does not extend to equal-PD slicing, the pruning path,
degenerate-assemblage NA patterns, or the three beta families of item 1, and the four
rate families of the A2 normalisation convention are reported as divergences rather than
as passes.

## Reproducibility of the report itself

`validation/reference/equivalence_report.csv` is a byte function of the reference
directory plus the two Newick/matrix inputs per scenario, and nothing else:

- re-running `Rscript validation/r/generate_reference.R <dir>` (the 2026-09-21 run under R 4.5.3 / ape 5.8.1 / treesliceR 1.1.0)
  rewrites all **105** data files of `validation/reference/` byte-for-byte - the
  three trees, six matrices, six adjacencies, twelve CpD/CpE rates, twelve rPD/rPE
  profiles, and the whole beta block (CpB, CpB_RW, rPB).  That is what makes the
  beta comparison a comparison with treesliceR rather than with these CSVs.
- `compare.py` pins its BLAS thread counts (`OMP_NUM_THREADS` and friends, set
  before numpy loads) because without the pin the report is *not* byte-reproducible:
  measured over consecutive runs, 7 of the 162 rows move in their last significant
  digits (a profile error going 9.161e-14 -> 9.141e-14, a fitted rate going
  1.47158405e-06 -> 1.47158404e-06) since multi-threaded reductions sum in a
  scheduling-dependent order.  No status changes, but the documentation quotes those
  digits, so the tolerance question is not the point - reproducibility is.
  `tests/test_validation_chain.py::TestReportDeterminism` runs the comparator twice
  in fresh processes and diffs the bytes.
- byte identity holds *within* an interpreter, not across them.  Measured against
  the numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3 environment: all 162 statuses come
  out the same and 7 rows move in their last digits (9.161e-14 -> 9.141e-14), so
  the claim the shipped report supports is "every row has the same status and the
  same error to the digits that matter", with `written_in` above naming the
  environment the archived figures were produced in.
