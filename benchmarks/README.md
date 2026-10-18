# Benchmark tables in this directory

Four independent protocols live here.  They are **not** interchangeable, two
column names mean different things in different files, and every row now names
the input it was measured on.  Every documented speed-up is traceable to one of
the tables below, and `recompute_ratios.py` prints that trace.

## Protocols

| files | what is timed | statistic | repeats | cross-language? |
|---|---|---|---|---|
| `results.csv` | `bench_core.py`: vectorised kernel vs the naive in-process baseline | **min** | 3 (naive baselines: 1) | **no** - Python vs Python |
| `results_py.csv`, `results_r.csv` | slicing kernel + CpD/CpE fits on the validation-scenario inputs | median | 5 | yes |
| `results_scaling_py.csv`, `results_scaling_r.csv` | **slicing kernel only** (`slice_pieces(...).contribution_matrix()` vs `phylo_pieces(...)`), every tree size at each of k = 25/50/100 | median | 15 (+1 discarded warm-up) | yes |
| `results_py_large.csv`, `results_r_large.csv` | end-to-end CpD on 250 assemblages, 1000/2000 tips, k = 100 | median | 5 | yes |
| `results_memory_py.csv` | `benchmark_memory.py`: tracemalloc peak of the `E x k` array vs `k` array-level copies | single run, no repeats | 1 | **no** - Python-side emulation only |

`slice_s` means **the E x k contribution matrix** in `results_py.csv`,
`results_scaling_py.csv` and `results_py_large.csv`, but **window construction
alone** in `testing_suite/benchmarks/benchmark_results.csv`.  Do not compare the
two across directories.

Every timing table carries `reps`, the number of timed runs the row is the median
of; `recompute_ratios.py` refuses a table whose `reps` is missing or disagrees
with the protocol documented above.

## What each row is - and what it is not

`results_py.csv` and `results_r.csv` hold **12 rows for 4 (tips, n_sites,
n_slices) cells**: the grid is 3 trees (`bd100`, `bd20`, `coal100`) x 2 matrix
structures (`clustered`, `random`) x k = 25/50, and `bd100` and `coal100` both
have 100 tips.  Those extra rows are **not** the five repeats.  Each row is
already the median of its own five timed runs; the rows differ because they were
measured on a *different tree and a different matrix*.  So there is no honest
"median over the four rows at (100, 30, 25)" - that would average two birth-death
and two coalescent sessions into a statistic nobody measured - and the pairing is
one-to-one on `(tips, n_sites, n_slices, source, structure)` instead.  Earlier
revisions of this directory indexed by `(tips, n_sites, n_slices)` and let the
**last** row win, which is why the archived grid ratios read from
`coal100__random` rather than from any declared scenario.

`results_scaling_r.csv` holds **21 rows for 18 (tree size, k) cells** for the same
reason in a different form: `tips = 100` is genuinely measured twice, on `bd100`
and on `coal100`, and only `bd100` exists in the Python table.

## Provenance of the `source` column

Both timing drivers write a `source` column naming the input tag, and the two
spellings are identical (`benchmark_scaling.py` used to write `bd100__tree`, the
R driver `bd100`; both now write `bd100`).  The archived tables predate it, so the
tags recorded here were **attributed from the drivers' enumeration order**, which
reproduces the archived row order exactly:

| tips | source in `results_scaling_r.csv` | source in `results_scaling_py.csv` | evidence |
|---|---|---|---|
| 1000 | `large1000` (rows 1-3) | `large1000` | `benchmarks/large/large1000__tree.tre`, written by `make_large_inputs.py` (`yule_tree(1000, age=1, seed=1100)`) |
| 2000 | `large2000` (rows 4-6) | `large2000` | `benchmarks/large/large2000__tree.tre`, `seed=2100` |
| 200 | `large200` (rows 7-9) | `large200` | `benchmarks/scaling/large200__tree.tre`, written by `make_scaling_inputs.py` (`seed=300`) |
| 500 | `large500` (rows 10-12) | `large500` | `benchmarks/scaling/large500__tree.tre`, `seed=600` |
| 100 | **`bd100` (rows 13-15) and `coal100` (rows 19-21)** | `bd100` (rows 4-6) | `validation/reference/bd100__tree.tre` and `coal100__tree.tre`, both 100 tips / 198 edges |
| 20 | `bd20` (rows 16-18) | `bd20` (rows 1-3) | `validation/reference/bd20__tree.tre` |

* The R driver walks `sort(list.files(dirs, "__tree\\.tre$"))` over the three
  directories, which orders the full paths `benchmarks/large/large1000`,
  `benchmarks/large/large2000`, `benchmarks/scaling/large200`,
  `benchmarks/scaling/large500`, `validation/reference/bd100`,
  `validation/reference/bd20`, `validation/reference/coal100` - the archived
  `tips` sequence 1000, 2000, 200, 500, **100, 20, 100** is that order.
* The Python driver sorts by `(edge count, file name)` and keeps one row per
  (tips, k).  `bd100` and `coal100` have the same 198 edges, so the name
  tie-break puts `bd100` first, and the archived `tips` sequence 20, 100, 200,
  500, 1000, 2000 is that enumeration after the duplicate is dropped.  The
  tie-break is now explicit in `benchmark_scaling.py` (`PUBLISHED_DUPLICATE`)
  instead of a side effect of `drop_duplicates`.

Nothing here needed an `unknown` tag: every rung is attributable from the drivers
plus the trees on disk.  `recompute_ratios.py` still treats `unknown` (or an
absent tag) as a hard error, so a future table cannot lean on that tolerance.

## Which published number comes from which table

| published claim | table pair | statistic | value |
|---|---|---|---|
| "76-556x" headline speed-up | `results_scaling_py.csv` / `results_scaling_r.csv` | `slice_s`, all 18 cells, `source` = `bd100` at 100 tips | 75.6 - 556.0 |
| "148-302x at 100 slices" fixed-k scaling | same pair, `n_slices = 100` (6 cells) | `slice_s` | 148.2 - 302.2 |
| IQR factor 4.3 on the widest endpoint | `results_scaling_r.csv`, cell `source=large1000`, tips=1000, k=25 | `slice_q75 / slice_q25` | 4.29 |
| "6.7-8.1x" end-to-end rate fits | `results_py_large.csv` / `results_r_large.csv` | `cpd_s`, tips=1000/2000 | 6.72 - 8.06 |
| peak memory 1.93 / 3.49 / 6.69 MB vs 2.81 / 5.62 / 11.23 MB | `results_memory_py.csv` | `array_peak_bytes`, `copies_peak_bytes`, **MB = 10^6 bytes** | see "Units" below |
| no published claim | `results_py.csv` / `results_r.csv` | `slice_s` over 12 one-to-one cells (45-161x) | reported, not quoted |
| no published claim | `results.csv` | Python-vs-Python baselines | not a cross-language number |

`recompute_ratios.py` prints `source: <python table> / <R table>` in front of
every claim and names both tables again on every ratio line, so the mapping above
is machine-enforced rather than a prose promise.

## Units: `results_memory_py.csv`

`array_peak_bytes` and `copies_peak_bytes` are **bytes**, as `tracemalloc`
reports them.  The documentation prints them in **MB = 10^6 bytes**, which is what
`benchmark_memory.py` itself prints:

| tips | k | array, bytes | array MB (10^6) | array MiB (2^20) | copies, bytes | copies MB (10^6) | copies MiB (2^20) |
|---|---|---|---|---|---|---|---|
| 2000 | 25 | 1933888 | 1.93 | 1.84 | 2808112 | 2.81 | 2.68 |
| 2000 | 50 | 3492618 | 3.49 | 3.33 | 5616224 | 5.62 | 5.36 |
| 2000 | 100 | 6691010 | 6.69 | 6.38 | 11232672 | 11.23 | 10.71 |

Read on a MiB basis the same rows are 1.84 / 3.33 / 6.38 and 2.68 / 5.36 / 10.71;
the published 1.93 / 3.49 / 6.69 and 2.81 / 5.62 / 11.23 are the 10^6 reading.
`copies_over_array` (1.45 / 1.61 / 1.68) is unit-free.

**`copies_peak_bytes` is an emulation, not a measurement of R.**  It is the
Python-side structural cost of the reference design - `k` copies of the tree's
edge and length arrays - traced inside this package's process.  There is no R
memory column in this directory, and none was taken; `benchmark_memory.py`'s
docstring says so.  The comparison these rows support is "one (E, k) array vs k
copies inside one Python process", i.e. structural, not cross-language.

## Reproducing a table

```sh
python make_scaling_inputs.py && python benchmark_scaling.py
Rscript benchmark_scaling_r.R validation/reference benchmarks/large benchmarks/scaling \
  benchmarks/results_scaling_r.csv
python benchmark_phyloslicer.py                     # validation grid, Python side
Rscript benchmark_r.R validation/reference benchmarks/results_r.csv
python make_large_inputs.py    && python benchmark_large.py
Rscript benchmark_large_r.R benchmarks/large benchmarks/results_r_large.csv
python benchmark_memory.py                          # Python, peak memory
python recompute_ratios.py                          # check every documented ratio
```

`recompute_ratios.py` needs no R and no timing run: it reads the archived tables,
recomputes each range the documentation quotes (76-556, 148-302 at 100 slices,
6.7-8.1 for the rate fits, and the quartile spread of the widest endpoint), prints
the table, the scenario and the statistic behind every cell, and exits nonzero if
a quoted range no longer follows from the tables, if a row's scenario tag is
missing or `unknown`, if a table's `reps` is not the documented count, or if the
two tables disagree about which cells exist.

## Caveats that must travel with these numbers

1. **The R columns are archived outputs.**  They were produced in an earlier
   session with treesliceR 1.1.0 under R 4.5.3; the R interquartile range is
   stored next to each median in `results_scaling_r.csv`.  Re-runs of the R side
   varied by up to a factor of two between sessions on the same machine, which
   is why every published ratio is quoted to at most two significant figures.
2. **Session drift affects both engines, not just treesliceR.**  The same
   `E x k` contribution-matrix call at 2000 tips and k = 100 is archived twice:
   **Python 0.95158 ms** in `results_scaling_py.csv` and **1.4513 ms** in
   `results_py_large.csv`, while R reads 140.98 ms and 156.17 ms in the matching
   tables.  That gives 148.2x from the scaling pair and 107.6x from the
   end-to-end pair (at 1000 tips, 176.6x vs 104.7x).  Both numbers are real
   measurements of the same operation taken in different sessions; neither table
   is a correction of the other, which is why each published ratio names the
   table it came from.
3. **`tips = 100` has two input trees** (`bd100`, a birth-death tree, and
   `coal100`, a coalescent tree), so `results_scaling_r.csv` carries 21 rows for
   18 (tree size, k) cells.  Ratios are computed against the `bd100` row,
   matching the Python side, and the selection is by the `source` tag rather than
   by row order.  The choice does not move the quoted ranges - both pairings give
   76-556 overall and 148-302 at k = 100 - but it does move the rung itself:
   the 100-tip ratios are 93.8 / 148.9 / 210.8x under `bd100` and 92.1 / 144.9 /
   209.4x if the Python `bd100` row is paired against the R `coal100` row, and
   the Spearman correlation of ratio against tree size at k = 50 is +0.43 versus
   +0.54.  `recompute_ratios.py` prints both readings, uses only the declared
   one, and fails if a table loses the tag that makes them distinguishable.
4. **No BLAS/thread pinning.**  Only `validation/compare.py` pins
   `OMP/OPENBLAS/MKL/VECLIB` threads; these benchmark scripts do not, so
   sub-millisecond medians move in their last digits between runs.
5. **Peak memory is measured on the Python side only**, and the copy side of that
   comparison is an emulation of the reference design rather than a measurement of
   the R process (see "Units" above).
6. **The raw repeat vectors are only half archived.**  `benchmark_scaling.py`,
   `benchmark_scaling_r.R`, `benchmark_phyloslicer.py` and `benchmark_r.R` now
   write every individual run (`repeat_times_s`, or `repeat_slice_s` /
   `repeat_cpd_s` / `repeat_cpe_s` where a row carries several statistics), so a
   median can be recomputed from its inputs.  The archived tables predate those
   columns: they store median, `q25`, `q75`, minimum and `reps` per cell, but no
   repeat vector.  A documented ratio is therefore auditable to its quartiles, not
   to its individual repeats; that is the residual limit on the endpoint claims,
   and `recompute_ratios.py` states it in its own output.  No median here was
   re-derived from repeats that were not archived.
7. **`results_r_large.csv` has `cpe_s = NA`** and `results_py_large.csv` has no
   `cpe_s` column, so no CpE comparison exists at large tree sizes.

## Machine

The archived tables were produced on an Apple silicon laptop.  The Python
environment recorded for the shipped validation artefacts
(`validation/reference/manifest.txt`) is Python 3.14.6 with
numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3; the R side used R 4.5.3
(aarch64-apple-darwin20.0.0) with ape 5.8-1 and treesliceR 1.1.0.
