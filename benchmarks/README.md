# Benchmark tables in this directory

Four independent protocols live here.  They are **not** interchangeable, and two
column names mean different things in different files.  Every documented speed-up
is traceable to one of the tables below.

## Protocols

| files | what is timed | statistic | repeats | cross-language? |
|---|---|---|---|---|
| `results.csv` | `bench_core.py`: vectorised kernel vs the naive in-process baseline | **min** | 3 (naive baselines: 1) | **no** - Python vs Python |
| `results_py.csv`, `results_r.csv` | slicing kernel + CpD/CpE fits on the validation-scenario inputs | median | 5 | yes |
| `results_scaling_py.csv`, `results_scaling_r.csv` | **slicing kernel only** (`slice_pieces(...).contribution_matrix()` vs `phylo_pieces(...)`), every tree size at each of k = 25/50/100 | median | 15 (+1 discarded warm-up) | yes |
| `results_py_large.csv`, `results_r_large.csv` | end-to-end CpD on 250 assemblages, 1000/2000 tips, k = 100 | median | 5 | yes |

`slice_s` means **the E x k contribution matrix** in `results_py.csv`,
`results_scaling_py.csv` and `results_py_large.csv`, but **window construction
alone** in `testing_suite/benchmarks/benchmark_results.csv`.  Do not compare the
two across directories.

## Reproducing a table

```sh
python make_scaling_inputs.py && python benchmark_scaling.py       # Python, scaling
Rscript benchmark_scaling_r.R                                      # R, scaling
python make_large_inputs.py    && python benchmark_large.py        # Python, end-to-end
Rscript benchmark_large_r.R                                        # R, end-to-end
python benchmark_memory.py                                         # Python, peak memory
python recompute_ratios.py                                         # check every documented ratio
```

`recompute_ratios.py` needs no R and no timing run: it reads the archived
tables, recomputes each range the documentation quotes (76–556, 148–302 at 100
slices, 6.7–8.1 for the rate fits, and the quartile spread of the widest
endpoint), prints the cells behind it and exits nonzero if a quoted range no
longer follows from the tables.

Both sides read the same `*__tree.tre` inputs, and both writers emit a `source`
column naming the input tag, so a row can be traced to the exact tree it was
measured on.

**The archived `results_scaling_*.csv` tables predate that column.** They are
older output, so their duplicated 100-tip rung is two unlabelled runs of the same
grid point. No published claim depends on which of the two rows is which, and
every ratio is quoted to at most two significant figures for that reason. To
restore provenance, re-run:

```bash
python make_scaling_inputs.py && python benchmark_scaling.py
Rscript benchmark_scaling_r.R
```

## Caveats that must travel with these numbers

1. **The R columns are archived outputs.** They were produced in an earlier
   session with treesliceR 1.1.0 under R 4.5.3; the R interquartile range is
   stored next to each median in `results_scaling_r.csv`.  Re-runs of the R side
   varied by up to a factor of two between sessions on the same machine, which
   is why every published ratio is quoted to at most two significant figures.
2. **`tips = 100` has two input trees** (`bd100`, a birth-death tree, and
   `coal100`, a coalescent tree).  `results_scaling_r.csv` therefore carries 21
   rows for 18 (tree size, k) cells.  Ratios are computed against the
   `bd100` row, matching the Python side.
3. **No BLAS/thread pinning.** Only `validation/compare.py` pins
   `OMP/OPENBLAS/MKL/VECLIB` threads; these benchmark scripts do not, so
   sub-millisecond medians move in their last digits between runs.
4. **Peak memory is measured on the Python side only.**
   `benchmark_memory.py` records the tracemalloc peak of the `E x k`
   contribution matrix against the same number of serialised `TreeArray`
   copies, per tree size, in `results_memory_py.csv`.  The reference side is not
   measured, because copying the tree once per slice is a property of the
   published R implementation rather than something to time against it; the
   comparison drawn here is array-versus-copies within this package.
5. **The raw repeat vectors are only half archived.**  `benchmark_scaling.py`
   and `benchmark_phyloslicer.py` now write a `repeat_times_s` column holding
   every individual run, so a median can be recomputed from its inputs.  The R
   driver writes median, `q25`, `q75` and minimum per cell and no repeat vector,
   and the archived Python tables predate the new column.  A documented ratio is
   therefore auditable to its quartiles, not to its individual repeats; that is
   the residual limit on the endpoint claims, and `recompute_ratios.py` states it
   in its own output.
6. **`results_r_large.csv` has `cpe_s = NA`** and `results_py_large.csv` has no
   `cpe_s` column, so no CpE comparison exists at large tree sizes.

## Machine

The archived tables were produced on an Apple silicon laptop.  The Python
environment recorded for the shipped validation artefacts
(`validation/reference/manifest.txt`) is Python 3.14.6 with
numpy 2.4.4 / scipy 1.17.1 / pandas 3.0.3; the R side used R 4.5.3
(aarch64-apple-darwin20.0.0) with ape 5.8.1 and treesliceR 1.1.0.
