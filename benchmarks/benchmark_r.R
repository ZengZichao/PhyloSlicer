#!/usr/bin/env Rscript
# Time treesliceR on the validation scenarios.
# Reads the exact tree/mat inputs exported by generate_reference.R so that
# the R and Python timings refer to byte-identical inputs.
#
# Each row names the tree it was measured on (source) and the matrix it was
# measured with (structure): the grid is 3 trees x 2 matrix structures x 2 slice
# counts, and (tips, n_sites, n_slices) alone does not identify a row, because
# bd100 and coal100 both have 100 tips.  Every value is the median of REPS timed
# runs and the runs themselves are archived beside it.
#
# Usage: Rscript benchmark_r.R <ref_dir> <out_csv>

suppressPackageStartupMessages({
  library(ape)
  library(treesliceR)
})

args <- commandArgs(TRUE)
ref_dir <- if (length(args) >= 1) args[1] else "validation/reference"
out_csv <- if (length(args) >= 2) args[2] else "benchmarks/results_r.csv"

tree_files <- list.files(ref_dir, pattern = "__tree\\.tre$", full.names = TRUE)
REPS <- 5L
# NV-08: name the box these rows came from, in the table itself; read once, before any
# timing loop.  Column names match machine_info.py (see benchmarks/machine_info.R).
.mi_self <- grep("^--file=", commandArgs(FALSE), value = TRUE)
.mi_dir <- if (length(.mi_self)) {
  dirname(normalizePath(sub("^--file=", "", .mi_self[[1]])))
} else {
  getwd()
}
source(file.path(.mi_dir, "machine_info.R"))
MI <- machine_fields()
cat(describe_machine(MI), "\n")
rows <- list()
ri <- 1

# one CSV cell holding every timed run, so an archived median can be recomputed
vec_s <- function(x) paste(sprintf("%.9g", x), collapse = ";")

for (tf in tree_files) {
  tag <- sub("__tree\\.tre$", "", basename(tf))
  tree <- read.tree(tf)
  for (mf in list.files(ref_dir, pattern = paste0("^", tag, "__.*__mat\\.csv$"), full.names = TRUE)) {
    # the middle token of <tag>__<structure>__mat.csv: clustered / random
    structure <- sub("__mat\\.csv$", "", basename(mf))
    structure <- sub(paste0("^", tag, "__"), "", structure)
    mat <- read.csv(mf, row.names = 1)
    mat <- as.matrix(mat > 0)
    storage.mode(mat) <- "numeric"
    colnames(mat) <- colnames(read.csv(mf, row.names = 1))

    for (n_slices in c(25L, 50L)) {
      # warm-up, then median of REPS timed runs (the protocol in benchmarks/README.md)
      invisible(suppressMessages(phylo_pieces(tree, n = n_slices, criterion = "my")))
      ts <- replicate(REPS, {
        t0 <- Sys.time()
        invisible(suppressMessages(phylo_pieces(tree, n = n_slices, criterion = "my")))
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      })
      t_slice <- stats::median(ts)

      invisible(suppressMessages(CpD(tree, n = n_slices, mat = mat, criterion = "my")))
      tc <- replicate(REPS, {
        t0 <- Sys.time()
        invisible(suppressMessages(CpD(tree, n = n_slices, mat = mat, criterion = "my")))
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      })
      t_cpd <- stats::median(tc)

      invisible(suppressMessages(CpE(tree, n = n_slices, mat = mat, criterion = "my")))
      te <- replicate(REPS, {
        t0 <- Sys.time()
        invisible(suppressMessages(CpE(tree, n = n_slices, mat = mat, criterion = "my")))
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      })
      t_cpe <- stats::median(te)

      rows[[ri]] <- data.frame(tips = Ntip(tree), n_sites = nrow(mat),
                               n_slices = n_slices,
                               slice_s = t_slice, cpd_s = t_cpd, cpe_s = t_cpe,
                               source = tag, structure = structure, reps = REPS,
                               repeat_slice_s = vec_s(ts), repeat_cpd_s = vec_s(tc),
                               repeat_cpe_s = vec_s(te))
      ri <- ri + 1
      cat(sprintf("%s/%s tips=%d sites=%d slices=%d: slice %.3fs cpd %.3fs cpe %.3fs\n",
                  tag, structure, Ntip(tree), nrow(mat), n_slices, t_slice, t_cpd, t_cpe))
    }
  }
}

out <- do.call(rbind, rows)
for (nm in names(MI)) out[[nm]] <- MI[[nm]]  # constant machine columns (NV-08)
write.csv(out, out_csv, row.names = FALSE)
cat("wrote", out_csv, "\n")
