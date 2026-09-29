#!/usr/bin/env Rscript
# Time treesliceR on the validation scenarios.
# Reads the exact tree/mat inputs exported by generate_reference.R so that
# the R and Python timings refer to byte-identical inputs.
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
rows <- list()
ri <- 1

for (tf in tree_files) {
  tag <- sub("__tree\\.tre$", "", basename(tf))
  tree <- read.tree(tf)
  for (mf in list.files(ref_dir, pattern = paste0("^", tag, "__.*__mat\\.csv$"), full.names = TRUE)) {
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
                               slice_s = t_slice, cpd_s = t_cpd, cpe_s = t_cpe)
      ri <- ri + 1
      cat(sprintf("%s tips=%d sites=%d slices=%d: slice %.3fs cpd %.3fs cpe %.3fs\n",
                  tag, Ntip(tree), nrow(mat), n_slices, t_slice, t_cpd, t_cpe))
    }
  }
}

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
cat("wrote", out_csv, "\n")
