#!/usr/bin/env Rscript
# Time treesliceR on the large benchmark inputs (supplement to benchmark_r.R).
# Usage: Rscript benchmark_large_r.R <large_dir> <out_csv>

suppressPackageStartupMessages({
  library(ape)
  library(treesliceR)
})

args <- commandArgs(TRUE)
in_dir <- if (length(args) >= 1) args[1] else "benchmarks/large"
out_csv <- if (length(args) >= 2) args[2] else "benchmarks/results_r_large.csv"
REPS <- 5L

rows <- list()
ri <- 1
for (tf in list.files(in_dir, pattern = "__tree\\.tre$", full.names = TRUE)) {
  tag <- sub("__tree\\.tre$", "", basename(tf))
  tree <- read.tree(tf)
  for (mf in list.files(in_dir, pattern = paste0("^", tag, "__mat\\.csv$"), full.names = TRUE)) {
    mat <- read.csv(mf, row.names = 1)
    hdr <- read.csv(mf, nrow = 1)
    colnames(mat) <- colnames(hdr)[-1]
    mat <- as.matrix(mat > 0)
    storage.mode(mat) <- "numeric"
    for (n_slices in c(100L)) {
      # warm-up, then median of REPS timed runs (the protocol in benchmarks/README.md)
      ts <- replicate(REPS, {
        t0 <- Sys.time()
        invisible(suppressMessages(phylo_pieces(tree, n = n_slices, criterion = "my")))
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      })
      t_slice <- stats::median(ts)
      tc <- replicate(REPS, {
        t0 <- Sys.time()
        invisible(suppressMessages(CpD(tree, n = n_slices, mat = mat, criterion = "my")))
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      })
      t_cpd <- stats::median(tc)
      rows[[ri]] <- data.frame(tips = Ntip(tree), n_sites = nrow(mat),
                               n_slices = n_slices,
                               slice_s = t_slice, cpd_s = t_cpd, cpe_s = NA)
      ri <- ri + 1
      cat(sprintf("%s tips=%d sites=%d k=%d: slice %.2fs cpd %.2fs\n",
                  tag, Ntip(tree), nrow(mat), n_slices, t_slice, t_cpd))
    }
  }
}
write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
cat("wrote", out_csv, "\n")
