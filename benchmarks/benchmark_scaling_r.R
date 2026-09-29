#!/usr/bin/env Rscript
# Fixed-k scaling benchmark, R side: treesliceR slicing-kernel time for the
# SAME number of slices across tree sizes.
#
# The earlier scripts (benchmark_r.R for 20/100 tips at k = 25/50 and
# benchmark_large_r.R for 1000/2000 tips at k = 100) varied tree size and k
# together, so "the speed-up grows with tree size" could not be separated
# from "the speed-up grows with the number of slices".  Here every tree size
# is timed at every k in {25, 50, 100}, median of five runs.
#
# Usage: Rscript benchmark_scaling_r.R <input_dirs...> <out_csv>
# Input directories must contain <tag>__tree.tre files.

suppressPackageStartupMessages({
  library(ape)
  library(treesliceR)
})

args <- commandArgs(TRUE)
out_csv <- args[length(args)]
in_dirs <- args[-length(args)]
REPS <- 15L  # median of REPS timed runs; warm-up run excluded.
# R timings on this laptop varied by up to a factor of two between sessions, so the
# reported median is accompanied by the interquartile range (see benchmarks/README.md).

tree_files <- unlist(lapply(in_dirs, function(d)
  list.files(d, pattern = "__tree\\.tre$", full.names = TRUE)))

rows <- list()
ri <- 1
for (tf in sort(tree_files)) {
  tag <- sub("__tree\\.tre$", "", basename(tf))
  tree <- suppressWarnings(read.tree(tf))
  for (k in c(25L, 50L, 100L)) {
    times <- numeric(REPS)
    invisible(suppressMessages(phylo_pieces(tree, n = k, criterion = "my")))
    for (rep in 1:REPS) {
      t0 <- Sys.time()
      invisible(suppressMessages(phylo_pieces(tree, n = k, criterion = "my")))
      times[rep] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    }
    qs <- stats::quantile(times, c(0.25, 0.75))
    rows[[ri]] <- data.frame(tips = Ntip(tree), n_slices = k, source = tag,
                             slice_s = stats::median(times),
                             slice_q25 = qs[[1]], slice_q75 = qs[[2]],
                             slice_min = min(times), reps = REPS)
    ri <- ri + 1
    cat(sprintf("%s tips=%d k=%d: slice median %.4fs (IQR %.4f-%.4f, min %.4f)\n",
                tag, Ntip(tree), k, stats::median(times), qs[[1]], qs[[2]], min(times)))
  }
}

out <- do.call(rbind, rows)
write.csv(out, out_csv, row.names = FALSE)
cat("wrote", out_csv, "\n")
