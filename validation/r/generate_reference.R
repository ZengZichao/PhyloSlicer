#!/usr/bin/env Rscript
# Generate treesliceR reference values for the cross-validation protocol
# (cross-validation step 4).  treesliceR is treated strictly as a black box.
#
# Usage: Rscript generate_reference.R <out_dir>
# Outputs per-scenario CSV files with CpD / CpE / r_phylo(PD|PE) reference values
# and the beta-diversity set - the exported queen adjacency plus CpB, CpB_RW and
# r_phylo(index = "PB") - plus the criterion = "PD" variants queued for the
# equal-PD criterion.
#
# RUN LOG.  R 4.5.3 (aarch64-apple-darwin20.0.0) with ape 5.8.1 and treesliceR
# 1.1.0, 2026-09-21.  Three facts came out of that run and all are load-bearing:
#
#   1. Every file this script writes for the *time* criterion reproduced the
#      committed reference byte-for-byte - the three Newick trees, the six
#      presence matrices and all twelve CpD/CpE CSVs (21 files, `cmp` clean),
#      which pins the seed, the RNG and the package versions.  The twelve
#      *__rPD.csv / *__rPE.csv files did NOT: the committed ones are
#      slice-by-site with 30 columns, this script writes site-by-slice with 50,
#      so the committed pair came from an older revision of this file.  That
#      settles the profile provenance, and the regenerated profiles - whose companion
#      inputs are provably the registered ones - are what
#      validation/reference/manifest.txt now registers.
#   2. The criterion = "PD" block below failed on every scenario (18 warnings,
#      0 files written), so no *_PD.csv exists, none can exist for this
#      treesliceR version as spelled in its own documentation, and the equal-PD
#      criterion stays an internal-invariant claim rather than a
#      cross-implementation one - now for a demonstrated reason.  CORRECTED BY
#      THE 2026-09-21 RE-RUN, which prints the message where it happens
#      (options(warn = 1) below): every criterion = "PD" call fails with
#        Error: object 'cutted_tree' not found
#      not with the "length of 'dimnames' [2] not equal to array extent" this log
#      used to record.  The cause is visible in phylo_pieces.R:42 / :152 - the
#      function branches on criterion == "my" and criterion == "pd", while CpD /
#      CpE / r_phylo document and forward "PD", so neither branch runs and
#      `cutted_tree` is never assigned.  Verified directly: r_phylo(criterion =
#      "PD") errors with 'cutted_tree' not found while r_phylo(criterion = "pd")
#      returns.  Equal-PD slicing is therefore *reachable* through the
#      undocumented lower-case spelling, which is a separate question from the
#      structural divergence described under "equal-PD slicing" below; this
#      script still only attempts the documented spelling, and nothing about
#      criterion = "PD" is claimed numerically here.
#   3. The beta block below reaches treesliceR, which is new: it used to be that
#      no CpB reference was ever generated, and every beta statement could only
#      be an internal invariant.  Per scenario it writes 1 adjacency input, 5 of
#      the 6 CpB(component, method) files, both CpB_RW files and 4 of the 6
#      r_phylo(index = "PB") profiles - 72 files over the six scenarios.  Three
#      families produce nothing, because the reference itself cannot produce
#      them (measured warnings, identical on all six scenarios):
#        CpB(comp = "sorensen", method = "pairwise") ->
#          "task 1 failed - 'qr' and 'y' must have the same number of rows",
#          because that branch holds `pbd` as a length-P vector of per-pair
#          Sorensen values (CpB.R:479) and recycles it against the length-k slice
#          vector (CpB.R:525) before nls ever sees it;
#        r_phylo(index = "PB", method = "pairwise") for sorensen and turnover ->
#          a ragged list whose element lengths are 300 / 750 / 1800 (= 50 slices
#          x the 6 / 15 / 36 pairs of a 6x5 queen neighbourhood), i.e. a
#          flattened k x P block rather than one curve per assemblage.  Its
#          nestedness branch does collapse to length 50 and is registered.
#      write_beta_rates() / write_beta_profile() refuse any other shape rather
#      than forcing one, so a family that treesliceR gets wrong stays absent
#      from the reference set instead of becoming a reference.
#
# Re-run with:
#   micromamba run -p ~/.local/share/mamba/envs/r-4.5.3 \
#     Rscript validation/r/generate_reference.R <out_dir>

suppressPackageStartupMessages({
  library(ape)
  library(treesliceR)
})

# Print every warning where it happens.  With the default (warn = 0) a redirected
# run only reports "There were N warnings" at the end and the text - which for the
# beta block below *is* the finding - is lost, so the guard clauses here would
# fail silently in exactly the situation they exist to document.
options(warn = 1L)

out_dir <- if (length(commandArgs(TRUE)) >= 1) commandArgs(TRUE)[1] else "reference"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260912)

# ---- scenario trees (20 / 100 tips; mixed birth-death + posterior-like) ----
make_bd <- function(n, b, d, tmax = 1) {
  repeat {
    tr <- rbdtree(b, d, Tmax = tmax)
    if (Ntip(tr) == n) {
      tr$edge.length <- tr$edge.length / max(node.depth.edgelength(tr))
      return(tr)
    }
  }
}

make_trees <- function() {
  trees <- list()
  trees$bd20 <- make_bd(20, 3.2, 0.2, tmax = 1)   # E[n] ~ exp(3.0)                   # birth-death, 20 tips
  trees$bd100 <- make_bd(100, 4.8, 0.2, tmax = 1)  # E[n] ~ exp(4.6)                 # birth-death, 100 tips
  # posterior-like: coalescent tree re-scaled to unit age
  t100 <- rcoal(100); t100$edge.length <- t100$edge.length / max(node.depth.edgelength(t100))
  trees$coal100 <- t100
  trees
}

trees <- make_trees()

# ---- presence matrices: random + clustered structure ----
make_mats <- function(tree, n_sites) {
  S <- length(tree$tip.label)
  mats <- list()
  random <- matrix(sample(c(0, 1), n_sites * S, replace = TRUE,
                          prob = c(0.7, 0.3)), nrow = n_sites)
  # guarantee every species occurs
  for (s in which(colSums(random) == 0)) random[sample(n_sites, 1), s] <- 1
  colnames(random) <- tree$tip.label
  for (r in which(rowSums(random) == 0)) random[r, sample(S, 1)] <- 1
  mats$random <- random

  clustered <- matrix(0, nrow = n_sites, ncol = S)
  side <- ceiling(sqrt(n_sites))
  coords <- cbind(rep(1:side, times = side), rep(1:side, each = side))[1:n_sites, , drop = FALSE]
  for (s in seq_len(S)) {
    cx <- runif(1, 1, side); cy <- runif(1, 1, side); r <- runif(1, 1, side / 2)
    inside <- (coords[, 1] - cx)^2 + (coords[, 2] - cy)^2 <= r^2
    if (!any(inside)) inside[sample(n_sites, 1)] <- TRUE
    clustered[inside, s] <- 1
  }
  colnames(clustered) <- tree$tip.label
  for (r in which(rowSums(clustered) == 0)) clustered[r, sample(S, 1)] <- 1
  mats$clustered <- clustered
  mats
}

# ---- queen adjacency on a small grid (matches treesliceR semantics) ----
queen_adj <- function(n_sites) {
  side <- ceiling(sqrt(n_sites))
  cells <- cbind(rep(1:side, times = side), rep(1:side, each = side))[1:n_sites, , drop = FALSE]
  adj <- matrix(0, n_sites, n_sites)
  for (i in seq_len(n_sites)) {
    dx <- abs(cells[, 1] - cells[i, 1])
    dy <- abs(cells[, 2] - cells[i, 2])
    nb <- which(dx <= 1 & dy <= 1)
    adj[i, nb] <- 1
  }
  adj
}

n_slices <- 50L

# ---- beta-diversity references (CpB / CpB_RW / r_phylo(index = "PB")) -------
# treesliceR's own names for the two arguments PhyloSlicer spells `component`
# and `approach` are `comp` and `method` (CpB.R:51, CpB_RW.R:50, r_phylo.R:70),
# and its adjacency argument is `adj` in all three.  Everything is wrapped in
# try(): a family the reference implementation cannot produce is reported below
# and simply writes no file, which validation/compare.py then reports as
# `no reference` rather than as an agreement.
beta_comps <- c("sorensen", "turnover", "nestedness")
beta_methods <- c("multisite", "pairwise")

write_beta_rates <- function(obj, path, tag, want_cols) {
  # A treesliceR CpB / CpB_RW result is a length(asb) x 3 matrix whose column
  # names are set by the function itself (CpB.R:753-758).  Anything else - a
  # vector, a matrix whose columns are the recycled per-pair Sorensen vector, a
  # short block because some neighbourhoods bailed out - is not a registrable
  # reference, so refuse it here instead of forcing it into three columns.
  if (!is.matrix(obj) || ncol(obj) != length(want_cols) ||
      !identical(colnames(obj), want_cols)) {
    warning(paste0(tag, ": unexpected ", want_cols[1], " result (",
                   paste(dim(obj), collapse = "x"), " cols=",
                   paste(colnames(obj), collapse = "/"), "), not written"))
    return(invisible(FALSE))
  }
  write.csv(cbind(site_id = as.integer(rownames(obj)), obj), path, row.names = FALSE)
  invisible(TRUE)
}

write_beta_profile <- function(obj, path, tag) {
  # r_phylo returns one vector of length n per assemblage.  The pairwise beta
  # branches return a flattened k x P block instead (measured: element lengths
  # 300 / 750 / 1800 for the 50-slice, 6-36-pair neighbourhoods of the 30-site
  # grid), which has no site-by-slice shape to register.
  lens <- lengths(obj)
  if (!is.list(obj) || length(obj) == 0 || !all(lens == n_slices)) {
    warning(paste0(tag, ": r_phylo(PB) is ragged (element lengths ",
                   paste(unique(lens), collapse = ","), " of ", n_slices,
                   " slices), not written"))
    return(invisible(FALSE))
  }
  df <- matrix(unlist(obj), nrow = length(obj), byrow = TRUE)
  colnames(df) <- paste0("slice_", seq_len(ncol(df)))
  write.csv(cbind(site_id = seq_len(nrow(df)), df), path, row.names = FALSE)
  invisible(TRUE)
}

for (tname in names(trees)) {
  tree <- trees[[tname]]
  if (!is.ultrametric(tree)) {
    tree$edge.length <- tree$edge.length * (1 / max(node.depth.edgelength(tree)))
  }
  # export the exact inputs so the Python side runs on identical data
  writeLines(write.tree(tree, digits = 12), file.path(out_dir, paste0(tname, "__tree.tre")))
  mats <- make_mats(tree, n_sites = 30L)
  adj <- queen_adj(30L)
  for (mname in names(mats)) {
    mat <- mats[mname][[1]]
    tag <- paste0(tname, "__", mname)
    write.csv(mat, file.path(out_dir, paste0(tag, "__mat.csv")))
    tag <- paste0(tname, "__", mname)

    cpd_ref <- treesliceR::CpD(tree, n = n_slices, mat = mat, criterion = "my")
    cpe_ref <- treesliceR::CpE(tree, n = n_slices, mat = mat, criterion = "my")
    rpd_ref <- treesliceR::r_phylo(tree, n = n_slices, mat = mat, index = "PD")
    rpe_ref <- treesliceR::r_phylo(tree, n = n_slices, mat = mat, index = "PE")

    write.csv(cbind(site_id = rownames(cpd_ref), cpd_ref),
              file.path(out_dir, paste0(tag, "__CpD.csv")), row.names = FALSE)
    write.csv(cbind(site_id = rownames(cpe_ref), cpe_ref),
              file.path(out_dir, paste0(tag, "__CpE.csv")), row.names = FALSE)

    cat(" rPD lengths:", paste(sapply(rpd_ref, length), collapse = ","), "\n")
    rpd_df <- matrix(unlist(rpd_ref), nrow = length(rpd_ref), byrow = TRUE)
    colnames(rpd_df) <- paste0("slice_", seq_len(n_slices))
    write.csv(cbind(site_id = seq_len(nrow(mat)), rpd_df),
              file.path(out_dir, paste0(tag, "__rPD.csv")), row.names = FALSE)

    rpe_df <- matrix(unlist(rpe_ref), nrow = length(rpe_ref), byrow = TRUE)
    colnames(rpe_df) <- paste0("slice_", seq_len(n_slices))
    write.csv(cbind(site_id = seq_len(nrow(mat)), rpe_df),
              file.path(out_dir, paste0(tag, "__rPE.csv")), row.names = FALSE)

    # ---- beta diversity: the exported adjacency, then CpB / CpB_RW / rPB -----
    # The adjacency is written out as an *input*, not left as an R-side detail,
    # so that the Python side is scored on exactly the neighbourhood sets the
    # reference consumed.  queen_adj(30) is the 8-neighbour queen relation
    # (edges and corners) on a 6 x 5 grid of unit cells, self included
    # (`adj[i, nb] <- 1` with nb containing i), which is what
    # phyloslicer.spatial.adjacency_from_grid(method = "queen") computes for the
    # same grid - rectangle bounds sharing an edge or a corner touch, diagonal
    # forced to 1.  tests/test_validation_chain.py::TestBetaAdjacencySemantics
    # asserts the two matrices are element-wise equal, so the coupling is
    # checked rather than asserted in a comment.
    adj_out <- adj
    dimnames(adj_out) <- list(paste0("s", seq_len(nrow(adj_out))),
                              paste0("s", seq_len(ncol(adj_out))))
    write.csv(adj_out, file.path(out_dir, paste0(tag, "__adj.csv")))

    for (bm in beta_methods) {
      for (bc in beta_comps) {
        label <- paste0(tag, "__CpB__", bc, "_", bm)
        cpb_ref <- try(suppressWarnings(suppressMessages(
          treesliceR::CpB(tree, n = n_slices, mat = mat, adj = adj,
                          comp = bc, method = bm, criterion = "my")
        )), silent = TRUE)
        if (inherits(cpb_ref, "try-error")) {
          warning(paste0(label, ": CpB failed - ",
                         conditionMessage(attr(cpb_ref, "condition"))))
        } else if (!write_beta_rates(cpb_ref, file.path(out_dir, paste0(label, ".csv")),
                                     label, c("CpB", "PB", "pBO"))) {
          next
        } else {
          cat(" wrote", label, "\n")
        }
      }
      # r_phylo(index = "PB") returns the curve CpB fits *before* dividing by
      # the whole-tree index (r_phylo.R:463 `return(r_pbd)`, vs CpB.R:536
      # `return(r_pbd/pbd)`), so it is registered as its own profile family.
      # Only the multisite branch returns one vector per assemblage; the
      # pairwise branch is ragged and write_beta_profile() refuses it.
      for (bc in beta_comps) {
        label <- paste0(tag, "__rPB__", bc, "_", bm)
        rpb_ref <- try(suppressWarnings(suppressMessages(
          treesliceR::r_phylo(tree, n = n_slices, mat = mat, adj = adj,
                              index = "PB", comp = bc, method = bm, criterion = "my")
        )), silent = TRUE)
        if (inherits(rpb_ref, "try-error")) {
          warning(paste0(label, ": r_phylo(PB) failed - ",
                         conditionMessage(attr(rpb_ref, "condition"))))
        } else if (!write_beta_profile(rpb_ref, file.path(out_dir, paste0(label, ".csv")),
                                       label)) {
          next
        } else {
          cat(" wrote", label, "\n")
        }
      }
    }
    for (bm in beta_methods) {
      label <- paste0(tag, "__CpB_RW__", bm)
      cpbrw_ref <- try(suppressWarnings(suppressMessages(
        treesliceR::CpB_RW(tree, n = n_slices, mat = mat, adj = adj,
                           method = bm, criterion = "my")
      )), silent = TRUE)
      if (inherits(cpbrw_ref, "try-error")) {
        warning(paste0(label, ": CpB_RW failed - ",
                       conditionMessage(attr(cpbrw_ref, "condition"))))
      } else if (!write_beta_rates(cpbrw_ref, file.path(out_dir, paste0(label, ".csv")),
                                   label, c("CpB_RW", "PB_RW", "pBO"))) {
        next
      } else {
        cat(" wrote", label, "\n")
      }
    }

    # ---- equal-PD slicing: queued, never yet run ---------------------------------
    # EXPECT A MISMATCH on any non-binary input, and that is the point of this
    # block, not a defect to fix.  treesliceR hard-codes the number of active
    # branches as 2, 3, 4, ... (squeeze_root.R:98, phylo_pieces.R:159:
    # `nBranch = 2:(length(unique(nodes$YearBegin)) + 1)`), which presumes a
    # strictly binary tree with all-distinct node depths.  On a multifurcating
    # tree that arithmetic does not even conserve total PD; PhyloSlicer counts
    # the branches actually alive at each depth instead.  Comparing these files
    # is therefore diagnostic documentation of a known divergence, NOT an
    # equivalence test: do not feed them to validation/compare.py's pass-rate
    # logic without reading validation/NOTES.md item 3 first.
    # All trees above are ultrametric binary (rbdtree/rcoal rescaled), so on
    # THIS scenario set the two implementations may well agree; the
    # multifurcating case that actually separates them is not generated here.
    # ANSWERED, 2026-09-21: the script was run under R 4.5.3 / ape 5.8.1 /
    # treesliceR 1.1.0 and the criterion = "PD" call failed on all six
    # scenarios (18 warnings, no *_PD.csv written).  The message is
    # `object 'cutted_tree' not found`, not the `length of 'dimnames' [2] not
    # equal to array extent` an earlier log recorded: phylo_pieces branches on
    # "my" / "pd" (phylo_pieces.R:42, :152) while the documented spelling is
    # "PD", so no slice is ever computed.  Binary ultrametric input is therefore
    # not even reaching the comparison.  Nothing may be claimed about criterion =
    # "PD" numerically, in either direction, from this block; whether the
    # undocumented lower-case spelling should be used instead is a separate
    # decision, and the structural warning above still applies to it.
    cpd_pd <- try(treesliceR::CpD(tree, n = n_slices, mat = mat, criterion = "PD"),
                  silent = TRUE)
    if (!inherits(cpd_pd, "try-error")) {
      write.csv(cbind(site_id = rownames(cpd_pd), cpd_pd),
                file.path(out_dir, paste0(tag, "__CpD__PD.csv")), row.names = FALSE)
    } else {
      warning(paste0("criterion='PD' CpD failed for ", tag, ": ", conditionMessage(attr(cpd_pd, "condition"))))
    }
    cpe_pd <- try(treesliceR::CpE(tree, n = n_slices, mat = mat, criterion = "PD"),
                  silent = TRUE)
    if (!inherits(cpe_pd, "try-error")) {
      write.csv(cbind(site_id = rownames(cpe_pd), cpe_pd),
                file.path(out_dir, paste0(tag, "__CpE__PD.csv")), row.names = FALSE)
    } else {
      warning(paste0("criterion='PD' CpE failed for ", tag, ": ", conditionMessage(attr(cpe_pd, "condition"))))
    }
    rpd_pd <- try(treesliceR::r_phylo(tree, n = n_slices, mat = mat,
                                      index = "PD", criterion = "PD"), silent = TRUE)
    if (!inherits(rpd_pd, "try-error")) {
      # NOTE: written site-by-slice with n_slices columns, matching the driver
      # above (NOT the slice-by-site, 30-column orientation of the committed
      # *__rPD.csv files - see validation/NOTES.md item 4).
      rpd_pd_df <- matrix(unlist(rpd_pd), nrow = length(rpd_pd), byrow = TRUE)
      colnames(rpd_pd_df) <- paste0("slice_", seq_len(ncol(rpd_pd_df)))
      write.csv(cbind(site_id = seq_len(nrow(mat)), rpd_pd_df),
                file.path(out_dir, paste0(tag, "__rPD__PD.csv")), row.names = FALSE)
    } else {
      warning(paste0("criterion='PD' r_phylo(PD) failed for ", tag))
    }
    cat("wrote reference:", tag, "\n")
  }
}

cat("done -> ", normalizePath(out_dir), "\n")
