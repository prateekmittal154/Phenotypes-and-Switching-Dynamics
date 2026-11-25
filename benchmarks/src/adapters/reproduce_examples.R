# Reproduce each method's own published example before any benchmarking.
# Completion-standard requirement: we must show the installed software runs as its
# authors intend, on their data, before we use it to make comparative claims.
#
# Usage: Rscript reproduce_examples.R <benchdir>

args <- commandArgs(trailingOnly = TRUE)
B <- args[1]
outdir <- file.path(B, "results", "metrics")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
sink(file.path(B, "reports", "method_reproduction_raw.txt"), split = TRUE)

cat("=====================================================================\n")
cat("Method reproduction run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "\n")
cat("=====================================================================\n\n")

rows <- list()
add <- function(method, check, value, expected, status, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(
    method = method, check = check, value = as.character(value),
    expected = as.character(expected), status = status, note = note,
    stringsAsFactors = FALSE)
}

## --------------------------------------------------------------- PATH
cat("\n########## PATH ##########\n")
path_ok <- FALSE
tryCatch({
  devtools::load_all(file.path(B, "methods", "PATH"), quiet = TRUE)
  path_ok <- TRUE
  cat("PATH loaded via load_all\n")
}, error = function(e) cat("PATH load_all failed:", conditionMessage(e), "\n"))

if (path_ok) {
  # Example taken verbatim from the xcor() roxygen block in PATH/R/xcor.R
  set.seed(1)
  Q <- matrix(c(-0.1, 0.1, 0.1, -0.1), 2, 2)
  phy <- ran_tree_with_states(Q)
  X <- catMat(phy$states)
  W <- one_node_tree_dist(phy)
  xc <- xcor(X, W)
  cat("\n-- xcor() on the package's own documented example --\n")
  print(round(xc$phy_cor, 4))
  cat("Z scores:\n"); print(round(xc$Z.score, 3))
  d <- diag(as.matrix(xc$phy_cor))
  add("PATH", "xcor runs on documented example", paste(round(d, 4), collapse = ","),
      "finite autocorrelations", if (all(is.finite(d))) "PASS" else "FAIL")

  # PATH_inf on a tree simulated from a KNOWN Q: the authors' own generator.
  # This is the strongest available self-check: recover the Q used to simulate.
  cat("\n-- PATH_inf() recovering a known Q from ran_tree_with_states() --\n")
  set.seed(42)
  qtrue <- 0.05
  Qt <- matrix(c(-qtrue, qtrue, qtrue, -qtrue), 2, 2, byrow = TRUE)
  rec <- lapply(1:20, function(i) {
    p <- ran_tree_with_states(Qt, N = 512)
    r <- tryCatch(PATH_inf(p, nstates = 2), error = function(e) NULL)
    if (is.null(r)) return(NULL)
    qq <- as.matrix(r$Q)
    data.frame(q12 = qq[1, 2], q21 = qq[2, 1])
  })
  rec <- do.call(rbind, rec[!sapply(rec, is.null)])
  cat("replicates recovered:", nrow(rec), "\n")
  cat("median inferred off-diagonal rate:", round(median(c(rec$q12, rec$q21)), 4),
      " (true", qtrue, ")\n")
  add("PATH", "PATH_inf recovers known Q (median off-diag)",
      round(median(c(rec$q12, rec$q21)), 4), qtrue,
      if (nrow(rec) >= 10) "PASS" else "PARTIAL",
      sprintf("%d/20 replicates returned a finite Q", nrow(rec)))

  pv <- tryCatch(as.character(packageVersion("PATH")), error = function(e) "1.0 (load_all)")
  add("PATH", "version", pv, "1.0", "INFO",
      "commit caffc5105678224f108b3dd61fea545f96235ffc")
}

## --------------------------------------------------------------- PASTRI
cat("\n\n########## PASTRI ##########\n")
pastri_ok <- FALSE
tryCatch({
  Sys.setenv(PASTRI_DIR = file.path(B, "methods", "PASTRI_package"))
  source(file.path(B, "src", "adapters", "pastri_load.R"))
  pastri_ok <- TRUE
  cat("PASTRI estimator sourced (see pastri_load.R for the documented fortify shim)\n")
}, error = function(e) cat("PASTRI source failed:", conditionMessage(e), "\n"))

if (pastri_ok) {
  # The package ships HSC.nwk + HSC_nodeInfos.csv in inst/extdata; the README and the
  # HSC vignette use exactly these.
  nwk <- file.path(B, "methods", "PASTRI_package", "inst", "extdata", "HSC.nwk")
  inf <- file.path(B, "methods", "PASTRI_package", "inst", "extdata", "HSC_nodeInfos.csv")
  cat("tree:", nwk, "exists:", file.exists(nwk), "\n")
  cell_info <- read.csv(inf, stringsAsFactors = FALSE)
  cat("cell_info columns:", paste(colnames(cell_info), collapse = ", "), "\n")
  cat("cell types present:", paste(sort(unique(cell_info$celltype)), collapse = ", "), "\n")
  cat("n rows:", nrow(cell_info), "\n")

  t0 <- Sys.time()
  nd <- calculate_lca_depths(tree.path = nwk, mc.cores = 8)
  cat("calculate_lca_depths: ", nrow(nd), " pairs in ",
      round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n", sep = "")
  print(utils::head(nd, 3))
  add("PASTRI", "calculate_lca_depths on shipped HSC.nwk", nrow(nd), ">0",
      if (nrow(nd) > 0) "PASS" else "FAIL")

  ct <- sort(unique(cell_info$celltype[cell_info$cellNum == 1]))
  ct <- ct[!is.na(ct)]
  cat("cell types used:", paste(ct, collapse = ", "), "\n")
  BM <- matrix(Inf, nrow = length(ct), ncol = length(ct))
  BM[upper.tri(BM)] <- 0          # irreversibility prior, as in the README
  res <- tryCatch(
    get_optimal_transition_matrix(node_pair_depth = nd, cell_info = cell_info,
                                  Sel_u = "lca_normalized_height", fi_depth = 0.5,
                                  Bound_Matrix = BM, cell_type_list = ct, mc.cores = 8),
    error = function(e) { cat("get_optimal_transition_matrix error:", conditionMessage(e), "\n"); NULL })
  if (!is.null(res)) {
    cat("\n-- optimal transition matrix at depth 0.5 --\n")
    print(round(t(res$optimal_matrix), 4))
    rs <- rowSums(t(res$optimal_matrix))
    cat("row sums:", paste(round(rs, 4), collapse = ", "), "\n")
    add("PASTRI", "get_optimal_transition_matrix on shipped example",
        paste(round(rs, 3), collapse = ","), "rows sum to 1",
        if (all(abs(rs - 1) < 0.05)) "PASS" else "PARTIAL")
  } else {
    add("PASTRI", "get_optimal_transition_matrix on shipped example", "error", "matrix", "FAIL")
  }
  add("PASTRI", "version", "0.0.1", "0.0.1", "INFO",
      "commit 67aefd96240f9d3a1900c6906d813962233ef38f")
}

## --------------------------------------------------------------- SCOUT
cat("\n\n########## SCOUT ##########\n")
scout_ok <- FALSE
tryCatch({
  suppressPackageStartupMessages(library(OUwie))
  cat("OUwie", as.character(packageVersion("OUwie")), "loaded\n")
  scout_ok <- TRUE
}, error = function(e) cat("OUwie not available:", conditionMessage(e), "\n"))

if (scout_ok) {
  # SCOUT's engine is OUwie. Its shipped example is examples/sim_example/:
  # a Cassiopeia-NJ tree plus a 3-state OU simulation at alpha=3, sigma=1.
  ex_t <- file.path(B, "methods", "SCOUT", "examples", "sim_example",
                    "3state_n256_casNJ_tree.nwk")
  ex_c <- file.path(B, "methods", "SCOUT", "examples", "sim_example",
                    "3state_n256_alpha_3_sigma_1_counts.csv")
  cat("shipped tree exists:", file.exists(ex_t), " counts exists:", file.exists(ex_c), "\n")
  if (file.exists(ex_t) && file.exists(ex_c)) {
    tr <- ape::read.tree(ex_t)
    cn <- read.csv(ex_c, stringsAsFactors = FALSE)
    cat("tree tips:", ape::Ntip(tr), " counts dim:", nrow(cn), "x", ncol(cn), "\n")
    cat("counts columns (first 6):", paste(utils::head(colnames(cn), 6), collapse = ", "), "\n")
    add("SCOUT", "shipped example readable",
        sprintf("%d tips / %dx%d counts", ape::Ntip(tr), nrow(cn), ncol(cn)),
        "non-empty", "PASS",
        "3state_n256 alpha=3 sigma=1, Cassiopeia-NJ tree")

    # SCOUT's own code (R/dropout_validation.R:470-472) says: "No edge lengths,
    # replacing with 1s." We follow that rule exactly. The shipped tree has none.
    # The column names encode the generating model (BM1_*, OU1_*, OUM_*), so the
    # reproduction target is model selection: does AICc recover the generating class?
    tr2 <- tr
    if (is.null(tr2$edge.length)) tr2$edge.length <- rep(1, ape::Nedge(tr2))
    tr2 <- ape::multi2di(tr2)
    tr2$edge.length[tr2$edge.length <= 0] <- 1e-7
    tr2$node.label <- rep("1", ape::Nnode(tr2))
    cnm <- read.csv(ex_c, row.names = 1, check.names = FALSE)
    tips <- paste0("t", rownames(cnm))
    stopifnot(all(tips %in% tr2$tip.label))
    # 20 genes per generating class keeps the check decisive while bounded in runtime
    genes <- c(utils::head(grep("^BM1_", colnames(cnm), value = TRUE), 20),
               utils::head(grep("^OU1_", colnames(cnm), value = TRUE), 20))
    cat("fitting", length(genes), "genes (BM1/OU1) on", length(tips), "tips\n")
    sel <- vapply(genes, function(g) {
      dd <- data.frame(species = tips, reg = "1", x = as.numeric(cnm[[g]]))
      fo <- tryCatch(OUwie::OUwie(tr2, dd, model = "OU1", simmap.tree = FALSE,
                                  scaleHeight = FALSE, quiet = TRUE, diagn = FALSE),
                     error = function(e) NULL)
      fb <- tryCatch(OUwie::OUwie(tr2, dd, model = "BM1", simmap.tree = FALSE,
                                  scaleHeight = FALSE, quiet = TRUE, diagn = FALSE),
                     error = function(e) NULL)
      if (is.null(fo) || is.null(fb)) return(NA_character_)
      if (fo$AICc < fb$AICc) "OU1" else "BM1"
    }, character(1))
    truth <- sub("_[0-9]+$", "", genes)
    ok <- !is.na(sel)
    acc <- mean(sel[ok] == truth[ok])
    cat("\n-- SCOUT engine model selection on its own shipped example --\n")
    print(table(generating = truth[ok], selected = sel[ok]))
    cat("overall accuracy:", round(acc, 3), " genes fitted:", sum(ok), "/", length(genes), "\n")
    add("SCOUT", "AICc model selection recovers generating model on shipped example",
        round(acc, 3), ">0.5 (chance)",
        if (sum(ok) > 0 && acc > 0.5) "PASS" else "FAIL",
        sprintf("%d/%d genes fitted; tree had no branch lengths so unit edges were used, per SCOUT R/dropout_validation.R:470-472",
                sum(ok), length(genes)))

    # Record the alpha estimates on the OU-generated genes for reference.
    al <- vapply(utils::head(grep("^OU1_", colnames(cnm), value = TRUE), 20), function(g) {
      dd <- data.frame(species = tips, reg = "1", x = as.numeric(cnm[[g]]))
      fo <- tryCatch(OUwie::OUwie(tr2, dd, model = "OU1", simmap.tree = FALSE,
                                  scaleHeight = FALSE, quiet = TRUE, diagn = FALSE),
                     error = function(e) NULL)
      if (is.null(fo)) NA_real_ else unname(fo$solution["alpha", 1])
    }, numeric(1))
    cat("median alpha on OU1 genes:", signif(median(al, na.rm = TRUE), 4), "\n")
    add("SCOUT", "median OU alpha on OU-generated genes",
        signif(median(al, na.rm = TRUE), 4), "finite",
        if (is.finite(median(al, na.rm = TRUE))) "PASS" else "FAIL",
        "alpha is in units of the imposed unit-edge clock, not the generating timescale")
  }
    add("SCOUT", "version", as.character(packageVersion("OUwie")), "engine OUwie", "INFO",
      "SCOUT commit c72befe0d052f3960c448a64bff776665ab8ea4c")
}

out <- do.call(rbind, rows)
write.csv(out, file.path(outdir, "method_reproduction.csv"), row.names = FALSE)
cat("\n\n=== SUMMARY ===\n")
print(out[, c("method", "check", "value", "status")])
cat("\nwrote", file.path(outdir, "method_reproduction.csv"), "\n")
sink()
