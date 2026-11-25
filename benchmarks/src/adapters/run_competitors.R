# Run PATH, PASTRI and SCOUT (OUwie engine) on a set of simulated replicates.
#
# Each replicate directory holds tree_true.nwk (branch lengths = generations),
# tree_unit.nwk (all edges 1) and states.csv (cell, clone, state, trait).
#
# Every method is run PER CLONE, because the quantities being compared -- across-clone
# variance in occupancy and in switching rate -- are defined across clones.
#
# Usage:
#   Rscript run_competitors.R <benchdir> <repdir_list_file> <tree_variant> <out_csv> [cores]

suppressPackageStartupMessages({ library(ape); library(dplyr) })

args  <- commandArgs(trailingOnly = TRUE)
B     <- args[1]
LIST  <- args[2]
TREEV <- args[3]            # "true" or "unit" or "recon"
OUT   <- args[4]
CORES <- if (length(args) >= 5) as.integer(args[5]) else 8

## ---------------------------------------------------------------- load methods
devtools::load_all(file.path(B, "methods", "PATH"), quiet = TRUE)
Sys.setenv(PASTRI_DIR = file.path(B, "methods", "PASTRI_package"))
source(file.path(B, "src", "adapters", "pastri_load.R"))
suppressPackageStartupMessages(library(OUwie))

## ---------------------------------------------------------------- helpers
stat_dist <- function(P) {
  # stationary distribution of a 2x2 stochastic matrix, by the standard closed form
  a <- P[1, 2]; b <- P[2, 1]
  if (!is.finite(a) || !is.finite(b) || (a + b) <= 0) return(c(NA_real_, NA_real_))
  c(b / (a + b), a / (a + b))
}

safe <- function(expr, default = NA_real_) tryCatch(expr, error = function(e) default,
                                                    warning = function(w) suppressWarnings(eval(expr)))

run_path_clone <- function(phy, states) {
  # A clone fixed in one state has no transitions to infer; PATH's catMat() then returns a
  # single column and PATH_inf() cannot form a 2x2 matrix. Record it rather than crash.
  if (length(unique(states)) < 2)
    return(list(k = NA_real_, pi = NA_real_, p_switch = NA_real_, note = "monomorphic clone"))
  # PATH_inf expects the state vector on the tree object.
  phy$states <- states
  r <- tryCatch(PATH_inf(phy, cell_states = "states", nstates = 2),
                error = function(e) NULL)
  if (is.null(r)) return(list(k = NA_real_, pi = NA_real_, p_switch = NA_real_,
                              note = "PATH_inf error"))
  Q <- suppressWarnings(as.matrix(r$Q))
  P <- suppressWarnings(as.matrix(r$P))
  if (nrow(Q) < 2 || ncol(Q) < 2 || nrow(P) < 2 || ncol(P) < 2)
    return(list(k = NA_real_, pi = NA_real_, p_switch = NA_real_, note = "Q not 2x2"))
  if (any(!is.finite(Q))) return(list(k = NA_real_, pi = NA_real_, p_switch = NA_real_,
                                      note = "non-finite Q"))
  # Q is already a rate matrix; our k is the TOTAL rate (sum of both off-diagonals), which
  # matches the generator of our two-state chain.
  kk <- Q[1, 2] + Q[2, 1]
  sd_ <- stat_dist(P)
  list(k = as.numeric(kk), pi = as.numeric(sd_[2]),
       p_switch = as.numeric(P[1, 2] + P[2, 1]), note = "")
}

# PASTRI returns a per-cell-cycle transition PROBABILITY matrix T (rows sum to 1), not a
# rate matrix. For a two-state chain with total rate k over one generation, the off-diagonal
# entries of the one-step probability matrix sum to 1 - exp(-k). We therefore convert with
# k = -log(1 - p_switch) so that PASTRI, PATH and Study D are all on the same per-generation
# rate scale. p_switch is also reported raw.
#
# fi_depth is swept rather than fixed. The authors use cumulative windows of 0.4-0.8 in
# their own applications and state that including more distal adjacencies overestimates
# rates, so a single fixed value -- especially 1.0 -- would not be a fair setting.
PASTRI_DEPTHS <- c(0.2, 0.4, 0.6, 0.8, 1.0)
PASTRI_HEADLINE <- 0.4

run_pastri_clone <- function(nwk_path, states_df, cores = 1) {
  blank <- list(k = NA_real_, pi = NA_real_, p_switch = NA_real_, note = "", sweep = NULL)
  nd <- tryCatch(calculate_lca_depths(tree.path = nwk_path, mc.cores = cores),
                 error = function(e) NULL)
  if (is.null(nd) || !nrow(nd)) { blank$note <- "lca_depths error"; return(blank) }
  ci <- data.frame(nodeLabel = states_df$cell, cell_ID = states_df$cell,
                   celltype = paste0("S", states_df$state), cellNum = 1,
                   stringsAsFactors = FALSE)
  ct <- c("S0", "S1")
  if (length(unique(ci$celltype)) < 2) { blank$note <- "monomorphic clone"; return(blank) }
  BM <- matrix(Inf, 2, 2)      # no irreversibility prior: both directions allowed

  sweep <- lapply(PASTRI_DEPTHS, function(fd) {
    res <- tryCatch(
      get_optimal_transition_matrix(node_pair_depth = nd, cell_info = ci,
                                    Sel_u = "lca_normalized_height", fi_depth = fd,
                                    Bound_Matrix = BM, cell_type_list = ct, mc.cores = cores),
      error = function(e) NULL)
    if (is.null(res)) return(data.frame(fi_depth = fd, p_switch = NA_real_,
                                        k = NA_real_, pi = NA_real_))
    Tm <- t(res$optimal_matrix)            # rows = from-state, per the README
    if (any(!is.finite(Tm))) return(data.frame(fi_depth = fd, p_switch = NA_real_,
                                               k = NA_real_, pi = NA_real_))
    ps <- as.numeric(Tm[1, 2] + Tm[2, 1])
    kk <- if (ps >= 1) Inf else -log(1 - ps)
    data.frame(fi_depth = fd, p_switch = ps, k = kk, pi = as.numeric(stat_dist(Tm)[2]))
  })
  sweep <- do.call(rbind, sweep)
  h <- sweep[sweep$fi_depth == PASTRI_HEADLINE, , drop = FALSE]
  list(k = if (nrow(h)) h$k[1] else NA_real_,
       pi = if (nrow(h)) h$pi[1] else NA_real_,
       p_switch = if (nrow(h)) h$p_switch[1] else NA_real_,
       note = "", sweep = sweep)
}

source(file.path(B, "src", "adapters", "scout_fn.R"))

## ---------------------------------------------------------------- main loop
reps <- readLines(LIST)
reps <- reps[nzchar(reps)]
cat("replicates:", length(reps), " tree variant:", TREEV, "\n")

rows <- list()
sweeps <- list()
tstart <- Sys.time()
for (ri in seq_along(reps)) {
  d <- reps[ri]
  nwk <- file.path(d, paste0("tree_", TREEV, ".nwk"))
  if (!file.exists(nwk)) { cat("missing", nwk, "\n"); next }
  tr_full <- read.tree(nwk)
  st <- read.csv(file.path(d, "states.csv"), stringsAsFactors = FALSE)
  clones <- sort(unique(st$clone))
  scen <- basename(dirname(file.path(d, "x")))

  for (cl in clones) {
    sub <- st[st$clone == cl, , drop = FALSE]
    if (nrow(sub) < 6) next
    phy <- tryCatch(keep.tip(tr_full, sub$cell), error = function(e) NULL)
    if (is.null(phy)) next
    sub <- sub[match(phy$tip.label, sub$cell), , drop = FALSE]

    # PATH
    rp <- list(k = NA_real_, pi = NA_real_, p_switch = NA_real_, note = "uncaught error")
    tp <- system.time(try(rp <- run_path_clone(phy, sub$state), silent = TRUE))[["elapsed"]]

    # PASTRI: needs a file path
    tmp <- tempfile(fileext = ".nwk"); write.tree(phy, tmp)
    rq <- list(k = NA_real_, pi = NA_real_, p_switch = NA_real_, note = "uncaught error",
               sweep = NULL)
    tq <- system.time(try(rq <- run_pastri_clone(tmp, sub, cores = 1), silent = TRUE))[["elapsed"]]
    unlink(tmp)

    # SCOUT engine on the continuous trait
    ro <- list(alpha = NA_real_, theta = NA_real_, sigma2 = NA_real_, sel = NA_character_,
               note = "uncaught error")
    ts <- system.time(try(ro <- run_ouwie_clone(phy, sub$trait, sub$cell), silent = TRUE))[["elapsed"]]

    if (!is.null(rq$sweep)) {
      sw <- rq$sweep
      sw$rep_dir <- d; sw$scenario <- scen; sw$tree_variant <- TREEV; sw$clone <- cl
      sweeps[[length(sweeps) + 1]] <- sw
    }
    rows[[length(rows) + 1]] <- data.frame(
      rep_dir = d, scenario = scen, tree_variant = TREEV, clone = cl, n_cells = nrow(sub),
      frac_state1 = mean(sub$state),
      PATH_k = rp$k, PATH_pi = rp$pi, PATH_pswitch = rp$p_switch, PATH_note = rp$note,
      PATH_secs = tp,
      PASTRI_k = rq$k, PASTRI_pi = rq$pi, PASTRI_pswitch = rq$p_switch,
      PASTRI_note = rq$note, PASTRI_secs = tq,
      SCOUT_alpha = ro$alpha, SCOUT_theta = ro$theta, SCOUT_sigma2 = ro$sigma2,
      SCOUT_model = ro$sel, SCOUT_note = ro$note, SCOUT_secs = ts,
      stringsAsFactors = FALSE)
  }
  if (ri %% 5 == 0 || ri == length(reps))
    cat(sprintf("  %d/%d replicates, %.1f min elapsed\n", ri, length(reps),
                as.numeric(difftime(Sys.time(), tstart, units = "mins"))))
}

out <- do.call(rbind, rows)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
write.csv(out, OUT, row.names = FALSE)
cat("wrote", OUT, "with", nrow(out), "clone-level rows\n")
cat("PATH finite k:", sum(is.finite(out$PATH_k)), "/", nrow(out), "\n")
cat("PASTRI finite k:", sum(is.finite(out$PASTRI_k)), "/", nrow(out), "\n")
cat("SCOUT finite alpha:", sum(is.finite(out$SCOUT_alpha)), "/", nrow(out), "\n")
if (length(sweeps)) {
  sw <- do.call(rbind, sweeps)
  swf <- sub("\\.csv$", "_pastri_depth_sweep.csv", OUT)
  write.csv(sw, swf, row.names = FALSE)
  cat("wrote", swf, "with", nrow(sw), "rows\n")
}
