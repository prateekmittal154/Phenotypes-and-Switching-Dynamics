# SCOUT (OUwie engine) on the real clone trees.
#
# Run separately because the first real-data pass failed for SCOUT on every clone: the
# published trees carry unary nodes that multi2di does not remove, so OUwie rejected them
# as "not rooted and fully dichotomous". run_competitors.R now calls collapse.singles()
# first; this script re-runs only the SCOUT column so the expensive PATH and PASTRI work
# is not repeated.
#
# Usage: Rscript run_scout_real.R <benchdir> <out_csv>
suppressPackageStartupMessages({ library(ape); library(OUwie) })
args <- commandArgs(trailingOnly = TRUE); B <- args[1]; OUT <- args[2]
source(file.path(B, "src", "adapters", "scout_fn.R"))

dirs <- list.dirs(file.path(B, "data", "real_trees"), recursive = FALSE)
cat("clone trees:", length(dirs), "\n")
rows <- list(); t0 <- Sys.time()
for (i in seq_along(dirs)) {
  d <- dirs[i]
  nwk <- file.path(d, "tree_unit.nwk")
  if (!file.exists(nwk)) next
  tr <- read.tree(nwk)
  st <- read.csv(file.path(d, "states.csv"), stringsAsFactors = FALSE)
  sub <- st[match(tr$tip.label, st$cell), , drop = FALSE]
  tsec <- system.time(
    ro <- tryCatch(run_ouwie_clone(tr, sub$trait, sub$cell),
                   error = function(e) list(alpha = NA_real_, theta = NA_real_,
                                            sigma2 = NA_real_, sel = NA_character_,
                                            note = paste("error:", conditionMessage(e))))
  )[["elapsed"]]
  base <- basename(d)
  rows[[length(rows) + 1]] <- data.frame(
    rep_dir = d, tree_variant = "unit",
    clone = sub("^.*__", "", base), animal = sub("__.*$", "", base),
    n_cells = nrow(sub),
    SCOUT_alpha = ro$alpha, SCOUT_theta = ro$theta, SCOUT_sigma2 = ro$sigma2,
    SCOUT_model = ro$sel, SCOUT_note = ro$note, SCOUT_secs = tsec,
    stringsAsFactors = FALSE)
  if (i %% 20 == 0)
    cat(sprintf("  %d/%d, %.1f min\n", i, length(dirs),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
out <- do.call(rbind, rows)
write.csv(out, OUT, row.names = FALSE)
cat("wrote", OUT, "with", nrow(out), "rows; finite alpha:",
    sum(is.finite(out$SCOUT_alpha)), "/", nrow(out), "\n")
print(table(out$SCOUT_note))
