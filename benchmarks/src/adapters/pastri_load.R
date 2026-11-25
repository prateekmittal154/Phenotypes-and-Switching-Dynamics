# Load PASTRI's own estimator source.
#
# DEVIATION, documented: the PASTRI package DESCRIPTION imports ggtree, but ggtree could
# not be installed in this environment -- the Bioconductor build for R 4.4 calls
# ggplot2:::check_linewidth (removed in ggplot2 4.x) and ggplot2::is_ggplot (added in
# ggplot2 4.x), so it fails against both 3.5.1 and 4.0.4; the GitHub build additionally
# needs ggiraph -> gdtools, whose system dependencies are absent.
#
# PASTRI uses ggtree for exactly one thing: fortify(tree), to turn a phylo object into a
# data frame, and it then reads only the `parent` and `node` columns
# (R/LCA_depths_calculation.R). We supply that single conversion from tidytree, which is
# by the same authors and is what ggtree's fortify.phylo is built on. The estimator itself
# -- adjacency construction and the Frobenius-norm optimisation -- is the authors' own
# unmodified code, sourced verbatim below.
suppressPackageStartupMessages({
  library(ape); library(dplyr); library(tidyr); library(parallel); library(tidytree)
})

fortify <- function(model, data, ...) {
  df <- as.data.frame(tidytree::as_tibble(model))
  # ggtree's fortify.phylo marks the root with parent == node; tidytree already does this.
  df
}

.pastri_dir <- Sys.getenv("PASTRI_DIR")
if (!nzchar(.pastri_dir)) stop("set PASTRI_DIR")
for (f in c("LCA_depths_calculation.R", "Transition_matrix_calculation.R",
            "Transition_matrix_optimization.R")) {
  source(file.path(.pastri_dir, "R", f), local = FALSE)
}
