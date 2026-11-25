options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 12, timeout = 1800)
lib <- .libPaths()[1]
cat("LIB:", lib, "\n")
need <- function(pkgs) {
  miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
  if (length(miss)) { cat("INSTALLING:", paste(miss, collapse=" "), "\n"); install.packages(miss) }
  for (p in pkgs) cat(sprintf("%-18s %s\n", p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING"))
}
cat("\n### PATH deps ###\n")
need(c("dplyr","magrittr","data.table","Matrix","expm","castor","qlcMatrix","phytools","ape"))
cat("\n### SCOUT deps ###\n")
need(c("OUwie","foreach","doParallel","reshape2","phylolm","corpcor","paleotree","nloptr",
       "future","future.apply","progressr","stringr","rlang","remotes"))
cat("\n### PASTRI deps ###\n")
need(c("tidyverse","argparse","patchwork","tidyr"))
cat("\nDONE_DEPS\n")
