options(repos=c(CRAN="https://cloud.r-project.org"), Ncpus=12, timeout=2400)
for (p in c("png","Rmpfr","geiger","paleotree","OUwie")) {
  if (!requireNamespace(p, quietly=TRUE)) { cat("## installing", p, "\n"); try(install.packages(p)) }
  cat(sprintf("%-12s %s\n", p, if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING"))
}
cat("DONE3\n")
