# Table B2: exact software provenance for every component used in the benchmark.
args <- commandArgs(trailingOnly = TRUE); B <- args[1]
gitinfo <- function(p) {
  if (!dir.exists(file.path(p, ".git"))) return(c(NA, NA))
  c(system(sprintf("git -C %s rev-parse HEAD", shQuote(p)), intern = TRUE)[1],
    system(sprintf("git -C %s log -1 --format=%%ai", shQuote(p)), intern = TRUE)[1])
}
pv <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
rows <- list()
add <- function(...) rows[[length(rows) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

for (m in c("PATH", "SCOUT", "PASTRI_package", "PASTRI_paper_scripts")) {
  gi <- gitinfo(file.path(B, "methods", m))
  add(component = m, kind = "method repository", version = "",
      commit = gi[1], commit_date = gi[2],
      source = switch(m,
        PATH = "https://github.com/landau-lab/PATH",
        SCOUT = "https://github.com/hrstuart/SCOUT",
        PASTRI_package = "https://github.com/yuxiaochen11/PASTRI",
        PASTRI_paper_scripts = "https://github.com/Yangwj01/PASTRI"))
}
add(component = "R", kind = "language", version = paste(R.version$major, R.version$minor, sep = "."),
    commit = "", commit_date = "", source = "conda-forge r-base")
for (p in c("ape", "OUwie", "phylolm", "paleotree", "corpcor", "nloptr", "expm", "castor",
            "qlcMatrix", "phytools", "Matrix", "data.table", "dplyr", "tidytree", "treeio",
            "fgsea", "ggplot2", "devtools"))
  add(component = p, kind = "R package", version = pv(p), commit = "", commit_date = "",
      source = "CRAN/Bioconductor")
out <- do.call(rbind, rows)
dir.create(file.path(B, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path(B, "results", "tables", "TableB2_software.csv"), row.names = FALSE)
print(out[, c("component", "kind", "version", "commit")], row.names = FALSE)
cat("\nggtree:", pv("ggtree"),
    "-- unavailable; PASTRI's fortify() is supplied from tidytree, see src/adapters/pastri_load.R\n")
