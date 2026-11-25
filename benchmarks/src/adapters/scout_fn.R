# The OUwie (SCOUT engine) per-clone fit, factored out so run_scout_real.R can reuse
# the identical function rather than a copy that could drift.
suppressPackageStartupMessages({ library(ape); library(OUwie) })

run_ouwie_clone <- function(phy, trait, tips) {
  if (is.null(phy$edge.length)) phy$edge.length <- rep(1, Nedge(phy))
  phy$edge.length[phy$edge.length <= 0] <- 1e-7
  # The published trees can carry unary (single-child) nodes, including at the root.
  # multi2di does not remove those, so OUwie still refuses the tree as "not rooted and
  # fully dichotomous". collapse.singles drops them first.
  phy <- collapse.singles(phy)
  phy <- multi2di(phy)
  phy$edge.length[phy$edge.length <= 0] <- 1e-7
  if (!is.rooted(phy) || !is.binary(phy))
    return(list(alpha = NA_real_, theta = NA_real_, sigma2 = NA_real_,
                sel = NA_character_, note = "tree not rooted/dichotomous after resolution"))
  # resolving the tree can reorder tips, so re-align the trait to the new order
  ord <- match(phy$tip.label, tips)
  if (anyNA(ord))
    return(list(alpha = NA_real_, theta = NA_real_, sigma2 = NA_real_,
                sel = NA_character_, note = "tip labels lost during resolution"))
  trait <- trait[ord]
  tips <- phy$tip.label
  phy$node.label <- rep("1", Nnode(phy))
  dd <- data.frame(species = tips, reg = "1", x = as.numeric(trait), stringsAsFactors = FALSE)
  fo <- tryCatch(OUwie::OUwie(phy, dd, model = "OU1", simmap.tree = FALSE,
                              scaleHeight = FALSE, quiet = TRUE, diagn = FALSE),
                 error = function(e) NULL)
  fb <- tryCatch(OUwie::OUwie(phy, dd, model = "BM1", simmap.tree = FALSE,
                              scaleHeight = FALSE, quiet = TRUE, diagn = FALSE),
                 error = function(e) NULL)
  if (is.null(fo)) return(list(alpha = NA_real_, theta = NA_real_, sigma2 = NA_real_,
                               sel = NA_character_, note = "OU1 error"))
  al <- unname(fo$solution["alpha", 1]); s2 <- unname(fo$solution["sigma.sq", 1])
  th <- tryCatch(unname(fo$theta[1, 1]), error = function(e) NA_real_)
  sel <- if (!is.null(fb)) { if (fo$AICc < fb$AICc) "OU" else "BM" } else NA_character_
  list(alpha = al, theta = th, sigma2 = s2, sel = sel, note = "")
}
