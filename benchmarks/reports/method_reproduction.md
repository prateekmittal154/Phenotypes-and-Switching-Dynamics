# Method reproduction

Before any comparative claim, each method was installed from its authors' own repository
and run on the authors' own example. Raw console output is in
`reports/method_reproduction_raw.txt`; the machine-readable summary is
`results/metrics/method_reproduction.csv`.

Environment: R 4.4.3 (conda-forge) in `environments/renv`, on the analysis server. Exact
package versions are in `results/tables/TableB2_software.csv`.

---

## PATH — reproduced

**Source** `https://github.com/landau-lab/PATH`, commit `caffc5105678224f108b3dd61fea545f96235ffc`,
tag v1.0, dated 2024-08-20. Loaded with `devtools::load_all()`.

**Check 1 — the documented `xcor()` example.** The example in the package's own roxygen
block (`R/xcor.R`) runs unchanged and returns a phylogenetic autocorrelation of **0.8317**
with z = **15.16** on a tree simulated from the package's own generator.

**Check 2 — recovering a rate matrix the package itself simulated.** Using
`ran_tree_with_states()` with a known symmetric Q at rate 0.05, 512 leaves, 20 replicates:
`PATH_inf()` returned a finite Q in **20 of 20** replicates, with a median off-diagonal rate
of **0.0469** against the true **0.05**. The installation therefore recovers a known
generating rate to within 6%.

**Install note.** PATH imports `fgsea` (used only by `xcor_gsea`). `fgsea` requires
ggplot2 >= 3.5.2; the environment was pinned to ggplot2 3.5.2 to satisfy it.

---

## PASTRI — reproduced, with one documented substitution

**Source.** Two repositories exist. The preprint's code-availability statement names
`https://github.com/Yangwj01/PASTRI` (commit `5ecc99dfc5b7ede11bd20bf45ae47a1883267d5e`),
which holds the paper's analysis scripts. The runnable R package is
`https://github.com/yuxiaochen11/PASTRI` (commit `67aefd96240f9d3a1900c6906d813962233ef38f`),
which the first repository's README directs users to install. Both were cloned; the package
is what we run. This is the genuine 2026 PASTRI of Yang, Li, Yu et al., not the older
identically named tumour-phylogeny method.

**Check — the shipped HSC example.** The package ships `inst/extdata/HSC.nwk` (778 tips) and
`HSC_nodeInfos.csv` (916 rows, 8 cell types). Following the README:

- `calculate_lca_depths()` returned **302,253** leaf pairs in ~36 s.
- `get_optimal_transition_matrix()` at `fi_depth = 0.5` with the README's
  upper-triangular `Bound_Matrix` returned an 8x8 matrix whose **rows sum to 1** and whose
  lower triangle is zero, i.e. the irreversibility constraint is respected.

**Substitution, and why.** The package imports `ggtree` solely to call `fortify(tree)`,
reading only the `parent` and `node` columns. `ggtree` could not be installed here: the
Bioconductor build for R 4.4 calls `ggplot2:::check_linewidth` (removed in ggplot2 4.x) and
`ggplot2::is_ggplot` (added in 4.x), so it fails against both 3.5.1 and 4.0.4, and the
GitHub build additionally requires `ggiraph` -> `gdtools`, whose system dependencies are
absent. We therefore supply `fortify` from `tidytree` (same authors; ggtree's
`fortify.phylo` is built on it) and source PASTRI's three estimator files verbatim. See
`src/adapters/pastri_load.R`. **No part of the estimator was reimplemented** — the
adjacency construction and the Frobenius-norm optimisation are the authors' code.

---

## SCOUT — engine reproduced; full pipeline not exercised

**Source** `https://github.com/hrstuart/SCOUT`, commit
`c72befe0d052f3960c448a64bff776665ab8ea4c`, dated 2026-10-21. Engine: `OUwie` 3.0.3.

**Check — the shipped simulation example.** `examples/sim_example/` contains a 256-tip
Cassiopeia-NJ tree and a counts matrix whose column names encode the generating model
(30 `BM1_*`, 30 `OU1_*`, 31 `OUM_*`). Tips map to rows as `t<rowname>` (256/256 matched).
**The shipped tree carries no branch lengths**; we substituted unit edges, which is SCOUT's
own documented rule (`R/dropout_validation.R:470-472`: *"No edge lengths, replacing with
1s."*).

Fitting OUwie `OU1` and `BM1` to 20 genes from each generating class and selecting by AICc:

| Generating model | selected BM1 | selected OU1 |
|---|---|---|
| BM1 | 5 | 15 |
| OU1 | 0 | 20 |

Overall accuracy **0.625** against a 0.5 chance baseline, with all 40 genes fitting
successfully. The engine identifies every OU-generated gene but calls 15 of 20
BM-generated genes OU as well.

**A scale check worth recording.** The shipped file is named
`3state_n256_alpha_3_sigma_1_counts.csv`, i.e. generated at alpha = 3. Fitting OU1 to the
20 `OU1_*` genes on the unit-edge tree returns a median alpha of **0.309** — about
ten-fold below the generating value. This is not a defect of the installation: the tree
carries no branch lengths, so the fitted alpha is expressed in units of the imposed unit
edges rather than the generating timescale. It is a concrete instance of the rate-time
confounding that Benchmark III examines, visible in a method's own shipped example.

**What this does and does not show.** It shows the installed OU machinery runs on the
authors' data and discriminates above chance. It is **not** a test of SCOUT's full
pipeline: `SCOUT()` adds Pagel's-lambda screening, an EM step, regime assignment and a
hybrid selection rule, none of which were used here. The directional bias toward OU should
therefore be read as a property of a bare OU1-vs-BM1 AICc comparison on a reconstructed
tree with imposed unit edges, not as a defect of SCOUT as published. SCOUT's own
`runSCOUT.identifiability()` exists precisely to probe this kind of question and was not
run.

**Install note.** SCOUT's `NAMESPACE` imports `TedSim` (GitHub, simulation only). We use
the OUwie engine directly on our own simulated data, so `TedSim` was not required.

---

## Summary

| Method | Installed from authors' source | Authors' own example reproduced | Caveat |
|---|---|---|---|
| PATH | yes, v1.0 `caffc51` | yes — recovers a known Q (0.0469 vs 0.05) | none |
| PASTRI | yes, `67aefd9` | yes — shipped HSC example, rows sum to 1, constraint respected | `fortify` supplied from tidytree; estimator unmodified |
| SCOUT | yes, `c72befe` (engine OUwie 3.0.3) | partially — OU/BM selection at 0.625 on the shipped example; median alpha 0.309 against a generating 3 | full SCOUT pipeline (EM, lambda screen, hybrid) not exercised; unit edges imposed because the shipped tree has none |
| Study D | n/a, this work | Milestone 6 verified version | none |
