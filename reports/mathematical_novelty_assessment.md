# Mathematical novelty assessment (Milestone 2)

An attempt to find prior work that already contains the Milestone 2 result, and an honest grading of what is new.

## The claim under assessment
For within-clone cell pairs at lineage distance L, with states centred at the compartment mean,
Cov = E_c[π_c(1 − π_c)e^{−2k_cL}] + Var_c(π_c); the asymptote identifies heritable switching propensity; it is
immune to organ induction and to stable-state selection; it is identified only when the decay is resolvable; and
founder bottlenecks inflate it by π(1 − π)/N_eff.

## Prior art, component by component

| Component | Status | Closest prior work |
|---|---|---|
| Cov = π(1 − π)e^{−2kL} for a two-state Markov process on a tree | **Classical** | Felsenstein 1985 (phylogenetic covariance); the binary analogue of the Ornstein–Uhlenbeck trait covariance used for phylogenetic signal (Pagel's λ, Blomberg's K) |
| Plateau = between-group variance | **Classical** | the law of total covariance; variance-components models |
| Phylogenetic autocorrelation of *cell states* in cancer | **Published** | PATH, Schiffman et al., *Nat Genet* 2024 — Moran's-I style autocorrelation of states on single-cell phylogenies, reported as heritability versus plasticity |
| Ornstein–Uhlenbeck modelling of expression on single-cell lineage trees | **Published** | SCOUT, Stuart et al. 2025 — fits OU processes per gene on lineage trees to separate drift from selection |
| Clone-level versus cell-level attribution of resistance | **Published** | Whiting et al., *Nat Commun* 2025; Goyal et al. (FateMap) 2023; Kim et al. 2018 |
| Using cell pairs at *varying* phylogenetic distance rather than one aggregate | **Published (2026)** | **PASTRI**, Yang et al. 2026 (bioRxiv) — "by leveraging cell pairs at varying phylogenetic distances, PASTRI accurately resolves stage-specific transition rates". Found only on the final novelty check; it removes the depth-resolved curve itself from the novelty claim. PASTRI's target is how transition rates change with developmental stage; it does not decompose the curve into a decaying term and an asymptote, and does not interpret an asymptote as between-clone propensity variance. |
| **Identifying the asymptote of that curve with Var_c(π_c), the heritable switching propensity** | **Not found** | PATH aggregates over the tree; SCOUT is per-gene and continuous; PASTRI targets stage-specific rates |
| **Proof that organ induction and stable-state selection contribute nothing to that asymptote** | **Not found** | — |
| **The identification condition (A and B separable iff the decay is resolvable) and the exact k = 0 equivalence** | **Not found** | analogous in spirit to Guan et al. 2025 (SDE drift and diffusion identifiable iff observed out of equilibrium) |
| **Founder-bottleneck inflation π(1 − π)/N_eff of the asymptote** | **Not found** | related to effective population size and the inverse Simpson index in population genetics, but not applied to this estimand |

## Grading
- **New theorem:** no. Lemma T1 and the total-covariance step are textbook.
- **New observation model and estimand:** partly. The depth-resolved curve itself is **not** new — PASTRI (Yang
  et al. 2026) already uses pairs at varying phylogenetic distance, which was found only on the final novelty
  check and is disclosed here rather than quietly dropped. What remains unfound is the *decomposition* of that
  curve into a decaying state-memory term and an asymptote, and the identification of the asymptote with
  Var_c(π_c), the heritable switching propensity.
- **New identification result:** yes, in the modest sense that T3 states exactly when the two are separable and
  exhibits the equivalence when they are not. This is the kind of statement that is obvious once posed and that
  nobody appears to have posed for this estimand.
- **New estimator:** yes. The two-step procedure is forced by T3 and is what makes the decomposition usable.
- **New negative methodological result:** yes, and this may be the most useful part. Counterexamples C6 and C8
  show that the two analyses a reader would naturally run — the variance of clone-level state fractions, and a
  nested model comparison favouring clone-specific states — both select "heritable propensity" under pure
  stable-state selection. Simulation puts numbers on it (naive statistic 0.040 when the truth is 0).

## Convergent biological prior art
Perelli et al., *Nature* 2025 ("Evolutionary fingerprints of epithelial-to-mesenchymal transition") argue for
"cell-state-restricted patterns of evolution, whereby cancer cell speciation is propagated to progeny within
restricted functional compartments", also in pancreatic cancer. That is independent biological support for the
re-analysis finding here — EMT state tracking with lineage rather than shifting within it — arrived at by
genetic ablation rather than by re-analysis. It strengthens the biological claim and further weakens any
suggestion that the idea is itself novel.

## Comparison with Milestone 1
Milestone 1's bound was judged a direct composition of maximal coupling and a density-band projection, and the
decision was REFORMULATE. Milestone 2 is strictly stronger on every axis that mattered there: it uses lineage
information rather than two marginals; it yields a decomposition with confidence intervals rather than a bound
under an unverifiable parameter Λ; it was validated against known ground truth; and it was applied to real data
that the Milestone 1 audit identified as the test bed.

## Honest venue assessment
This supports a methods paper — *PLoS Computational Biology*, *Cell Systems*, *Genome Biology* or *Bioinformatics*
— whose contribution is an estimator, an identification condition and a cautionary re-analysis. It does not
support a mathematics venue, and the brief's top-tier biology venues would require more animals and a second
experimental system. The strongest single sentence that is fully defensible is: *the metastatic EMT enrichment in
the macsGESTALT pancreatic model is attributable to clonal composition rather than to within-clone state change,
and the standard evidence for heritable switching propensity does not survive a decomposition that separates it
from frozen state heritability.*
