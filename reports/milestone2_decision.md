# Study D, Milestone 2: decision

**Date:** 2026-11-23. **Pre-registration** frozen 18:13:52 +0530 before any EMT state or metastasis-versus-primary
contrast was computed (`checksums/FREEZE_SHA256_m2_simeonov.txt`). **Deviations D5–D9** logged in
`protocols/deviations.md`; D6 and D8 changed conclusions and are flagged there.

## Decision: **GO, with the scope narrowed to a methods-and-re-analysis paper**

Milestone 1 ended in REFORMULATE because its bound was classical and no data could test it. Both of those are now
resolved: the estimand uses lineage structure, it was validated against known ground truth, and it was applied to
the dataset the Milestone 1 audit identified. The result does **not** support the brief's top-tier biology venues,
and that is stated rather than finessed.

| Criterion from the brief | Verdict | Evidence |
|---|---|---|
| **Mathematical progress** beyond the Milestone 1 scalar bound | **Yes** | Decomposition Cov = A(L) + B with B = Var_c(π_c); identification condition T3; founder-bottleneck inflation T4; two-step estimator. Proofs in `theory/proofs/T1_T2.md`. |
| **Computational progress** against ground truth and alternatives | **Yes** | Individual-cell simulator with full ancestry; 5 mechanisms × 15 replicates; **0 false positives in 45 null replicates**; stress tests for misclassification, label loss and sample size. |
| **Biological progress** — reproduce a published result and add to it | **Yes** | Reproduced the mouse-1 metastatic EMT enrichment, then showed it is accounted for by clonal composition rather than within-clone state change. |
| **Integration** | **Yes** | The same estimand, estimator and identification condition run through theory, simulation and data; the identification condition failed in the real data exactly where the theory says it would. |
| **Generalisability** | **Partial** | GSE161363 (Quinn) is downloaded and audited but not yet analysed; n = 2 animals in one model. |

## What was established

**1. The metastatic EMT enrichment in this experiment is a change in clonal composition.** In mouse 1 the
composition term (+0.40 to +0.47 EMT units) exceeds the total primary-to-metastasis difference at all three sites,
and it is computed from primary-tumour expression only, so no technical difference between compartments can
produce it. The metastases are dominated by clone 1, which was already the most mesenchymal clone in the primary
(+0.58 against +0.05, −0.81, +0.07). Mouse 2 shows no enrichment at all.

**2. Organ induction adds almost nothing.** A model with environment-dependent switching alone is 0.6 AIC better
than a single-state model in mouse 1 and 16 AIC *worse* in mouse 2. Where an environment effect is supported
(mouse 2, M4), its sign is **toward the epithelial state**, not away from it.

**3. Clone identity predicts an unseen metastatic site.** Held out the entire peritoneal compartment: mean
absolute error on per-clone high-EMT fraction falls from 0.288 to 0.182 (mouse 1) and 0.209 to 0.089 (mouse 2).

**4. But the standard evidence for heritable *switching propensity* does not survive the decomposition.** The
nested model comparison favours clone-specific states by ΔAIC 250–440, and the naive clone-variance statistic is
large — yet simulation shows both behave identically under pure stable-state selection (naive statistic 0.040 when
the truth is 0). Applying the decomposition, heritable propensity is supported in **1 of 8** mouse-by-site
combinations (mouse 2 primary, B = 0.042, CI [0.005, 0.081]). In the other 7 the decay rate was unresolvable, so
by Proposition T3 the split is not identified there at all.

**5. Within-clone state shifts are not interpretable in these data.** Many held-out programmes and a pure QC
covariate shift as much as or more than EMT between compartments (genes detected moves 1.05 SD in mouse 2). The
within-clone EMT residual is therefore withdrawn as a biological claim.

## What was not established
- Whether metastatic competence is inherited as a *capacity to switch* rather than as a *state*: supported in one
  compartment of one animal, not in the other seven.
- Any human relevance. No human data were analysed in Milestone 2; the Milestone 1 audit already showed the NSCLC
  atlas is unpaired.
- Absolute switching and growth rates: a fixed-tree likelihood conditions on the tree, so these need a second
  time point (T5.4).
- Transfer to a second experimental system: GSE161363 is downloaded but not analysed.

## Honest novelty position
Lemma T1 is textbook and the plateau follows from the law of total covariance. PATH (Schiffman 2024) already
measures phylogenetic autocorrelation of cell states, and **PASTRI (Yang 2026) already uses pairs at varying
phylogenetic distance** — found on the final novelty check and disclosed in
`mathematical_novelty_assessment.md`. What appears to be unclaimed is the decomposition of that curve into state
memory and a propensity asymptote, the proof that induction and stable-state selection cannot contribute to the
asymptote, the identification condition, the bottleneck inflation term, and the two cautionary counterexamples
(C6, C8) showing that the two analyses a reader would run by default both select "heritable propensity" when it
is absent. Perelli et al. (*Nature* 2025) reach a convergent biological conclusion by genetic ablation.

## Next experiments, in priority order
1. **Transfer to GSE161363 (Quinn).** Data are already on disk. This is the cheapest test of whether the
   composition-versus-within-clone result is specific to the macsGESTALT system. ~1 day.
2. **Propagate lineage-reconstruction uncertainty.** The authors' trees and clone calls are used as given. Bootstrap
   over indel characters to see whether the mouse-2 primary plateau survives.
3. **Design the measurement the theory asks for.** By T3 the binding constraint is the *spread of MRCA depths*, not
   cell number: a recorder tuned for a wider range of edit depths in the primary tumour would make B identifiable.
   This is a concrete, falsifiable experimental-design recommendation and the most valuable follow-up.
4. **Only then** consider human translation, and only for programme conservation, never for mechanism.

## Separation of claims
- **Mathematical:** T1–T5, proved. T1 classical.
- **Simulation:** E8, E9, with ground truth.
- **Real-data observation:** the composition decomposition, the model comparison, the held-out prediction, the
  single positive plateau.
- **Biological interpretation:** clone-level selection of an already-mesenchymal lineage, in one mouse model, n = 2.
- **Not claimed:** within-clone induction or its absence; human mechanism; absolute rates.
