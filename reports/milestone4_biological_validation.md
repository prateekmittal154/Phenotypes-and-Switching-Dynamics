# Milestone 4: biological validation across two lineage systems

Six independent animals in total: two from macsGESTALT/pancreatic (Simeonov) and four from
Cassiopeia/A549-lung-xenograft (Quinn). States are defined within each system; the Quinn states never use the
Simeonov EMT axis.

## 1. Simeonov reproduced exactly
The Milestone 2 and 3 results are unchanged and were re-verified against the frozen checksums. In mouse 1 the
clone-composition term (+0.40 to +0.47 EMT units) exceeds the total primary-to-metastasis difference at all
three metastatic sites, and is computed from primary-tumour expression reweighted by metastatic clone
abundances, so no cross-compartment technical effect can produce it. Mouse 2 shows no enrichment. Clone
occupancy is heterogeneous ($V_\pi$ 0.031–0.053, $p<10^{-10}$); switching rates are not assessable.

## 2. Quinn: clone-associated occupancy replicates in an independent system

Cell-cycle state (fixed external gene set), all sites, per-clone maximum likelihood:

| Mouse | clones | $V_\pi$ | shared-occupancy test | clones with rate off-bound | $V_{\log k}$ | shared-rate test |
|---|---|---|---|---|---|---|
| 5k | 30 | 0.0072 | $p<10^{-16}$ | 28/30 | 0.654 | $p=0.0007$ |
| 30k | 3 | 0.0043 | $p=0.030$ | 3/3 | 0.000 | $p=0.22$ |
| 10k | 11 | 0.0103 | $p<10^{-8}$ | 11/11 | 0.514 | $p=0.56$ |
| 100k | 15 | 0.0429 | $p<10^{-16}$ | 13/15 | 0.935 | $p=0.0044$ |

**Occupancy heterogeneity is present in all four animals.** It also replicates under the second, independent
state definition (unsupervised PC1): $V_\pi$ = 0.0053, 0.0182, 0.0047, 0.0114 with $p<10^{-4}$ in every mouse.
Two different state axes, four animals, one conclusion.

**More clone rates are estimable here than in Simeonov**: 55 of 59 clone fits returned rates away from the
optimiser bounds, against 18 of 35. It is tempting to credit the recorder — Cassiopeia resolves 8–14 coalescent
levels against macsGESTALT's 3–4 — and an earlier version of this report did so. **The empirical check does not
support that.** Computing the covered band from the realised distribution of pair distances rather than from the
idealised capacity law (`milestone4_design_optimization.md` §4), Cassiopeia's ceilings are 0.09–0.30 per
resolved branch, *below* macsGESTALT's 0.83–2.00, because its trees are star-like and short-distance pairs are
scarce. The fitted Quinn rates (median 0.72–1.09 per resolved branch) therefore lie **above** their own covered
band, by a factor of three to ten. Being away from an optimiser bound is not the same as being identifiable.

## 3. But the within-clone shifts are technically confounded, again
The same control applied to Simeonov in Milestone 2 was applied here
(`results/tables/m4_quinn_technical_control.csv`), measuring the within-clone primary-to-metastasis shift in
primary-SD units for the state and for quality covariates:

| Mouse / site | cell-cycle state | genes detected | library size |
|---|---|---|---|
| 5k M1 | −0.64 | **−1.39** | **−0.95** |
| 5k M2 | −0.48 | **−1.10** | **−0.85** |
| 5k Liv | +0.34 | **+0.59** | **+0.66** |
| 30k M2 | −0.34 | −0.22 | −0.21 |

In the largest mouse the quality covariates move **more than twice as far as the biological signal, and in the
same direction**. Cell-cycle scores are known to scale with library complexity. The within-clone cell-cycle
shifts are therefore **not interpretable as biology** and no claim is made from them. This mirrors the
Milestone 2 finding for EMT and suggests the problem is generic to cross-tissue single-cell comparisons rather
than specific to either dataset.

## 3b. Animal-level eligibility, stated explicitly

Not every animal supports every analysis, and the reasons are structural rather than post-hoc.

| Mouse | lineage groups | usable trees, all sites | usable trees, primary only | occupancy analysis | switching analysis | depth-adjusted analysis |
|---|---|---|---|---|---|---|
| 5k | 100 | 30 (cap) | 30 | yes | yes | yes |
| 10k | 13 | 11 | 7 | yes | yes | yes |
| 100k | 31 | 15 | 10 | yes | yes | yes |
| **30k** | **4** | **3** | **1** | all-sites only | all-sites only | **ineligible** |

**The 30k animal is excluded from the depth-adjusted analysis.** Its barcode library has only four lineage
groups, of which three give trees above the 25-cell threshold across all sites and just one does so within the
primary site. A variance across clones cannot be estimated from one clone, so the primary-only analysis covers
three animals, not four. Its all-sites estimate ($V_\pi = 0.0043$, $p = 0.030$ on three clones) is reported for
completeness and carries little weight.

Two eligibility rules were applied uniformly and before seeing results (deviation D20): only the authors'
primary `hybrid_priors.alleleThresh` reconstruction is used, and trees above 1,500 matched cells are skipped
with at most 30 fitted per mouse, to bound the likelihood-ratio optimisation. The 5k animal is the only one
where the 30-tree cap binds.

## 4. Does the occupancy result survive the technical confound?
The decisive test restricts to **primary-site cells only** (one tissue, so no cross-tissue contrast) and
regresses the state on log library size and genes detected before binarising:

| Mouse | clones | $V_\pi$ | shared-occupancy test | $V_{\log k}$ | shared-rate test |
|---|---|---|---|---|---|
| 5k | 30 | 0.0021 | $p=0.024$ | 0.161 | $p=0.65$ |
| 10k | 7 | 0.0288 | $p<10^{-9}$ | 0.000 | $p=0.98$ |
| 100k | 10 | 0.0321 | $p<10^{-16}$ | **1.664** | $p=0.0042$ |
| 30k | 1 usable tree | — | — | — | — |

**Clone-associated occupancy survives in all three analysable mice**, though attenuated roughly three-fold in
the 5k mouse (0.0072 → 0.0021, $p$ 0.024). This is the strongest evidence in the project that occupancy
heterogeneity is real rather than technical: it holds in a second tumour model, a second recorder chemistry,
two state definitions, and within a single tissue with sequencing depth adjusted.

**Switching-rate heterogeneity does not survive as a general claim.** After depth adjustment it is significant
in **1 of 3** mice (100k, $p = 0.0042$), down from 2 of 4. A single animal is not replication. Two further
cautions compound this: the fitted rates lie three- to ten-fold above the empirically covered band (above), and
the calibration work shows the shared-occupancy test is mildly anti-conservative under misspecification. The
defensible statement is therefore weaker than the one this project hoped for: **no dataset analysed here can
identify switching-rate heterogeneity, and the one nominally significant animal sits outside the regime where
the estimate is trustworthy.** Testing the hypothesis needs a recorder built for it.

## 5. Compositional mechanism differs by phenotype
A useful contrast falls out of the same decomposition. In Simeonov, the metastatic EMT difference is
**compositional** (which clones seeded). In Quinn, the primary-to-metastasis cell-cycle difference is
**within-clone** (composition contributes −0.01 to +0.03 of a total of −0.22 at 5k M1), i.e. the site changes
proliferation across clones rather than selecting slow- or fast-cycling clones. That is biologically sensible —
but given §3 the within-clone part is confounded here, so this contrast is reported as a hypothesis the
decomposition can test, not as an established result.

## 6. Gene-level validation
Not repeated for Quinn. In Simeonov, nine held-out programmes gave no association with clone occupancy that
replicated across the two animals, with three changing sign. Quinn's state is itself a gene programme
(cell cycle), so an independent-gene test would need a third programme set and, given §3, would face the same
depth confound. Deferred rather than run badly.

## 7. Spatial and histological validation
Not performed. No public dataset with matched lineage, expression and spatial coordinates was located
(see `milestone4_dataset_expansion.md`). No imaging component was forced.

## 8. What the depth controls mean for the science

The sequencing-depth confound is not a footnote; it shapes the main interpretation in three ways.

1. **It rules out the within-clone channel in both systems.** In Simeonov several unrelated programmes and a
   quality covariate moved as much as the EMT score between compartments; in Quinn, genes detected moved 1.39 SD
   where the cell-cycle state moved 0.64 SD, in the same direction. Any mechanism that acts *within clones
   across tissues* — environmental induction, organ-specific reprogramming — is therefore not measurable in
   either dataset. This is why the compositional channel carries the whole interpretation.
2. **It is why the composition term is the trustworthy statistic.** The composition term is built from
   primary-tumour expression reweighted by metastatic clone abundances, so it never compares expression across
   tissues and is structurally immune to the confound. That is not a convenient accident; it is the reason to
   prefer it.
3. **It sets the standard for the occupancy claim.** Because occupancy is estimated from clones that span
   tissues, it could in principle inherit the same confound. The primary-site-only, depth-adjusted analysis
   exists precisely to remove that route, and the claim survives it in all three analysable mice. Readers should
   weight that analysis, not the all-sites one, as the evidence for occupancy heterogeneity.

A generic lesson follows for the field: **single-cell comparisons across dissociated tissues should report
quality covariates on the same scale as the biological signal.** Had we not, we would have reported an
environmental cell-cycle effect in Quinn and a within-clone EMT reversal in Simeonov, both of which are better
explained by library complexity.

## 9. Summary of what is now established biologically
1. **Clone-associated state occupancy is real and replicates**: two tumour models, two recorders, six animals,
   two independent state definitions, and surviving a within-tissue depth-adjusted analysis.
2. **Metastatic EMT enrichment in the pancreatic model is compositional**, not a within-clone state change.
3. **Switching-rate heterogeneity remains unproven, and is not currently measurable**: significant in one
   animal only, and there the fitted rates lie outside the empirically covered band.
4. **Within-clone cross-tissue comparisons are confounded by sequencing depth in both datasets**, which is a
   methodological warning for the field as much as a limitation of this analysis.
