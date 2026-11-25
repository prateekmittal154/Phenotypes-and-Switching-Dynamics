# Re-analysis of Simeonov et al. 2021 (macsGESTALT, GSE173958)

Pre-registered in `protocols/m2_simeonov_prereg.md`, frozen at 2026-11-23T18:13:52+0530 (SHA-256 in
`checksums/FREEZE_SHA256_m2_simeonov.txt`) **before any EMT state or metastasis-versus-primary contrast was
computed**. Only tissue counts, clone sizes and the authors' coarse clusters had been inspected.

## Data
- Processed lineage and cell tables from the authors' Mendeley deposit (doi:10.17632/t98pjcd7t6.1); both archives
  verified against their published SHA-256.
- Expression matrices from GEO GSE173958 (1.67 GB). All **28,028** cancer cells in the authors' tables matched the
  expression matrices, and every gene in the frozen panel was detected.
- Cells with a position in the authors' lineage trees: **10,608** (mouse 1, 29 clones) and **13,834** (mouse 2,
  66 clones).
- Independent units: 2 mice. Clones are the bootstrap unit within mouse; cells are never treated as units.

**State definition validity.** The EMT score, built from a frozen gene list and z-scored against primary-tumour
cells only, correlates with the authors' own independently derived EMT pseudotime at Spearman **ρ = 0.755**. The
score was not tuned to achieve this.

## H1 — reproduction, and where the enrichment comes from

Mouse 1 reproduces the published enrichment: metastatic cells are more mesenchymal than primary cells
(peritoneal **+0.40**, liver **+0.16**, lung **+0.30** EMT units; Mann–Whitney p < 1e-4). Mouse 2 does not
(−0.05, −0.002, −0.34).

Decomposing the mouse-1 difference into a change in **which clones are present** versus a **shift within clones**:

| Mouse / site | total | clone composition | within clone |
|---|---|---|---|
| M1 peritoneal | +0.399 | **+0.471** | −0.072 |
| M1 liver | +0.163 | **+0.397** | −0.234 |
| M1 lung | +0.302 | **+0.467** | −0.165 |
| M2 peritoneal | −0.049 | +0.037 | −0.087 |
| M2 liver | −0.002 | +0.159 | −0.161 |
| M2 lung | −0.344 | −0.132 | −0.212 |

**The composition term alone exceeds the total difference at every mouse-1 site.** The composition term is
computed from *primary-tumour* expression values reweighted by metastatic clone abundances, so it cannot be
produced by any technical difference between the metastatic and primary samples.

The per-clone table makes the mechanism concrete. In mouse 1 the metastases are dominated by clone 1 (3,229
peritoneal, 2,160 liver, 120 lung cells), and clone 1 was **already the most mesenchymal clone in the primary
tumour** (EMT +0.58, against +0.05, −0.81 and +0.07 for the other seeding clones). In mouse 2 a single clone
(clone 2) dominates every compartment, and its EMT score barely moves (+0.12 primary → −0.01 peritoneal).

**Reading.** The metastatic EMT enrichment in this experiment is a change in clonal composition — which lineages
founded metastases — rather than a shift of cells into mesenchymal states after arrival. That is consistent with
the authors' "selection" framing and sharpens it: the selection acts at the level of clones that were already
mesenchymal, not on cell states within clones.

**The within-clone residual is not interpretable as biology.** It is negative at all six mouse-site combinations,
but a sign test over clones is not significant anywhere (p = 0.125–1.00; only 1–6 paired clones per site), and the
technical control below shows larger shifts in unrelated programmes.

## Technical control (deviation D5)
The within-clone primary-to-metastasis shift, in primary-SD units, for every held-out programme and three QC
covariates:

| Programme | M1 peritoneal | M1 liver | M2 peritoneal | M2 liver |
|---|---|---|---|---|
| **EMT score** | **−0.13** | **−0.31** | **−0.09** | **−0.15** |
| hypoxia | −0.55 | −0.35 | −0.38 | −0.56 |
| glycolysis | −0.36 | −0.20 | −0.56 | −0.67 |
| chromatin | −0.24 | −0.41 | −0.37 | −0.18 |
| S100 | +0.38 | +0.53 | +0.38 | +0.43 |
| proliferation | +0.34 | +0.42 | +0.10 | +0.24 |
| ribosomal fraction | +0.71 | +0.47 | +0.91 | +1.02 |
| genes detected | +0.01 | +0.23 | −1.05 | −0.75 |

Many unrelated programmes, and a pure QC covariate (genes detected, up to 1.05 SD in mouse 2), shift as much as or
more than EMT. The compartments differ technically — dissociation, ambient RNA and depth all differ between
pancreas, liver and peritoneum. **The within-clone EMT shift cannot be separated from these effects and is not
claimed as biology.** The composition term is unaffected, because it uses only primary-tumour expression.

## H2 — state heritability (clone-stratified permutation, the pre-registered test)
| Mouse / site | observed Spearman(MRCA depth, |Δ EMT|) | null mean | p |
|---|---|---|---|
| M1 primary | −0.044 | +0.016 | < 0.001 |
| M1 peritoneal | −0.092 | −0.001 | < 0.001 |
| M1 liver | −0.057 | −0.001 | < 0.001 |
| M2 primary | −0.032 | −0.041 | 0.71 |
| M2 peritoneal | −0.001 | −0.002 | 0.52 |
| M2 liver | +0.034 | 0.000 | 0.98 |

Closer relatives have more similar EMT states in mouse 1 at all three sites, confirming state heritability.
Mouse 2 shows none, consistent with its much shallower trees (most leaves at a single edit depth). The uncorrected
Spearman p-values in the first analysis pass were invalid because pairs are not independent; only the permutation
results above are used.

## H3 — clone propensity transfers from primary to metastasis
Only one contrast had enough shared clones (mouse 1, peritoneal, 4 clones with ≥ 20 cells in both compartments):
r = 0.986, parametric p = 0.014, **clone-permutation p = 0.082**. Suggestive, not significant; with four clones a
permutation test cannot go much below this. Underpowered.

## H4 — the strict test for heritable switching propensity
| Mouse / site | clones | plateau B | 95% CI | step-1 decay rate k |
|---|---|---|---|---|
| M1 primary | 9 | 0.000 | [0.000, 0.044] | 0.02 (at bound) |
| M1 peritoneal | 5 | 0.000 | [0.000, 0.034] | 0.02 (at bound) |
| M1 liver | 4 | 0.108 | [0.000, 0.169] | 0.02 (at bound) |
| **M2 primary** | **19** | **0.042** | **[0.005, 0.081]** | **0.33 (resolved)** |
| M2 peritoneal | 5 | 0.000 | [0.000, 0.067] | 0.02 (at bound) |
| M2 liver | 5 | 0.046 | [0.000, 0.128] | 0.02 (at bound) |

**In 7 of 8 mouse-by-site combinations the decay rate hit its lower bound**, which by Proposition T3 means A and B
are not separately identified there: the plateau fraction had a bootstrap CI of [0, 1]. This is counterexample C7
occurring in real data, and it is a property of the measurement, not of the estimator.

The single informative case is **mouse 2's primary tumour**, the only compartment with enough clones (19) and
enough depth variation: B = 0.042, CI [0.005, 0.081]. That is positive evidence for heritable switching propensity
in one compartment of one animal, with n = 2 animals.

## Model comparison on the real trees (M0–M4)
Two-state CTMC fitted by Felsenstein pruning on within-clone trees (7 trees / 8,747 leaves in mouse 1; 19 trees /
13,294 leaves in mouse 2):

| Model | ΔAIC mouse 1 | ΔAIC mouse 2 |
|---|---|---|
| M0 single state | 250.5 | 439.9 |
| **M1 clone-specific propensity** | **0.0** | 40.3 |
| M2 organ induction only | 249.9 | 456.2 |
| **M4 clone + organ** | 30.9 | **0.0** |

- Clone-specific propensity is decisive in both animals.
- **Organ induction on its own (M2) is worth essentially nothing** — 0.6 AIC better than M0 in mouse 1, and 16
  *worse* in mouse 2. Induction only helps once clone differences are in the model (M4 in mouse 2).
- The fitted metastatic offset is **negative** in both mice (−1.04 and −0.39 on the logit scale), i.e. metastatic
  terminal branches move toward the epithelial state, agreeing in sign with the within-clone residual.

**Held-out site prediction.** Fitting without any peritoneal cells and predicting each clone's peritoneal
high-EMT fraction: mean absolute error falls from **0.288 to 0.182** (mouse 1) and **0.209 to 0.089** (mouse 2)
when clone-specific propensity is included. Clone identity generalises to an unseen metastatic site.

**But by counterexample C8 this is not proof of Mechanism B.** M1's clone-specific fractions mix heritable
propensity, frozen state memory and founder sampling. The decomposition, which separates them, supports heritable
propensity in only one of eight compartments. Both results are reported, and they disagree; the decomposition is
the stricter test and the weaker claim is the defensible one.

## Exploratory
Across clones, mean EMT correlates with the S100 programme (the authors' own highlighted gene family) and with
hypoxia and glycolysis. These are exploratory, use the same cells that define the states, and are not corrected
for multiple testing: `results/tables/m2_real_programmes.csv`.

## Limitations
1. **n = 2 animals**, one tumour model, one laboratory.
2. **Clone sizes are extremely skewed** — one clone holds 9,629 of mouse 1's cells and 11,003 of mouse 2's — so the
   effective number of independent clonal units is far smaller than the nominal count.
3. **Compartments differ technically**, which rules out interpreting within-clone shifts.
4. **Tree resolution is shallow**, which is why the identification condition fails in 7 of 8 compartments.
5. **The primary tumour is sampled at endpoint, not at dissemination** — the Milestone 1 assumption that the
   observed primary represents the disseminating pool is untested here too.
6. Lineage reconstruction uncertainty is taken as given: the authors' published trees and clone calls are used
   unchanged, and their error is not propagated.
