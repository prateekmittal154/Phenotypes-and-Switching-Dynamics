# Milestone 2: Simeonov (GSE173958 + Mendeley t98pjcd7t6) pre-registration

**Frozen before any EMT state was computed or any metastasis-vs-primary contrast was examined.** Only the metadata
(tissue counts, clone sizes, coarse-cluster counts) had been viewed.

## Data and units
- **Cells:** the authors' cancer-cell table (`mouseN_cds_colData.txt`), restricted to cells with a cloneID. The tables
  hold 12,657 cells for mouse M1 and 15,371 for mouse M2.
- **Clones:** static-barcode clones (cloneID).
- **Within-clone trees:** `pdac_mouseN_full_edgelist.txt`, from the evolving barcodes.
- **Expression:** GEO matrices, joined on the sample prefix and the 10x barcode.

**Independent units.**
- Mouse: n = 2, the only biological replicates.
- Clone: within a mouse. Clones share a host, so they are *not* independent animals.
- Cells are never treated as units.
- All CIs are clone-bootstrap within mouse, and results are reported per mouse. A pooled estimate is only a
  clone-weighted descriptive.

## State definition (state gene set S; frozen)
- **E:** Cdh1, Epcam, Krt8, Krt18, Krt19, Krt7, Cldn4, Cldn7, Esrp1, Ocln.
- **M:** Vim, Zeb1, Zeb2, Snai2, Twist1, Cdh2, Fn1, Sparc, Tgfbi, Col1a1.
- **Expression scale:** x = log1p(1e4 · count / nCount_RNA).
- **Per-gene z-scores** use the mean and SD of **primary-tumour cells of the same mouse**.
- **EMT score:** e = mean_z(M) − mean_z(E).
- **States:** 4 states E, EH, LH and M, cut at the PT-cell quartiles of e for that mouse. The 2-state model uses
  "high" (e > PT median) vs "low".
- **Sensitivity analysis:** 3 states at PT tertiles.
- **Validation of the score:** Spearman correlation of e with the authors' M1 EMT pseudotime. This is **not** used to
  tune anything.

## Held-out gene programmes (V; never used for states)
| Programme | Genes |
|---|---|
| Proliferation | Mki67, Top2a, Mcm2, Pcna, Ccnb1, Cdk1 |
| Glycolysis | Eno1, Pgk1, Ldha, Pkm, Tpi1, Aldoa, Gapdh, Slc2a1, Hk2 |
| Hypoxia | Vegfa, Bnip3, Ndrg1, Car9, Egln3, P4ha1, Adm |
| S100 | S100a4, S100a6, S100a10, S100a11 |
| Stress / AP-1 | Fos, Jun, Atf3, Egr1, Hspa1a, Hspa1b |
| Stemness / progenitor | Cd44, Sox9, Prom1, Aldh1a1, Lgr5 |
| Interferon | Isg15, Ifit1, Ifit3, Irf7, Stat1 |
| Basement membrane / adhesion | Lamb3, Lamc2, Itgb4, Itga6, Col17a1 |
| Chromatin | Hmga2, Hmga1, Ezh2, Kdm5b |

## Confirmatory hypotheses (α = 0.05, two-sided unless stated; reported per mouse)
- **H1, reproduction.**
  - Each metastatic site (peritoneal, liver, and lung where n ≥ 50) has a higher mean e than the PT.
  - This holds both unstratified and as a **within-clone** shift (clone-weighted mean difference, clone bootstrap).
  - Clone aggression (adj_aggro) is related to the clone's PT-cell EMT profile: Spearman correlation between the
    clone's PT fraction in LH and aggression, clones with ≥ 10 PT cells.
- **H2, state heritability (replicates PATH-type findings).**
  - Within clone and site, |e_i − e_j| is smaller for cell pairs whose tree MRCA lies deeper in the evolving-barcode
    tree (more shared edits) than for pairs whose MRCA is the clone root.
  - Test: clone-stratified permutation.
- **H3, heritable propensity across sites (new).**
  - Across clones with ≥ 20 cells in the PT **and** ≥ 20 in a metastatic site, the clone's PT mean e correlates
    positively with its metastatic-site mean e, after removing the site effect.
  - Test: Pearson correlation, clone-permutation p, and clone-bootstrap CI.
- **H4, plateau of the within-clone covariance by MRCA depth (interpreted only under the clock-like editing
  assumption, per theory T2).**
  - C(h) is the covariance of the 2-state indicators of same-site cell pairs whose MRCA is at edit depth h.
  - Fit C(h) = A·exp(−kh) + B by nonlinear least squares, weighted by pair counts.
  - B > 0 is evidence for heritable propensity *conditional on the depth calibration*.
- **H5, erasure.**
  - In metastatic sites, the state autocorrelation between cells with MRCA at the clone root vs deeper is compared
    with the PT.
  - If it decays to its clone-level value within the shallowest depth class, founder-state imprint is short-lived,
    and the metastatic enrichment then requires post-founding dynamics (theory T3).

## Model comparison on real trees (pre-specified)
- **Method:** a 2-state CTMC on each within-clone tree (unit edge lengths), Felsenstein pruning, leaves in a given
  site.
- **Models:**
  - M0: common (a, c) across clones, one site dynamic.
  - M1: clone-specific stationary fraction, common total rate.
  - M2: site-specific rates on terminal edges.
  - M4: M1 + M2.
- **Model choice:** AIC plus parametric-bootstrap LRT (B = 200).
- **Predictive check:** fit on PT + liver, then predict the peritoneal-site per-clone fractions (held-out site).
- **M3** (differential growth / bottleneck) is not expressible in a fixed-tree CTMC likelihood. It is assessed by
  simulation (synthetic experiments) and site-occupancy statistics.

## Exploratory (labelled as such)
- Held-out programme scores vs clone propensity.
- Coarse clusters.
- CTC and surgical-site (M1 only) samples.
- 3-state sensitivity.
