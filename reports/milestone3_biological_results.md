# Milestone 3: biological results (Simeonov, GSE173958)

Analysis units unchanged from Milestone 2: 2 mice, clones as the resampling unit, cells never treated as
replicates. The Milestone 2 pre-registration governs the EMT state definition; the per-clone estimator is new and
is reported as such.

## 1. Milestone 2 reproduced
The composition decomposition is unchanged and is not affected by the Milestone 3 correction, because it never
used the plateau: in mouse 1 the clone-composition term (+0.40 to +0.47 EMT units) exceeds the total
primary-to-metastasis difference at all three metastatic sites, and mouse 2 shows no enrichment. The metastatic
EMT enrichment is a change in which clones are present.

## 2. Clone-associated occupancy is robust; clone-specific switching rates are not estimable

Per-clone (π_c, k_c) by maximum likelihood with profile intervals, then noise-corrected variance components:

| Mouse | Scope | clones | V_π (occupancy) | LRT shared π | clones with identified rate | V_log k | LRT shared k |
|---|---|---|---|---|---|---|---|
| 1 | all sites | 11 | **0.0314** (raw 0.0448 − noise 0.0135) | p < 1e-16 | 7 of 11 | 0.000 | p < 1e-10 |
| 1 | primary only | 8 | **0.0526** | p < 1e-12 | 5 of 8 | 0.000 | p < 1e-8 |
| 2 | all sites | 24 | **0.0397** | p < 1e-10 | 11 of 24 | 0.000 | p = 1e-4 |
| 2 | primary only | 18 | **0.0346** | p < 1e-20 | 8 of 18 | 1.197 | **p = 0.96** |

**Occupancy.** Clones differ in equilibrium high-EMT occupancy in both animals and under both scopes, by a wide
margin. Leave-one-clone-out confirms this is not driven by any single clone, including the dominant one
(mouse 1 range 0.019–0.036 across 11 drops; mouse 2 range 0.032–0.043 across 24).

**Switching rate.** The noise-corrected V_log k is zero in three of four analyses: the spread of per-clone rate
estimates is smaller than their own estimation noise. The one positive value (mouse 2 primary, 1.197) is
contradicted by the likelihood-ratio test for a shared rate, **p = 0.96**, so it is attributed to estimation
noise rather than biology.

**Why, and it is the predicted reason.** The identifiability window for these trees is k ∈ [0.17, 0.50] (mouse 1)
and [0.25, 0.50] (mouse 2) per unit branch length. The fitted rates are mostly 1–35, i.e. **above the window**:
17 of 35 clone fits hit a bound. The data say state memory decays faster than one resolvable branch, which is
exactly the failure mode Proposition U4 predicts for trees of maximum depth 3–4. Clone-specific switching
dynamics are therefore **not assessable in this dataset**, and no claim is made about them.

## 3. Predicting a held-out metastatic site
The entire peritoneal compartment was withheld, then each clone's high-EMT fraction there was predicted
(mean absolute error):

| Predictor | mouse 1 | mouse 2 | mouse 1, dominant clone removed | mouse 2, dominant clone removed |
|---|---|---|---|---|
| grand mean | 0.296 | 0.177 | 0.316 | 0.214 |
| primary-tumour clone mean | 0.168 | 0.218 | 0.195 | 0.242 |
| clone mean over all training sites | **0.077** | 0.185 | **0.093** | 0.229 |
| tree-based clone occupancy π̂_c | 0.105 | **0.101** | 0.099 | **0.125** |

**Clone identity carries the signal**; which clone-level summary is best differs by animal. The tree-based
estimate is clearly best in mouse 2 (0.101 against 0.177–0.218) and roughly ties the simple clone mean in
mouse 1. With two animals this does not establish that lineage-depth-aware estimation beats a plain clone mean,
and that is stated rather than finessed. Notably, the primary-tumour clone mean is **worse than the grand mean**
in mouse 2, so "measure the primary and extrapolate" is not reliable.

The dominant clone (3,229 of mouse 1's and 7,650 of mouse 2's held-out cells) does not drive these conclusions:
removing it leaves the ordering intact.

## 4. Independent-gene validation: nothing replicates
Nine held-out programmes, never used to define states, correlated across clones with estimated occupancy, before
and after adjustment for clone size, sequencing depth, mitochondrial fraction and compartment mix:

| Programme | mouse 1 adjusted ρ (p) | mouse 2 adjusted ρ (p) |
|---|---|---|
| glycolysis | −0.26 (0.45) | **+0.82 (<0.001)** |
| hypoxia | −0.05 (0.89) | **+0.82 (<0.001)** |
| stress / AP-1 | **−0.67 (0.023)** | **+0.43 (0.035)** |
| S100 | **−0.83 (0.002)** | +0.37 (0.073) |
| stemness | **−0.65 (0.032)** | +0.16 (0.46) |
| chromatin | +0.60 (0.051) | −0.07 (0.75) |
| proliferation | −0.07 (0.83) | −0.15 (0.49) |

**No programme replicates across the two animals, and three change sign** (stress, S100, stemness). Mouse 2's
hypoxia and glycolysis associations are strong and biologically plausible — a mesenchymal–hypoxic link is
well described — but they are absent in mouse 1. With 18 tests, no multiplicity correction and n = 2 animals,
the honest conclusion is that **these data do not supply independent gene-level validation of the inferred
occupancy differences.** This is reported as a negative result.

## 5. Comparison with the established one-number method
A PATH-style depth-weighted phylogenetic autocorrelation of the state, computed on the same cells:

| Mouse | primary | peritoneal | liver |
|---|---|---|---|
| 1 | **0.524** | 0.090 | 0.098 |
| 2 | 0.109 | 0.004 | 0.007 |

This agrees with the Milestone 2 permutation test — mouse 1 has clear phylogenetic state structure, mouse 2 very
little — and it is a useful summary. What it cannot do is say *what kind* of structure it is. In mouse 1's
primary the autocorrelation is 0.524 and the decomposition attributes a clone-occupancy component of
V_pi = 0.053 to it; the single number cannot separate between-clone occupancy from within-clone state memory,
which is the whole question. The two approaches are complementary, and the decomposition is the one that answers
the mechanistic question.

## 6. Independent model system: Quinn (GSE161363)
The lineage files and metadata were downloaded; the expression matrices were not, so no cell states were defined
and no biological re-analysis was performed there (deviation D13). What the Quinn trees do supply is a
**cross-technology test of the identifiability limit**, and they fail it in the same way: Cassiopeia's windows
span 0.09–0.24 decades, against 0.18 for macsGESTALT. Deeper Cassiopeia trees (max depth 8–14 against 3–4) shift
the window to slower rates without widening it. The limitation is a property of current recorders, not of this
dataset or this estimator.

## 7. What this milestone establishes biologically
1. **Clones differ in equilibrium EMT occupancy**, robustly, in both animals, and this survives removal of any
   single clone.
2. **Metastatic EMT enrichment is clonal composition**, reproduced from Milestone 2 and untouched by the
   correction.
3. **Clone-specific switching dynamics cannot be assessed** from these trees, for a reason the theory predicts
   and quantifies.
4. **No independent gene programme replicates** across the two animals as a correlate of clone occupancy.

The brief asked what to report if the data support only clone-associated state occupancy. That is exactly what
they support, and it is what is reported.
