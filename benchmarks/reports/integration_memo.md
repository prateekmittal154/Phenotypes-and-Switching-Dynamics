# Integration memo — Cell Systems manuscript

For the authors preparing `Study_D_Cell_Systems_Submission/`. This benchmark package is
separate and **nothing in the manuscript has been edited**. Below is exactly what should
change, with filenames and section references.

---

## 1. The strongest new finding

**Three independently written implementations of the same estimand disagree about the
switching rate by more than an order of magnitude on identical trees, and they bracket the
truth from both sides.** PATH underestimates by 1.6–2.3x, PASTRI overestimates by 11–16x,
Study D overestimates by 1.3–1.6x (Table B3, Figure B2B).

This is a stronger argument than the manuscript currently makes. At present the
non-identifiability of switching kinetics is argued from our own estimator's coverage and
from an information calculation. It can now be argued from the behaviour of other groups'
software: if the rate were identified, three correct implementations would not disagree
this much.

**The supporting finding, nearly as useful:** replacing true branch times with unit edges
shifts the rate by +1.37 (Study D), +1.56 (PATH) and +1.75 (SCOUT) log units while leaving
occupancy bit-identical for 96–100% of clones (Table B4, Figure B4). PASTRI's rate is
**exactly invariant** to branch times — 505 of 505 clones unchanged — because its depth is
topological. The manuscript's rate–time confounding claim (Proposition V4) is thereby
reproduced in three external packages and shown to be structural in a fourth.

---

## 2. The manuscript claim this changes

**`manuscript/main.tex`, Discussion, the "Two boundaries should be stated plainly"
paragraph.** It currently reads:

> We compared our estimator against competing approaches by estimand, and reimplemented a
> PATH-style statistic from its published definition \cite{schiffman2024}, but did not run
> the original PATH, SCOUT or PASTRI implementations \cite{stuart2025,yang2026pastri}; a
> benchmark against the authors' software remains to be done.

**This is now false and must be removed.** All three were installed from their authors'
repositories, reproduced on the authors' own examples, and run. Replace it with a short
statement of the benchmark result and a pointer to the new figure.

Suggested replacement, to be edited to taste:

> Running the original PATH, SCOUT and PASTRI implementations on the same trees makes the
> point more strongly than our estimator can alone. On simulated data with known rates the
> three recover clone occupancy as well as we do, but their switching-rate estimates
> disagree with ours and with each other by more than an order of magnitude, bracketing
> the truth from both sides (Figure 7). Replacing true branch times with unit edges moves
> every rate estimate by about four-fold while leaving occupancy unchanged, and PASTRI's
> rates are exactly invariant to branch times because its depth is topological. The
> asymmetry this paper reports is therefore a property of the data, not of our likelihood.

**Also revise `supplementary/supplementary.tex`, Section 5 "Scope of the claims", item 2**,
which currently says the original implementations were not installed or run. It should now
record what *was* run, and carry forward the remaining limits (see §6 below).

---

## 3. Which benchmark figure belongs in the main paper

**Recommendation: one new main figure, as Figure 7, built from Figure B2B and Figure B4A–B.**

The three panels that earn main-text space:

| Panel | Source | Why |
|---|---|---|
| A | Figure B2A | all four methods recover occupancy — the positive control |
| B | Figure B2B | the three-way rate disagreement — the headline |
| C | Figure B4A + B4B | the unit-edge penalty, and PASTRI's exact invariance |

An alternative that avoids adding a figure is to append B2B and B4A as panels F and G of
the existing **Figure 5** (`figures/Figure_5.pdf`, "Sparse sampling, not reconstruction
error, removes most of the information"), whose subject is exactly this. That figure
already has five panels and would become crowded; a separate Figure 7 reads better.

Source data for either option is in `results/metrics/clone_level_all.csv`, and the figure
script is `src/metrics/make_bench_figures.py`.

---

## 4. What belongs in supplementary information

| Item | File | Placement |
|---|---|---|
| Estimand compatibility matrix | `Figure_B1_estimand_matrix.pdf`, `TableB1_estimand_matrix.csv` | new supplementary figure, cited from Methods |
| Information and identifiability | `Figure_B3_information_identifiability.pdf` | new supplementary figure |
| Clock and reconstruction sensitivity, full | `Figure_B4_clock_sensitivity.pdf` | new supplementary figure |
| Real-data cross-method comparison | `Figure_B5_real_data.pdf` | new supplementary figure |
| Software and commits | `TableB2_software.csv` | supplementary table |
| Clone-level accuracy by mechanism | `TableB3_clone_level_true_tree.csv` | supplementary table |
| Variance components | `TableB3b_variance_components.csv` | supplementary table |
| Tree-geometry results | `TableB3c_geometry.csv`, `TableB3d_window_predicts.csv` | supplementary table |
| Held-out prediction | `TableB3e_holdout_prediction.csv` | supplementary table |
| Clock sensitivity | `TableB4_clock_sensitivity.csv` | supplementary table |
| No-estimate rates | `TableB4b_no_estimate_rates.csv` | supplementary table |
| Reconstruction | `TableB4c_reconstruction.csv`, `TableB4d_reconstruction_stability.csv` | supplementary table |
| Animal-level real data | `TableB5_real_data.csv`, `TableB5b_cross_method_agreement.csv` | supplementary table |
| Contribution summary | `TableB6_scientific_contribution.csv` | supplementary table |
| Method reproduction evidence | `reports/method_reproduction.md` | new supplementary section |

A new supplementary section, **"Comparison with existing methods"**, should hold the
estimand matrix, the reproduction evidence and these figures. It slots naturally after the
current Section 4 ("Simulation study").

---

## 5. Where to cite PATH, SCOUT and PASTRI

The bibliography already has `schiffman2024` (PATH), `stuart2025` (SCOUT) and
`yang2026pastri` (PASTRI). Required changes:

1. **Introduction**, where phylogenetic signal is introduced: already cites
   `schiffman2024`. No change needed.
2. **Introduction**, the sentence on methods that estimate switching rates
   (`\cite{gunnarsson2023,stuart2025}`): extend to `yang2026pastri`, since PASTRI is
   directly a switching-rate method.
3. **Results, new subsection** after "What limits recovery of switching kinetics":
   a short paragraph reporting the benchmark, citing all three.
4. **Discussion**: the replacement text in §2 above.
5. **Methods, new subsection "Comparison with existing methods"**: name the versions and
   commits, state that each was reproduced on its authors' own example, and record the
   PASTRI `fortify` substitution (`reports/deviations.md`, B1) and that SCOUT is
   represented by its OUwie engine rather than the full pipeline (B3).

**Add to the bibliography:** PASTRI is currently cited as a bioRxiv preprint. The DOI is
`10.64898/2026.08.16.745065`. The software repositories should be cited explicitly:
`github.com/landau-lab/PATH` (v1.0, `caffc51`), `github.com/hrstuart/SCOUT` (`c72befe`),
`github.com/yuxiaochen11/PASTRI` (`67aefd9`).

---

## 6. Claims that must be revised rather than removed

1. **"The single largest gap"** language about the missing benchmark, wherever it appears
   in `reports/` and in `submission/external_assets_needed.md` item 10, is now closed and
   should be struck. `external_assets_needed.md` currently lists "PATH, SCOUT and PASTRI
   installations" as a non-blocking gap; delete that entry.

2. **Do not claim Study D estimates rates better.** It does not. It produces the largest
   spurious rate heterogeneity of the three binary methods (`V_log k` of 2.2–3.3 where the
   truth is 0, against PATH 0.21–0.35 and PASTRI 0.20–0.37). The manuscript should say
   that no method recovers rate heterogeneity, that ours over-disperses while PATH and
   PASTRI compress, and that the contribution is the identifiability diagnostic and the
   interval, not a precision gain. Writing it any other way would be contradicted by
   `TableB3b_variance_components.csv` in the same package.

3. **Keep the scope statement in Supplementary §5, item 2, but rewrite it.** The original
   implementations have now been run; what remains unexercised is SCOUT's full pipeline
   (EM, lambda screening, hybrid selection), PASTRI's permutation test, and a GESTALT-like
   recorder in the reconstruction benchmark. Those are listed in `reports/deviations.md`
   under "Not done".

---

## 7. One sentence for the abstract, if space allows

> Three independently developed methods, run on the same trees, recover clone occupancy
> alike but disagree about switching rates by more than an order of magnitude.

This belongs after the sentence reporting that occupancy heterogeneity replicates across
five animals, and before the design recommendation.
