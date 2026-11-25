# Milestone 5: reconstruction-aware inference and its calibration

Three inferential strategies were compared, as the brief sets out:

- **A. Plug-in** — estimate one tree, condition on it as if known. What every published analysis does.
- **B. Tree ensemble** — average the likelihood over several plausible reconstructions.
- **C. Joint / integrated** — infer dynamics from the barcode likelihood without committing to a tree.

## 1. What was built, and what was not

A and B are implemented and run, on both synthetic (E17) and real (E18, E19) data. **C was not implemented.**
That is a deliberate scoping decision rather than an oversight, and the reason is Proposition V1: the joint
approach is bounded above by the information in $(B, X_S)$, which is itself bounded by the true sampled
genealogy. Since Result 2 of `milestone5_information_loss.md` shows the informative pairs are **absent from
the sample itself**, a correct joint method cannot recover them. Building one would be a large undertaking
whose ceiling we have already measured. This is stated as a limitation, not as a claim that C is worthless —
it would still beat A and B on precision.

## 2. Ensemble versus plug-in, synthetic ground truth (480 replicates)

| Strategy | median paired log-ratio vs the true-tree fit | IQR | mean \|log error\| |
|---|---|---|---|
| plug-in, greedy | +1.009 | 2.73 | — |
| plug-in, NJ | +0.643 | 1.32 | — |
| plug-in, UPGMA | +0.708 | 1.47 | — |
| average single tree | — | — | **1.663** |
| **ensemble over three trees** | **+0.642** | 1.51 | **1.368** |

The ensemble is closer to the truth than the greedy tree in **70.2%** of replicates, than NJ in 65.6% and than
UPGMA in 64.0%; it beats the average single tree in **69.6%**, reducing the mean absolute log error from 1.663
to 1.368 (−18%).

**But the ensemble does not remove the bias.** Its median paired ratio (+0.642) is essentially identical to
NJ's (+0.643). Ensembling averages away *which tree you happened to pick*; it cannot touch the clock
miscalibration that every tree shares. This is the practical reading of Proposition V2: a more stable estimate
is not a less biased one.

Bias depends strongly on the true rate:

| true k | L2 (clock error only) | NJ plug-in | ensemble |
|---|---|---|---|
| 0.3 | +0.719 | +0.851 | +0.850 |
| 1.0 | +0.355 | +0.443 | +0.373 |
| 2.0 | +0.290 | +0.247 | **−0.180** |

Slow switching is where the damage is worst — the same regime in which the naive pair-sum overstates
information most (Milestone 4) and in which interval coverage collapses (below). The three failures coincide.

## 3. Ensemble versus plug-in on real data
Quinn, 25 clones of the 5k mouse with three *published* reconstructions
(`milestone5_tree_reconstruction.md`): $V_\pi$ moves by ±6% across trees (0.0083–0.0094) while $V_{\log k}$
ranges 0.28–1.81, a 6.5-fold spread, and the ensemble returns the **smallest** value. Simeonov, with three
reconstructions we built from the authors' allele records: $V_\pi$ is 0.0464–0.0467 (mouse 1) and
0.0563–0.0564 (mouse 2) — identical to three decimal places.

**Occupancy is reconstruction-proof in both systems. Rate heterogeneity is not.**

## 4. Interval calibration
The authoritative interval calibration is the 300-replicate Milestone 4 experiment
(`reports/milestone4_calibration.md`), which establishes, with Wilson intervals:

- the χ² reference for the **shared-occupancy** test is anti-conservative, 0.143 [0.108, 0.188] against a
  nominal 0.05, so the bootstrap threshold of 32.59 must be used instead of 19.68;
- the **shared-rate** test is conservative, 0.013 [0.005, 0.034];
- coverage of the occupancy-variance interval is 0.88–0.99 in normal regimes but **0.253 [0.207, 0.305]** below
  the identifiability window;
- per-clone occupancy intervals cover at 0.91–0.97 except under misspecification, where they fall to 0.719.

Re-checked against the calibrated threshold, the real occupancy results survive everywhere they were claimed:
4 of 4 Simeonov statistics, 3 of 4 Quinn all-sites (the 3-clone animal fails, as expected) and 3 of 3
Quinn primary-only depth-adjusted.

**A gap, stated plainly.** Interval coverage *under reconstruction error* was not measured: E17 reports bias
and spread of point estimates, not coverage of intervals computed on reconstructed trees. Given that
reconstruction roughly doubles the spread without shifting the centre, fixed-tree intervals are very likely
to under-cover, but we have not quantified by how much. This is the first thing to compute next.

## 5. Recommendations
1. **Use an ensemble whenever more than one defensible reconstruction exists.** It is free, it improves
   accuracy in ~70% of cases, and on real data it halves an estimate that would otherwise look like biology.
2. **Do not report fixed-tree intervals for switching rates.** They ignore a variance component that is
   comparable to the signal.
3. **Report occupancy from a single tree without concern** — it is stable to reconstruction in both systems.
4. **Treat any switching-rate estimate in the slow regime as untrustworthy**, where clock bias, information
   inflation and coverage collapse all coincide.
