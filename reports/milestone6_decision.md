# Study D, Milestone 6: decision

**Date:** 2026-11-24. Deviations D27–D29 logged. Milestones 1–5 preserved; frozen protocols verify.

## Decision: **GO to submission preparation** — the validation milestone succeeded on its own terms

The brief's success criterion was not a newly identified switching rate but "a mathematically and empirically
validated explanation of when switching rates can be inferred, how uncertainty should be quantified, and what
realistic experimental choices improve that inference." All three are now established.

| Experiment | Outcome |
|---|---|
| 1. Reconstruction-aware interval calibration | **Done.** Fixed-tree inference covers at 0.494 [0.449, 0.538]; reconstruction-aware at 0.675 [0.632, 0.715] pooled, and **0.967–1.000 within the identifiable window** against 0.883 for fixed-tree. Below the window every variant fails (0.03–0.38). |
| 2. Branch-time uncertainty | **Done, with a negative result reported.** A single global scale *worsens* calibration (0.362). Ultrametric times plus ensembling recover the bias to +1.07, matching true-topology inference. No external clock was imposed. |
| 3. Sampling design | **Done.** Two-timepoint harvesting triples informative pairs (0.052 vs 0.017) and recovers about half the gap to an unattainable oracle, at no cost. Balanced clone sampling is worst for kinetics. Coverage gains are suggestive, not established, at this replicate count. |
| 4. Biological robustness | **Done.** Dataset-specific bootstrap: occupancy heterogeneity significant in **8 of 8** eligible combinations across **5 animals**. Negative findings preserved. |
| 5. Manuscript and submission package | **Done.** Compiled PDF, claim-to-evidence audit, journal assessment, submission checklist. |

## The three results that close the project's central question

**1. Inference is sound; the observation process is the problem.** Given true ancestry and true branch times,
profile intervals cover at 0.927 against a nominal 0.95. Everything below that is observational loss, not
estimator failure.

**2. Uncertainty can be quantified honestly, and the honest answer is regime-dependent.** Reconstruction-aware
inference attains near-nominal coverage *within* the identifiable window and fails outside it. That is the
correct behaviour: by the information-processing inequality no method can manufacture coverage where the
recorder wrote nothing. A method reporting nominal coverage below the window would be wrong.

**3. Rate and time are confounded, and the obvious fix is wrong.** Proposition V4: scaling all branches by $c$
and dividing $k$ by $c$ leaves the likelihood unchanged, so unit edges estimate a rate in unknown units. The
natural repair — one global scale factor — makes calibration worse. The valid construction is ultrametric and
uses only what the experiment already knows.

## What the project now concludes, in one sentence
**Phenotypic state occupancy is detectable from lineage-resolved cells and robust to reconstruction choice,
while the dynamics generating those states are not identifiable in any published dataset examined — limited
by the genealogy actually sampled, by temporal calibration, and by recorder capacity, in that order.**

## Honest ledger of what remains
| Gap | Severity |
|---|---|
| Original PATH / PASTRI / SCOUT implementations not run | **Highest.** A referee will ask. Estimands are compared instead, and where they differ this is stated. |
| Sampling-design coverage comparison underpowered (40 observations/strategy; 3.2 h for 5 replicates) | Moderate. The pair-fraction result is solid; the coverage ordering is labelled suggestive. |
| Joint barcode-likelihood inference not implemented | Low. Its ceiling is bounded by Proposition V1 and by the sampling result. |
| Code is scripts, not an installable package; no archived DOI | Blocks submission; mechanical to fix. |
| Author, affiliation, funding and data-availability statements | Blocks submission; mechanical. |
| One tree size and generation count throughout the simulations | Low, but worth a sensitivity sweep. |

## Recommended next actions, in order
1. Package the code with tests, archive, cite the DOI.
2. Fill the administrative sections of the manuscript.
3. Run the 30-replicate sampling-design repeat (~1 day) to firm up the coverage claim.
4. Benchmark against the original competing implementations.
5. Submit to *Genome Biology* (see `manuscript/journal_targets.md`); escalate to *Nature Methods* only after
   (1) and (4).

## Separation of claims
- **Mathematical:** Propositions T1–T2, U1–U6, V1–V4, proved. Lemma T1 and the data-processing inequality are
  classical and attributed.
- **Simulation:** E8–E21, with ground truth and Monte Carlo uncertainty on every error rate, including two
  documented failed corrections (naive global scale; the Milestone 5 RF bug).
- **Observational:** occupancy heterogeneity in 5 animals; reconstruction stability; rate unidentifiability.
- **Not claimed:** any mechanism for why clones differ; human relevance; that better algorithms would recover
  the rates; any gene-level or spatial validation.
