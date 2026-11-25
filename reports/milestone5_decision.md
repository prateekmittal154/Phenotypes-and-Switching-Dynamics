# Study D, Milestone 5: decision

**Date:** 2026-11-24. Deviations D22–D26 logged. Milestones 1–4 preserved; the outstanding Milestone 4
calibration completed untouched and is incorporated.

## Decision: **GO** — the reconstruction question is answered, and the answer is negative for algorithmic fixes

The central question was how much switching information is lost to sampling, to recording and to
reconstruction, and whether uncertainty-aware inference can recover it. All three components are now measured.
The answer changes the project's recommendation but not its conclusions.

| Success criterion from the brief | Verdict |
|---|---|
| 1. Quantitative decomposition of information loss | **Met.** Branch-length miscalibration costs 1.83× in $\hat k$; topology error adds only 1.04–1.10× for NJ, UPGMA or an ensemble (2.74× for the simplified greedy). Reconstruction roughly doubles-to-triples the spread of the error without shifting its centre. |
| 2. Reconstruction-aware inference with better calibration | **Partly.** The ensemble beats the average single tree in 69.6% of replicates and cuts mean absolute log error by 18%, but leaves median bias unchanged. Interval coverage under reconstruction error was not measured (D26). |
| 3. Validated explanation of why reconstructed lineages fail to identify rapid switching | **Met, and it is the main result.** The true sampled genealogy already contains only 1.6–1.9% close-relative pairs; reconstructions give 1.2–3.2%. The missing information was never sampled. |
| 4. A more effective design at fixed cost | **Partly.** Dense within-clade sampling is identified as the only lever that addresses the missing pairs, and costs nothing in sequencing; no credible cost-per-precision figure is given, because the pair-based power calculation is unreliable. |
| 5. Occupancy robust across reconstruction strategies | **Met, in both systems.** Quinn: $V_\pi$ = 0.0083–0.0094 across three *published* reconstructions whose heights differ twofold. Simeonov: 0.0464–0.0467 and 0.0563–0.0564 across three reconstructions built from the authors' allele records. |

## The three findings that matter

**1. Star-likeness is biological, not an inference artefact.** This was the open question from Milestone 4 and
the most damaging unanswered criticism of the manuscript. Simulating the full pipeline against known
genealogies shows the scarcity of informative close-relative pairs is already present in the true sampled tree.
Sampling ~100 cells from a large expanding population yields mostly distant relatives. **No reconstruction
method can recover pairs that are absent from the sample**, so Milestone 4's 0.70-decade penalty stands and is
attributable to sampling.

**2. The clock assumption costs more than topology error.** Fitting on the *true* topology with unit edges
inflates $\hat k$ by 1.83×. Reconstruction adds little beyond that. The field's effort goes into tree search;
this says the larger return is in calibrating resolved branches to time.

**3. Occupancy is reconstruction-proof; rate heterogeneity is not.** Across three published Quinn
reconstructions, $V_\pi$ moves ±6% while $V_{\log k}$ spans 6.5-fold (0.28–1.81), with the ensemble smallest.
In Simeonov the same split appears: $V_\pi$ identical to three decimals, while median $\hat k$ ranges 1.8–19.8
in mouse 2 depending on algorithm. **The Milestone 4 conclusion that switching rates are unidentifiable
survives reconstruction uncertainty and is strengthened by it.**

## What this does not establish
- **Joint barcode-likelihood inference was not implemented** (D24). Proposition V1 bounds it by the sampled
  genealogy and finding 1 shows that bound is the binding one, but it would still improve precision and its
  absence is a real gap in the method comparison.
- **Interval coverage under reconstruction error is unmeasured** (D26). Given that reconstruction doubles the
  spread, fixed-tree intervals very likely under-cover; by how much is unknown.
- **The recorder simulation is not a model of either published chemistry.** Synchronous divisions, independent
  per-site editing, an edit rate not estimated from data, and reimplemented reconstruction algorithms. It
  supports relative statements only.
- **The greedy result is not a verdict on Cassiopeia.** Our greedy lacks the priors and refinement of the
  published solver.
- **Real-data ensembles rest on one Quinn animal** (D25); the other three deposits contain a single tree.

## Effect on the manuscript
The central argument is now the one the brief anticipated: recoverable information depends on the genealogy
actually sampled, the information the recorder preserves, and the uncertainty reconstruction introduces. A new
section presents the information-processing chain and the three findings above, with two new figures.
Idealised capacity claims are qualified by the empirical geometry throughout. The compiled PDF is 10 pages.

## Next experiments, in priority order
1. **Interval coverage under reconstruction error** — the cheapest remaining gap, and it bears directly on
   every published switching-rate interval.
2. **A clock-calibration arm in simulation**: quantify how much of the 1.83× bias an independently measured
   edit rate removes. This is the highest-value design claim and is currently argued rather than measured.
3. **Dense-within-clade versus broad sampling at equal budget**, to put a number on finding 1's recommendation.
4. **Joint barcode-likelihood inference** for one small system, to confirm it approaches the sampled-genealogy
   ceiling rather than exceeding it.
5. Only then: benchmark against the original PATH, PASTRI and SCOUT implementations, which remains the main
   obstacle to a methods-journal submission.

## Separation of claims
- **Mathematical:** Propositions V1–V3, proved; V1 is the data-processing inequality applied to a
  parameter-independent chain, not a new theorem.
- **Simulation:** E17, 480 replicates with known ground truth, including a documented and corrected bug (D23).
- **Real-data observation:** reconstruction stability of occupancy in two systems; instability of rate
  heterogeneity.
- **Design recommendation:** argued from the simulation, with cost figures deliberately withheld.
- **Not claimed:** that any reconstruction is biologically correct; that better algorithms would recover the
  rates; any human relevance.
