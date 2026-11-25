# Figure-by-figure evidence audit

Every figure proposed for the manuscript, with the claim it is used to support, the evidence actually behind
it, and whether that evidence is sufficient. Figures that do not survive this audit are dropped from the
manuscript rather than reworded.

| Fig | Claim in the manuscript | Source | Evidence class | Sufficient? |
|---|---|---|---|---|
| 1 (`m3_fig15_audit`) | The plateau measures occupancy and is invariant to the switching rate; a common induction shift rescales but cannot create it | Closed-form evaluation of Proposition 2 | **Mathematical identity**, plotted | **Yes.** Panel (a) is a direct evaluation of the formula; panel (b) is a two-line delta-method consequence. No data involved, no inference. |
| 2 (`m3_fig18_separation`) | The per-clone estimator separates occupancy from switching heterogeneity, and fails below the identifiability window while the tests stay calibrated | E10-focus, 20 clones × 10 replicates, four ground truths | **Simulation with known truth** | **Yes, with the caveat shown.** 10 replicates is thin for the rejection rates; the 300-replicate calibration (Fig 5) is the authoritative version and they agree. |
| 3 (`m4_fig22_fisher_audit`) | kL\* ≈ 0.43, mildly occupancy-dependent; estimating π costs ≤6%; the optimum is flat over a ten-fold band | E12, numerical Fisher information over a grid | **Numerical computation** | **Yes.** Deterministic, no Monte Carlo error. Replaces the Milestone 3 value 0.405. |
| 4 (`m4_fig25_capacity_law`) | W = 1 + log₁₀(resolved levels); published recorders cannot reach fast switching | E15b, idealised balanced-tree law | **Algebraic identity given a numerical band** | **Only with Fig 4b.** The idealised law over-predicts real coverage by a median 0.70 decades (E16). The figure must be shown alongside the empirical verification or it misleads. **Action: add the empirical band to the right-hand panel.** |
| 5 (`m4_fig24_calibration`) | χ² reference is adequate for the shared-rate test, mildly anti-conservative for occupancy, and interval coverage collapses in the frozen regime | E14, 300 replicates × 6 scenarios, Wilson intervals | **Simulation with Monte Carlo uncertainty** | **Yes.** This is the authoritative calibration evidence. |
| 6 (`m4_fig23_full_tree_information`) | Treating pairs as independent overstates information by 1–221×, worst for slow rates and deep trees | E13, exact enumeration at 8 leaves + validated Monte Carlo | **Exact computation plus validated Monte Carlo** (6.0% max relative error) | **Yes.** |
| 7 (`m2_fig10_composition_vs_within_clone`) | Metastatic EMT enrichment in the pancreatic model is compositional | Simeonov re-analysis, 2 mice | **Observational, one dataset, n = 2 animals** | **Yes for the compositional term**, which uses primary-tumour expression only and so cannot be produced by cross-tissue technical effects. **Not for the within-clone bars**, which Fig 8 shows are confounded. The caption must say so. |
| 8 (`m4_fig26_quinn`) | Clone occupancy is heterogeneous in four independent animals and survives depth adjustment; within-clone shifts track sequencing depth | Quinn re-analysis, 4 mice, two state definitions | **Observational, independent system** | **Yes.** Panel 3 is the negative control and is load-bearing, not decorative. |
| — (`m3_fig16_optimal_depth`) | Information collapses onto one function of kL | E11 | Numerical | **Superseded by Fig 3**, which audits the same quantity properly. **Dropped.** |
| — (`m3_fig17_identifiability_windows`) | Recorders span 0.09–0.24 decades | E12 predecessor | Numerical, **wrong criterion** | **Dropped.** The claim is superseded (D17). |
| — (`m3_fig19_real_perclone`) | Per-clone occupancy estimates in Simeonov | M3 | Observational | **Supplement only.** Informative but not load-bearing. |
| — (`m3_fig20_heldout_ladder`) | Clone identity predicts a held-out metastatic site | M3 | Observational, 2 mice | **Supplement only.** The ordering differs between the two animals, so it cannot carry a main claim. |
| — (`m3_fig21_equal_budget`) | Occupancy and switching want opposite designs | E11 | Simulation | **Supplement.** Superseded in the main text by the capacity law, which is the sharper statement. |

## Actions taken
1. **Fig 4 revised** to show the empirical covered band next to the idealised law.
2. **Figs on the 0.405 optimum and the 0.09–0.24-decade windows removed** from the manuscript.
3. **Fig 7 caption** restricted to the compositional term.
4. Five figures moved to supplement; the main text carries eight.

## Claims in the manuscript with no figure
- Proposition 3 (identifiability window) — stated and proved in text; Fig 2 and Fig 4 show its consequences.
- The Quinn 30k animal's ineligibility — a count, reported in text and in
  `milestone4_dataset_expansion.md`.
- The competing-method comparison — a table, not a figure, because the estimands differ and a performance plot
  would imply a benchmark we did not run.
