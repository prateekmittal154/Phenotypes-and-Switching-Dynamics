# Study D, Milestone 3: decision

**Date:** 2026-11-23. Deviations D10–D15 logged in `protocols/deviations.md`; D10 is a correction to Milestone 2
and D15 is cosmetic. All Milestone 1 and 2 outputs preserved; the frozen protocols and the Simeonov
pre-registration still verify.

## Decision: **GO**, on a narrower and better-supported claim than Milestone 2 made

The brief asked not to lead with venue. Judged on what was learned, Milestone 3 clears three of the four bars it
sets, and the honest headline is now **clone-associated state occupancy**, not inherited switching dynamics.

| Standard from the brief | Verdict |
|---|---|
| 1. A more informative and **correctly identified** distinction between inherited state and inherited switching dynamics | **Met.** The two are now separate estimands (V_π, V_log k) with an explicit identifiability window; simulation separates OCC from RATE with 0% false positives on the rate test. |
| 2. A substantially improved inference method with demonstrated advantages | **Partly.** The per-clone estimator does what the pooled plateau cannot, but it fails outside the window (FROZE), and on the real data it beats a plain clone mean in only one of two animals. |
| 3. A biological result supported by independent lineage **and** gene-expression information | **Not met.** The lineage result is solid; no gene programme replicates across the two animals, and three change sign. Reported as a negative result. |
| 4. A mathematically justified experimental-design improvement resolving an identification limitation | **Met, and it is the strongest result.** k·L\* ≈ 0.405; equal-budget trade-off quantified; both published recording technologies shown to miss the window. |

## The correction that reframes the project
Milestone 2 called the plateau "heritable switching propensity". It is invariant to every switching rate. It
measures **clone-associated equilibrium occupancy** (U1). Two further claims were narrowed: immunity to
stable-state selection holds only when clones share occupancy (U2), and a common induction shift cannot *create*
a plateau but does *rescale* one (U3). The Milestone 2 plateau numbers stand as occupancy estimates; the
interpretation placed on them does not. The Simeonov composition finding is untouched, since it never used the
plateau.

> **Superseded in part by Milestone 4.** The window figures below use an ad-hoc criterion
> (1/(2L_max), 1/(2L_min)). On a Fisher-information criterion (>=50% of peak information, k*L in
> [0.10, 1.0]) the same recorders cover **1.1-1.2 decades**, not 0.09-0.24. The qualitative
> conclusion is unchanged -- the fitted Simeonov rates still sit above the covered range -- but the
> quantitative claim was too pessimistic. See `reports/milestone4_design_optimization.md`.

## What is now established

**Mathematics.**
- Joint identifiability of (π_c, k_c) per clone, with the explicit window **1/(2L_max) ≲ k_c ≲ 1/(2L_min)** (U4).
- Fisher information for a lineage pair about log k peaks at **k·L\* ≈ 0.405**, stable over two decades of k and
  across occupancy; information about π comes from cell counts instead (U5). Hence: *occupancy is bought with
  cells, switching rate with pairs about one switching time apart.*

**Method, with its failure mode.** On synthetic ground truth with 20 clones: the rate test has 0% false positives
under both null mechanisms and 100% power under rate heterogeneity. But under frozen states (k below the window)
the variance-component point estimates break badly — V̂_π five-fold too high, V̂_log k spuriously 3.14 against a
truth of 0 — while the likelihood-ratio tests stay calibrated. **Lead with the tests, treat the point estimates
as descriptive.** Magnitudes are under-recovered (11–18% of true V_log k) even when the mechanism is detected.

**Biology (Simeonov, 2 mice).**
- Clones differ robustly in equilibrium EMT occupancy: V_π ≈ 0.031–0.053, LRT p < 1e-10 in both animals and both
  scopes, stable to leaving out any clone including the dominant one.
- Clone-specific switching rates are **not assessable**: 17 of 35 clone fits hit a bound, fitted rates (1–35) sit
  above the window these trees support ([0.17, 0.50]), and the one positive V_log k is contradicted by its own
  LRT (p = 0.96).
- Clone identity predicts a held-out metastatic site; the tree-based estimate wins clearly in mouse 2 and ties a
  plain clone mean in mouse 1.
- The metastatic EMT enrichment remains a change in clonal composition (Milestone 2, unchanged).

**Cross-technology.** macsGESTALT and Cassiopeia trees both give identifiability windows spanning only
**0.09–0.24 decades**. The limitation is a property of current lineage recorders, not of this dataset or this
estimator, and it caps what PATH, PASTRI, SCOUT or this method can extract about dynamics.

## What is not established
- Any claim about inherited switching *dynamics* in real data — the measurement cannot support it.
- Gene-level corroboration: nothing replicates across two animals.
- Calibrated interval coverage for V_π and V_log k: recovery was measured, coverage was not. This is a real gap.
- Quinn was used for tree geometry only; no states were defined there (D13). Park's DNA Typewriter data were not
  located (D14), so the time-resolved question is still open.
- Founder-bottleneck effects were characterised in Milestone 2 and not re-run here under the new estimator.

## Next experiments, in priority order
1. **Interval calibration.** Coverage of the profile and bootstrap intervals for V_π and V_log k across clone
   counts and depth spectra. Cheap, and it is the main methodological gap.
2. **A recorder designed to the window.** Stagger editing, or mix edit rates, so resolvable depths span the
   plausible range of k rather than one narrow band. By U5 this is worth far more than extra cells, and it is a
   concrete, falsifiable prediction the field can test.
3. **Define states in the Quinn data** and repeat the occupancy analysis, for a genuine second system.
4. **Locate the Park DNA Typewriter data**; its timing information is the only route to absolute rates, which a
   fixed-tree likelihood cannot give (T5.4).

## Separation of claims
- **Mathematical:** U1–U6, proved; U1–U3 are corrections, U4–U5 elementary once posed.
- **Simulation:** E10-focus, E11, with known ground truth, including a documented failure case.
- **Real-data observation:** clone occupancy heterogeneity; rate heterogeneity not assessable; held-out
  prediction; depth spectra.
- **Biological interpretation:** clones differ in EMT occupancy in one pancreatic model, n = 2 animals.
- **Not claimed:** inherited switching dynamics; gene-programme correlates; human relevance.
