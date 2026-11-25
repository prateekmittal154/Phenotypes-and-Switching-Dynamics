# Study D, Milestone 1: decision

**Date:** 2026-11-23.
**Protocols:** frozen 17:36:35 +05:30, before any experiment was run (`checksums/FREEZE_SHA256_m1.txt`).
**Deviations:** D1–D4, all logged in `protocols/deviations.md`. None of them changes a gate input.

## Decision: **REFORMULATE**

The frozen rule says: "REFORMULATE if G1 fails but G2–G4 pass".
- G1 failed.
- G2, G3 and G4 passed.
- G7 passed, which rules out the STOP branch "G1 ∧ G6 ∧ G7 fail".

| Gate | Result | Evidence |
|---|---|---|
| G1 Originality | **FAIL** | The central bound P2 is maximal coupling composed with a TV projection onto a density-ratio band. It has the partial-identification logic of the marginal sensitivity model (Tan 2006†; Dorn et al. 2024). It is new only as an explicit statement for this estimand. |
| G2 Validity | PASS | Complete proofs of P1–P12. Closed form vs independent LP: max error 1.2e-12 over 1,000 instances. |
| G3 Nontriviality | PASS on the frozen family | Strictly inside TV and SM (gap ≥ 0.05) in 76%, 87% and 83% of instances at Λ = 2, 5 and 10. **Critical review:** in a biologically shaped family it is positive in only 4% of instances at a correctly specified Λ = 10, and in 0% at Λ = 100. |
| G4 Feasibility | PASS | d = 1,000: closed form 0.4 ms, robust LP 0.5 s |
| G5 Robustness | PASS (within the frozen scenarios) | Robust bound coverage ≥ 0.998; false attribution ≤ 0.002. It **fails outside the frozen scenarios**: under-stated Λ (S8), single metastases with few founders (P8) and regional bias (P9) all give near-certain false claims. |
| G6 Biological relevance | **FAIL** | No patient-matched primary/metastatic single-cell data in the audited human cohorts. Λ is not calibratable. |
| G7 Validation route | PASS | GSE173958 and GSE161363 have DNA-barcode lineage measured independently of expression. Both are public. |

## The 9 required items

**1. The exact unresolved question.** Given primary and metastatic single-cell states, how much of the metastatic
state distribution *must* be due to state transitions? The quantity is τ, the lineage-level fraction of metastatic
cells whose state differs from their founder's state. The question is asked under explicit, quantitative limits on how
strongly selection can favour any pre-existing state. Which measurements remove the dependence on those limits?

**2. Closest competitors.**
- Weinreb 2018: snapshot non-identifiability, no bounds.
- The unbalanced/dynamic OT family (WOT, TIGON, SUDO): the growth penalty *chooses* the selection/transport split.
- Gunnarsson 2023 and Guan 2025: dynamics identifiability.
- PATH 2024, Simeonov 2021, Kim 2018, FateMap: point claims, no sharp sets.
- Marginal-sensitivity-model bounds (Dorn 2024 and others): the same mathematics in causal inference.
- Our difference: an explicit sharp set for τ, a valid finite-sample LCB, and a design calculus (P8–P10, P12). The
  core mathematics is classical. Details: `novelty_and_prior_art.md`.

**3. Assumptions and estimands.** See `theory/definitions.md`.
- τ under M(Λ).
- Mean-field metastasis.
- P measured at dissemination and representatively.
- Clone labels independent of expression (P7).

**4. Proofs and counterexamples.**
- Complete proofs: `theory/proofs/`.
- Explicit counterexamples C1–C5: `theory/counterexamples/C1_to_C5.md`.

**5. Executed simulations.** S1–S8, equivalence twins, design, fluctuation and dynamics, all on the server. See
`simulation_results.md`.

**6. Comparison with established bounds.**
- τ_min interpolates exactly between TV (Λ = 1) and support mismatch (Λ = ∞).
- The upper end of the identified set is vacuous (τ_max median 1.0).
- The naive TV reading and the plug-in sharp bound both falsely claim induction 56–100% of the time under pure
  selection. The robust LCB does so ≤ 0.2% of the time.

**7. Human-data feasibility.** G6 fails.
- The NSCLC atlas is unpaired: 11 primary-tumour and 31 brain-metastasis patients, zero overlap.
- Its genotype data are pooled per sample.
- The other two cohorts are bulk.
- See `human_data_feasibility.md`.

**8. What has and has not been learned.**

*Learned (mathematics and simulation):*
- Cross-sectional data can **never exclude induction**. They exclude pure selection only through an assumed selection
  bound Λ, which no (P, Q) can check.
- The sharp bound is **rarely informative** unless a metastatic state is absent from, or extremely rare in, the
  primary.
- The common analysis readings are invalid:
  - "a novel state ⇒ induction" and "P differs from Q ⇒ induction" claim induction under pure selection with high
    probability;
  - so does the plug-in sharp bound, because unseen rare states are mistaken for absent ones.
- Three design traps:
  - **regional sampling bias multiplies Λ**, so more cells from one region make false claims *more* certain;
  - a **single oligoclonal metastasis** supports no deterministic bound;
  - the **multi-metastasis "fluctuation test" cannot detect heritable founder-level induction**, because all
    metastases share one primary.
- **Clone labels** independent of expression turn a vacuous bound into an informative one (S5: 0 → 0.08 against a true
  0.19), but only if within-clone selection is bounded. Label error erodes power fast.
- **Only founder-linked lineage recording** identifies τ without Λ. It needs ~9–60 metastatic cells for τ = 0.05–0.3.

*Not learned:*
- Any biological fact about human metastasis.
- The realistic magnitude of Λ for transcriptional states.
- Whether the clone-resolved bound is informative on real lineage data.

**9. Recommended next experiment (Milestone 2, to be pre-registered before data inspection).**

*Apply the clone-resolved robust bound P7 to GSE173958 (macsGESTALT).*

Design:
- Primary tumour vs liver, lung and peritoneal metastases, per mouse.
- States: the paper's EMT-continuum bins, fixed *a priori* from the primary tumour only.
- Clone labels: the static barcodes.

Pre-specified question: for which Λ_in does τ_L > 0? In other words, is the published "selection of late-hybrid EMT
states" *forced*, or is within-clone induction of the same magnitude equally compatible?

Outcomes:
- **Positive for the method:** τ_L(Λ_in ≤ 10) > 0 in at least one mouse-site pair.
- **Negative but publishable as a design principle:** τ_L = 0 for all Λ_in ≥ 2. The published selection claim is then
  not distinguishable from induction even with lineage data.

Also:
- Repeat on GSE161363.
- Verify the data availability of Park 2026 (DNA Typewriter). If its data are public, they give the first direct,
  Λ-free check of τ.

**Publication framing if continued.** A methods note or design-principle paper, e.g. *PLoS Comput Biol*, *Cell
Systems* or *Bioinformatics*. It is **not** a new-theorem paper. The prompt's top-tier venues are not justified by the
mathematics.

## Separation of claims
- **Mathematical:** P1–P12. P11 is known. P2 is classical in substance.
- **Simulation:** E1–E7.
- **Human data:** feasibility only. No inference was performed.
- **Biological interpretation:** none claimed.
- **Unverifiable assumptions:**
  - the value of Λ;
  - P at dissemination equalling P at surgery;
  - representative regional sampling;
  - the clone-size law;
  - state discretisation fixed in advance.
