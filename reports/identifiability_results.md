# Identifiability results (Milestone 1, Tasks 2–4 and 6)

**Status labels:**
- **[M]** mathematical result, proved in `theory/proofs/`;
- **[Sim]** simulation evidence;
- **[A]** unverifiable assumption.

## Formal model and estimand
The model is (P, w, K) → Q = πK, with π ∝ Pw. The estimand is τ = γ(X ≠ Y): the fraction of metastatic cells whose state
differs from their founder's state. Full definitions are in `theory/definitions.md`.

**What the starting equation conflates.**
- (i) Dissemination × survival × colonisation × expected expansion are multiplied into w. They are not separable from
  (P, Q), and they do not need to be separated for τ.
- (ii) Post-founding growth of *induced* cells is part of K, not w. τ counts cells, not events.
- (iii) Finite founders break the mean-field equation for a single metastasis (P8).
- (iv) The primary is observed at surgery, not at dissemination, and from a region, not the whole tumour (P9) [A].

## Results

| # | Result | Status | One-line content |
|---|---|---|---|
| P1 | Observational equivalence | [M], known in spirit | Pure selection and pure induction reproduce the same (P, Q). Pure induction fits *every* (P, Q). |
| P2 | Sharp τ_min(Λ) | [M] | τ_min(Λ) = min_m max(E_{Λm}(Q ‖ P), m·E_{1/m}(P ‖ Q)). Closed form matches an independent LP to **1.2e-12** (1,000 instances, E1). |
| P3 | Limits | [M] | τ_min(1) = TV (max error 4e-16), τ_min(∞) = SM (1e-16), monotone (max violation 1.7e-11, bisection tolerance). τ_max is vacuous: median 1.0, minimum 0.77 over 200 instances. |
| P4 | Asymmetry | [M] | Cross-sectional data never exclude induction. They exclude pure selection iff Λ < Λ\* = max(Q/P) / min(Q/P). |
| P5 | Finite samples | [M] + [Sim] | The plug-in is invalid: P(false) ≥ (1 − p)^{n_P}. The robust Clopper–Pearson LP bound is valid by construction. |
| P6 | Measurement error | [M] | 1-Lipschitz in Q, (3Λ − 1)/2-Lipschitz in P. |
| P7 | Clone-resolved bound | [M] | Exact clone decomposition. Without labels the same assumptions are vacuous. Simpson reversal (C2); state-pure clones uninformative (C3). |
| P8 | Founder bottleneck | [M] + [Sim] | A single metastasis with few founders defeats any deterministic bound. |
| P9 | Regional sampling ≡ selection | [M] + [Sim] | Bias ρ turns Λ into Λρ. More cells from one region do not help. |
| P10 | Fluctuation test fails | [M] + [Sim] | Founder-level heritable induction ≡ selection across any number of metastases from one primary. |
| P11 | Dynamics | [M], known | Two snapshots never separate the mechanisms. ≥ 3 non-equilibrium timepoints identify A up to cI. Equilibrium loses the split. |
| P12 | Design law | [M] + [Sim] | n_P ≈ ln(1/α′)(Λ(1 − q) + q)/q. Founder recording needs ≈ ln(1/α)/τ cells for any Λ. |

## Nontriviality against established bounds (Task 6)

**What P2 reduces to.**
- At Λ = 1 it equals **TV**; at Λ = ∞ it equals **support mismatch**.
- In between, it is the TV projection of Q onto the density-ratio band {mP ≤ π ≤ ΛmP}. That is a classical object
  (density-ratio classes; Huber–Strassen band geometry†).
- Combined with the **maximal coupling** theorem, it gives P2.
- **Conclusion:** P2 is *not* a new theorem. It is a new, explicit statement for this estimand.

**Frozen family (G3).** Q and P are independent Dirichlet(0.5) draws. τ_min lies ≥ 0.05 inside both reference bounds in:

| Λ | Share of instances |
|---|---|
| 2 | 76% |
| 5 | 87% |
| 10 | 83% |
| 100 | 48% |

**Biologically shaped family (critical review, deviation D2).**
- Setup: P ~ Dirichlet(2) over 8 states, log-normal selection (σ = 1, median Λ_true = 16), and 15% transitions.
- τ_true has a median of 0.13.
- With Λ **correctly specified** (Λ ≥ Λ_true), τ_min > 0 in only 18% of instances at Λ = 5 (n = 49), 4% at Λ = 10
  (n = 272) and **0%** at Λ = 100 (n = 972).
- The bound is non-trivial *as mathematics*. It is **almost never informative** when induction moves cells between
  states that already exist in the primary and selection is moderate.
- It bites mainly when a metastatic state is **absent or extremely rare** in the primary (S2), or when Λ is assumed
  small.

**Discretisation.**
- Merging states can only lower TV-type bounds (data processing), so τ_min is discretisation-dependent.
- Splitting a primary state into rare sub-states can create apparent novelty. States must be fixed before the
  metastatic data are seen [A].

**Hidden rare states.**
- They are handled by P5b. The bound is valid but needs n_P ∝ Λ/q (P12).

**Measurement error.**
- It is handled additively in Q (P6).
- In S7 (ε_b ≤ 0.1 with Λ = 50), the robust bound produced no false positives even without an allowance. The
  misclassified mass is absorbed by permitted selection.
- Primary-side distortion is far more damaging, with Lipschitz constant (3Λ − 1)/2.

## What changes interpretation
- **(i)** "State X pre-exists in the primary and is enriched in the metastasis" (Tagore 2025, CIN-high cells):
  - this is consistent with selection;
  - it is **not evidence against induction** (P4);
  - it requires patient-matched P and Q to even define τ.
- **(ii)** "State Y is absent from the primary, so it was induced":
  - this is valid only if Λ < (q/p_upper)(1 − p_upper)/(1 − q);
  - for q = 0.3, ruling out pure selection at Λ = 100 needs ≥ 1,444 primary cells, and at Λ = 1,000 ≥ 14,409 cells,
    from **representative** regions.
- **(iii)** A shared genotype between primary and metastasis is not evidence of a state transition (C3).
- **(iv)** Several metastases from one primary do not substitute for independent cultures (P10).
