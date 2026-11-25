# Simulation results (Milestone 1, Task 5)

**Run details.**
- Canonical run on the vedas server, 2026-11-23: `experiments/run_all.py`, seed 20261123, R = 400 replicates per cell.
- Tables are in `results/tables/`, metrics in `results/metrics/milestone1_metrics.json`, figures in
  `results/figures/m1_fig1–7`.
- Everything here is **simulation evidence**, not biology.

## Structural tests
**Equivalence twins (E4).**
- For S1, S3 and S4, the induction-only twin has the *same* (P, Q) (max |ΔQ| = 0) but τ = 0.79, 0.74 and 0.81, against
  0, 0.12 and 0.
- Every estimator output was **bit-identical** between twin pairs.
- The robust bound covered both members.

**Counterexamples C1–C5 (E3).**
- All five reproduce numerically (see `theory/counterexamples/C1_to_C5.md`).

## Coverage and false attribution (E4, cells per site n_P = n_Q)
| Scenario (truth) | Analyst Λ | Plug-in τ_min > 0 | TV > 0 | Novel-state mass > 0 | Robust τ_L > 0 | Robust coverage |
|---|---|---|---|---|---|---|
| S1 pure selection, p = 0.005 (τ = 0) | 50 | 0.65–0.77 | 1.00 | 0.43 (n = 200) → 0 | ≤ 0.002 | ≥ 0.998 |
| S4, p = 1e-4 (τ = 0) | 4,285 | 0.60–1.00 | 1.00 | 1.00 → 0.15 (n = 20k) | 0 | 1.00 |
| S4, p = 5e-4 (τ = 0) | 857 | 0.56–0.89 | 1.00 | 0.89 → 0 | ≤ 0.002 | ≥ 0.998 |
| S7 pure selection + misclassification ε_b ≤ 0.1 | 50 | 0.61–0.64 | 1.00 | 0 | 0 | 1.00 |
| S8 pure selection (Λ_true = 50), Λ assumed 5 / 10 / 20 | — | 1.00 | 1.00 | 0 | **1.00 / 1.00 / 0.41** | **0 / 0 / 0.60** |

**G5 rule.** The frozen rule (coverage ≥ 0.93 in S1, S3, S4, S5, S6; false attribution ≤ 0.05 under S1 and S4) is
**met**:
- coverage ≥ 0.998 everywhere;
- false attribution ≤ 0.002.

**What else the table shows.**
- The naive readings are badly anti-conservative. "Difference = induction" (TV) always claims induction. The plug-in
  sharp bound claims it 56–100% of the time under pure selection.
- **S8 is the central weakness.** Under-stating Λ produces near-certain false induction claims, and nothing in (P, Q)
  can detect the misspecification.

## Power and tightness (E4)
| Scenario | τ_true | τ_min (population) | Mean robust τ_L at n = 1k / 2k / 5k / 20k | P(τ_L > 0) |
|---|---|---|---|---|
| S2 pure induction into an absent state | 0.30 | 0.30 | 0.02 / 0.14 / 0.22 / 0.28 | 0.97 / 1 / 1 / 1 |
| S6 reversible switching + growth (Λ_eff = 19.9) | 0.37 | 0.115 | 0.01 / 0.03 / 0.06 / 0.09 | 0.38 / 0.84 / 1 / 1 |
| S3 mixed (Λ = 10) | 0.12 | 0.024 | ≈ 0 / 0 / 0 / 0.004 | 0 / 0 / 0.04 / 0.80 |

Under Λ-bounded selection, even the *population* bound recovers only 20–30% of τ when the induced states already exist
in the primary (S3, S6).

## Clone-resolved data (E4-S5; n = 3,000 cells per site, 3 clones)
- **Truth.** τ_true = 0.190 and the clone-resolved τ_min = 0.174.
- **Population bound.**
  - Valid at Λ = 5,000, since clone selection is unconstrained: it is **0**, i.e. vacuous.
  - Misspecified at Λ = 5: mean 0.077. In the **null twin** (same clones, no induction) it is positive in **100%** of
    replicates. That is false.
- **Clone-resolved robust bound.**
  - Mean 0.082 / 0.033 / 0.014 at clone-label error 0 / 0.05 / 0.2.
  - Coverage 1.00 at every error level.
  - Null twin: 0% false positives.
- **Label noise.** It **erodes power quickly** (0.082 → 0.014 at 20% error). In this design it did not create false
  positives.

## Design calculus (E5)
**Primary cells needed to exclude pure selection** when a metastatic state at frequency q is absent from the primary
(analytic, Bonferroni over 6 states):

| q | Λ = 2 | 10 | 100 | 1,000 | 10,000 |
|---|---|---|---|---|---|
| 0.05 | 238 | 1,177 | 11,734 | 117,306 | 1.17 M |
| 0.10 | 115 | 559 | 5,560 | 55,568 | 556 k |
| 0.30 | 32 | 148 | 1,444 | 14,409 | 144 k |

**Checks of the analytic law.**
- The LP simulation (q = 0.3) agrees: power switches from 0 to 1 between 100–300 cells (Λ = 10), 1k–3k (Λ = 100) and
  10k–30k (Λ = 1,000).
- Regional depletion ρ multiplies every entry by ρ.
- **Founder-linked lineage recording** needs 9 / 29 / 59 metastatic cells for q = τ = 0.3 / 0.1 / 0.05, for **any** Λ.
  This assumes error-free founder-state recording.

**Regional sampling bias.** The truth is pure selection; the primary is sampled from one region where the selected
state is ρ-fold depleted.

| ρ | False induction at n = 2k | at n = 20k |
|---|---|---|
| 3 | 0.06 | **1.00** |
| 10 | 0.61 | **1.00** |

More cells make the wrong answer more certain.

## Founder bottleneck and multi-metastasis designs (E6)
**A single metastasis under neutral, transition-free seeding (τ = 0):**

| N_f | False-induction rate |
|---|---|
| 1–3 | 0.78–1.00 at every assumed Λ |
| 10 | 0.89 at Λ = 10 |
| 100 | 0.95 at Λ = 1 |
| 1,000 | 0.21 at Λ = 1 |

The mean-field bound is only meaningful for polyclonal metastases or pooled metastases.

**Eight metastases from one primary:**
- The between-metastasis variance under selection and under heritable founder-level induction is indistinguishable
  (KS p = 0.9997, 0.58, 0.99 and 0.82 at N_f = 1, 3, 10 and 100). This is P10.
- Per-cell induction gives near-binomial variance (≈ 7e-5, against 3e-3 to 0.21 under selection; p ≈ 0).
- So the fluctuation design detects only non-heritable, per-cell induction.

## Dynamics (E7; 2-state, truth a = 0.08, c = 0.02, Δb = 0.6)

| Design | 10–90% range of the fitted parameters |
|---|---|
| 2 snapshots, noise-free | Converges to (0.049, 0.009, 0.96): a different, equally exact solution |
| 2 snapshots, 2,000 cells per timepoint | Δb 0.03–0.97 |
| 4 snapshots, noise-free | Exact |
| 4 snapshots, 2,000 cells per timepoint | a 0.065–0.097, c 0–0.078, Δb 0.52–0.77 |
| 8 snapshots, 20,000 cells per timepoint | Δb 0.58–0.64 |

c, the back-switching rate, is the practically weakest parameter, consistent with Gunnarsson et al. 2023.

## Computational cost (G4)
- Closed form at d = 1,000: 0.4 ms.
- Robust LP at d = 1,000: 0.50 s.
- The whole experiment suite takes about 80 s on one core.
