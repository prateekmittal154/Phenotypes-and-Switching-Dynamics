# Explicit counterexamples

All five are verified numerically in E3 (`results/metrics/milestone1_metrics.json`, key `E3`).

## C1: selection-only and induction-containing mechanisms with identical data
P = (0.90, 0.09, 0.01) and Q = (0.40, 0.30, 0.30).

| Mechanism | w | K | τ | Λ needed |
|---|---|---|---|---|
| Pure selection | (0.444, 3.33, 30) | I | 0 | 67.5 |
| Pure induction | (1, 1, 1) | rows = Q | 0.61 | 1 |
| Minimal-transition mechanism at Λ = 1 | (1, 1, 1) | maximal coupling | 0.50 = TV | 1 |

The data cannot distinguish these mechanisms. τ_min(Λ) = 0.50, 0.42, 0.26, 0.23, 0.125 and 0 at Λ = 1, 2, 5, 10, 30
and 67.5. τ_max = 1 at Λ = 100.

**Biological reading.** A 30% metastatic state seen at 1% in the primary needs at least 12.5% induced cells only if no
state can be ≥ 30× more metastatically fit than another. Large differences in metastatic efficiency between subpopulations are
classically reported in transplantation assays. Their magnitude for *transcriptional states* in human tumours is
unknown, and it was **not verified in Milestone 1**. That gap is exactly the problem.

## C2: Simpson reversal (clone labels raise the bound from 0)
Two clones and two states, with Λ_in-bounded within-clone selection and unconstrained clone selection.

P(c, x) = [[0.4, 0.1], [0.1, 0.4]] and Q(c, x) = [[0.1, 0.4], [0.4, 0.1]].

**Population.** P = Q = (0.5, 0.5), so τ_min^pop = 0 even at Λ = 1.

**Clone-resolved.** τ_min^clone = 0.467 at Λ_in = 2, 0.30 at Λ_in = 4, and 0 at Λ_in = 16.

The clone labels reveal that each clone flipped its dominant state, which is invisible in aggregate.

## C3: state-pure clones (clone labels add nothing)
P(c, x) = diag(0.7, 0.2, 0.1) and Q(c, x) = diag(0.1, 0.3, 0.6).
- τ_min^clone = 0 even at Λ_in = 1, because clone selection explains everything.
- A shared genotype between primary and metastasis is not evidence against induction.
- A changed clone composition is not evidence for selection *of states*.

## C4: two snapshots (2-state dynamics)
x₀ = 0.05 and x₁ = 0.30 at t = 1.
- Pure selection: Δb = 2.10, a = c = 0.
- Pure induction: a = 0.305, Δb = c = 0.
- Both reproduce x₁ exactly. RK4 gives 0.300 and 0.300.
- In E7, the noise-free 2-snapshot fit converges to (a, c, Δb) = (0.049, 0.009, 0.96) instead of the true
  (0.08, 0.02, 0.60).

## C5: equilibrium
Both parameter sets hold x\* = 0.3 as a fixed point (RK4 verifies 0.300 after t = 50):
- **Selection-dominated:** a = 0.01, c = 0.373, Δb = 0.5;
- **Switching-dominated:** a = 0.3, c = 0.7, Δb = 0.

A metastasis sampled at its steady-state composition carries no information on the selection/switching split. This is
consistent with Fassoni (2026), who shows that fitness is asymptotically irrelevant under uniform competition.
