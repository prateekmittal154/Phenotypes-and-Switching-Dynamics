# Simulation protocol (frozen before execution)

## General settings
- **State space:** finite, d states (d = 6 unless stated).
- **Model:** the primary P, founder selection w, and kernel K give the joint founder–metastasis law
  γ(x, y) = π(x) K(y | x), with π ∝ P·w and Q = πK.
- **Truth:** τ_true = γ(X ≠ Y), and Λ_true = max w / min w over supp P.
- **Seed:** 20261123 for all experiments.
- **Replicates:** R = 400 per cell of each design grid.
- **Cell sampling:** multinomial, n_P primary cells and n_Q metastatic cells (default 2,000 / 2,000; grid
  {200, 1,000, 5,000, 20,000}).

## Scenarios
- **S1, pure selection.** A rare state s has p = 0.005 and w_s / w_other = 50; K = I. Analyst Λ = 50.
- **S2, pure induction.** w ≡ 1. A fraction 0.3 of cells moves into a state s absent from P. Analyst Λ = 50.
- **S3, mixed.** Selection ratio 10 on a rare state plus 15% induction into a neighbouring state. Analyst Λ = 10.
- **S4, rare-state undersampling.** p_s ∈ {1e-4, 5e-4, 2e-3}, w ratio chosen so that Q(s) = 0.3, K = I. Analyst Λ =
  Λ_true. This is the plug-in failure test.
- **S5, clone-specific effects.**
  - 3 clones with unconstrained clone selection (clone weight ratio up to 1,000).
  - Within-clone state selection ratio ≤ 5.
  - Clone 2 has 20% within-clone induction.
  - Clone-label error ε ∈ {0, 0.05, 0.2}.
  - Compare population and clone-resolved bounds.
- **S6, reversible plasticity during expansion.**
  - 2 states with switching rates a, c and growth rates b1, b2 over time t.
  - The founder–descendant joint law is π(x)[exp(t(B+T))](x, y), normalised.
  - Λ_eff includes expansion.
- **S7, confounded measurement.** Pure selection (S1 truth) plus organ-specific state misclassification ε_b ∈ {0.02,
  0.05, 0.1}. Evaluated with and without an ε_b allowance.
- **S8, misspecification.** S1 truth with Λ_true = 50, analysed with Λ ∈ {5, 10, 20}.

## Observational-equivalence pairs (structural tests)
For S1, S3 and S4 truths, the induction-only twin (w ≡ 1, K(y | x) = Q(y)) reproduces (P, Q) exactly. Both members of
each pair are passed through every estimator, and the outputs must be identical.

## Estimators and baselines
1. **Plug-in τ_min(Λ)** from empirical P̂ and Q̂.
2. **Robust lower confidence bound τ_L(Λ):** the LP over simultaneous Clopper–Pearson boxes for P and Q (Bonferroni,
   α = 0.05).
3. **Baseline TV(P̂, Q̂)**, the "difference = induction" reading (the Λ = 1 bound).
4. **Baseline support mismatch** Q̂(states with P̂ = 0), the "novel state = induction" reading.
5. **Clone-resolved τ_L** (S5).

## Metrics
- **Coverage:** P(τ_L ≤ τ_true).
- **False induction attribution:** P(estimate > 0) when τ_true = 0, reported for every estimator.
- **Bound gap:** τ_true − τ_L.
- **Run time.**
- **Exact agreement** between the closed form and an independent LP on 1,000 random instances (Dirichlet(0.5) P and Q,
  d ∈ {3, …, 12}, Λ log-uniform on [1, 1,000]).

## Design experiment (Q5)
Required n_P for P(τ_L > 0) ≥ 0.8, as a function of (q, Λ) under:
- pure induction into an absent state;
- regional sampling bias ρ;
- clone-resolved data;
- founder-linked lineage recording.

## Multi-metastasis fluctuation design
- M = 8 metastases, each from N_f ∈ {1, 3, 10} founders.
- Heritable vs fast-switching states.
- Test the between-metastasis overdispersion of the state-s fraction.
