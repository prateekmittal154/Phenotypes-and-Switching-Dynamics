# The lineage-resolved model and the covariance decomposition (Milestone 2)

This extends Milestone 1, where the data were two state *distributions* (P, Q) and the only estimand was the scalar
τ. Here the data are *lineage-resolved*: each cell carries a clone label, a position in a reconstructed tree, a
compartment and a state. The question becomes which of the five mechanisms A–E of the Milestone 2 brief are
distinguishable.

## 1. The process

**Individuals.** A cell is (c, s, o, ν): clone c ∈ {1..C}, state s ∈ {0, 1} (low / high EMT), compartment
o ∈ {0 primary, 1..O metastatic}, and node ν in the lineage.

**Dynamics.** In compartment o a cell of clone c and state s
- divides at rate β(o, s),
- dies at rate δ(o, s),
- switches 0 → 1 at rate a_c g_o and 1 → 0 at rate b_c ℓ_o,
- if in the primary, disseminates at rate μ κ_s.

A disseminating cell picks compartment o with probability q_o and founds successfully with probability φ_o.

**Mechanism map.**
| Mechanism | Parameter signature |
|---|---|
| A stable-state selection | κ₁ ≫ κ₀, and a_c, b_c small (state nearly frozen) |
| B heritable switching propensity | Var_c(a_c / (a_c + b_c)) > 0 |
| C environmental induction | g_o ≠ 1 or ℓ_o ≠ 1 for o ≥ 1 |
| D differential expansion and bottleneck | β, δ depend on o and s; φ_o small so founders are few |
| E coupled | several of the above |

**Clone-level propensity.** Write π_c = a_c / (a_c + b_c) for the stationary high-state fraction of clone c, and
k_c = a_c + b_c for its total switching rate. (π_c, k_c) is the natural parameterisation: π_c is the *propensity*,
k_c the *speed*. Mechanism B is exactly Var_c(π_c) > 0.

The expectation-level equation in the brief, dn/dt = [diag(r) + Kᵀ]n + u, is the mean of this process. It is used
only for intuition: Proposition T4 shows it is not adequate for a metastasis with few founders.

## 2. The central result

**Lemma T1 (pair covariance under a two-state CTMC).** Let two cells of clone c have their most recent common
ancestor (MRCA) at lineage distance L from each, with the MRCA state drawn from the clone's stationary law and the
two post-MRCA paths conditionally independent given that state. Then

  Cov(X_i, X_j | c, L) = π_c (1 − π_c) e^{−2 k_c L}.

**Proposition T2 (decomposition into state memory and heritable propensity).** Draw a within-clone pair at lineage
distance L, the clone being drawn with weight w_c, and centre states at the compartment mean μ = Σ_c w_c π_c. Then

  **Cov(X_i, X_j) = E_c[π_c (1 − π_c) e^{−2 k_c L}]  +  Var_c(π_c)
                  = A(L)                             +  B.**

- **A(L)** is *state memory*: how long a lineage remembers its ancestor's state. It decays with lineage distance.
- **B = Var_c(π_c)** is *heritable switching propensity*: Mechanism B, and it does **not** decay with distance.

Writing the covariance as a function of MRCA depth h (deeper MRCA = closer relatives = smaller L) gives a curve that
decays from close to distant pairs and **flattens at B**. The plateau is the quantity of interest.

Proofs: `proofs/T1_T2.md`. Numerical verification: E9 (`results/tables/e9_mechanism_recovery.csv`).

## 3. Consequences

**Corollary T2.1 (induction does not inflate B).** Compartment-dependent rates g_o, ℓ_o shift the compartment mean
but are common to all clones within a compartment. Centring within compartment therefore removes them, and
Mechanism C contributes nothing to B. Verified: in E9 the induction scenario gives B̂ = 0 with no false positives.

**Corollary T2.2 (stable-state selection does not inflate B either, provided the decay is resolved).** Mechanism A
is k_c → 0 with common π_c. It puts all weight in A(L) and none in B. Verified: B̂ = 0 in E9, although the naive
statistic Var_c(empirical clone fraction) is 0.028 and would have wrongly indicated heritable propensity.

**Proposition T3 (identification condition, and the exact ambiguity).** A and B are separately identified from the
covariance curve **iff the decay of A(L) is resolvable over the observed range of lineage distances**, that is iff
A > 0 and 2k (L_max − L_min) is bounded away from 0. In the limit k → 0 the curve is flat and

  C(h) = A + B  for every h,

so only the sum is identified: **perfectly heritable state (Mechanism A with k = 0) and heritable propensity
(Mechanism B) are observationally equivalent.** This is the lineage-resolved analogue of the Milestone 1
selection/induction equivalence, and it is sharp — see `counterexamples/C7.md`.

*Estimator consequence (the two-step procedure).* Centring **within clone** removes B by construction, so the
clone-centred curve identifies k alone. Fixing k from that curve makes the compartment-centred curve linear in
(A, B), so B follows by non-negative least squares. Fitting A, k and B jointly to one curve is unstable exactly
because of T3. Implementation: `src/lineage_inference/ctmc.py::two_step_plateau`.

**Proposition T4 (founder bottlenecks inflate B̂ when founders are unresolved).** Suppose the cells of clone c in
compartment o descend from N founders with relative expansion sizes Z₁..Z_N, and that the lineage tree does not
resolve founder identity. The realised high-state fraction of that clone-compartment then has variance
π_c(1 − π_c) / N_eff with N_eff = (Σ Z_i)² / Σ Z_i², the inverse Simpson number of founders. This adds to B̂:

  E[B̂] ≈ Var_c(π_c) + E_c[π_c (1 − π_c) / N_eff,c].

So B must be estimated where founders are many — in practice the primary tumour — or founder identity must be
resolved by the tree. This is the lineage-resolved form of Milestone 1's Proposition P8. Verified: in E9 the
bottleneck scenario gives B̂ = 0 in the primary, while its metastases retain too few clones to estimate B at all.

**Proposition T5 (what each measurement buys).** Collected in `identifiability_extensions.md`:
clone labels alone identify Var_c of the *realised* fraction, which confounds A, B and founder sampling; clone
labels plus a resolved tree identify A and B separately, under T3; compartment labels identify induction only as a
compartment main effect, never as a per-clone effect; a second time point is needed for the absolute rates a_c and
b_c, since a single harvest identifies only π_c and k_c.

## 4. Relation to existing theory
- **Lemma T1 is classical.** It is the two-state Markov model on a tree (Felsenstein 1985 lineage covariance; the
  binary analogue of the Ornstein–Uhlenbeck covariance used for phylogenetic signal).
- **The plateau as a between-group variance** is the law of total covariance, also classical.
- **PATH** (Schiffman et al. 2024) measures phylogenetic autocorrelation of cell states — essentially C at one
  aggregated distance — and reports heritability versus plasticity. It does not separate the decaying and
  non-decaying parts, so a clone-level propensity difference and a slowly switching state are not distinguished.
- **The contribution here** is therefore not Lemma T1. It is: (i) resolving the covariance as a *function* of MRCA
  depth and identifying its asymptote with Var_c(π_c); (ii) Corollaries T2.1–T2.2, which show the plateau is immune
  to induction and to stable-state selection; (iii) the sharp identification condition T3 and the two-step
  estimator it forces; (iv) Proposition T4, the founder-bottleneck inflation. This is "classical ingredient, new
  observation model and estimator", which is what the Milestone 2 brief permits — and it is strictly richer than
  the Milestone 1 scalar bound, which used no lineage information at all.
