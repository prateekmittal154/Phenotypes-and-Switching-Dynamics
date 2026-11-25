# Milestone 2: theory

Full statements in `theory/branching_model.md`, proofs in `theory/proofs/T1_T2.md`, counterexamples in
`theory/counterexamples/C6_C8.md`.

## What was added beyond Milestone 1
Milestone 1 had two state distributions and one scalar, τ, bounded under an unverifiable selection limit Λ.
Milestone 2 uses lineage structure, and the estimand becomes a **decomposition** rather than a bound.

## The result
For within-clone cell pairs at lineage distance L, with states centred at the compartment mean,

  **Cov(X_i, X_j) = E_c[π_c(1 − π_c) e^{−2 k_c L}] + Var_c(π_c) = A(L) + B,**

where π_c is clone c's stationary high-state fraction and k_c its total switching rate.
- **A(L)**, which decays with lineage distance, is memory of the ancestor's *state*.
- **B = Var_c(π_c)**, the asymptote, is heritable switching *propensity* — Mechanism B of the brief.

Plotted against MRCA depth, the covariance decays from close to distant relatives and **flattens at B**.

## Why this is useful
| Claim | Status |
|---|---|
| Organ induction does not contribute to B, once states are centred within compartment (T2.1) | proved; 0/15 false positives in simulation |
| Stable-state selection does not contribute to B either, when the decay is resolved (T2.2) | proved; 0/15 false positives, although the naive clone-variance statistic reads 0.040 |
| A and B are separately identified **iff** the decay is resolvable; at k = 0 they are observationally equivalent (T3) | proved; counterexample C7; realised in 7 of 8 real site-by-mouse combinations |
| Founder bottlenecks inflate B̂ by π(1 − π)/N_eff when founder identity is unresolved (T4) | proved; 0/7 false positives in the primary, but metastases are uninformative |
| Clone labels alone bound B from above but do not estimate it (T5.1, C6) | proved; the naive statistic exceeds the truth in every null scenario |

## The estimator the theory forces
Because of T3, fitting A, k and B jointly to one curve is unstable. The two-step procedure is:
1. **Centre within clone.** This removes B by construction, so the resulting curve identifies the decay rate k.
2. **Centre within compartment, with k fixed.** The curve is then linear in (A, B) and both are non-negative, so a
   non-negative least-squares fit is exact.
Confidence intervals come from a clone bootstrap. Implementation: `src/lineage_inference/ctmc.py::two_step_plateau`
and `src/lineage_inference/propensity.py::decompose`.

## Honest placement
Lemma T1 — the e^{−2kL} covariance of a two-state Markov process on a tree — is classical, and the plateau follows
from the law of total covariance. PATH (Schiffman et al. 2024) already measures phylogenetic autocorrelation of
cell states. The contribution here is the *resolution of that autocorrelation as a function of MRCA depth*, the
identification of its asymptote with heritable propensity, the proof that induction and stable-state selection
cannot contribute to it, the sharp identification condition, the founder-bottleneck inflation term, and the
estimator these force. This is a new observation model and estimator built from classical ingredients, which the
brief permits; it is not a new theorem, and it is not presented as one.
