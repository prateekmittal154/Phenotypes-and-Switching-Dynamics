# Statistical analysis plan (Milestone 1)

1. **Inference target.** τ is partially identified. The identified set under Λ-bounded selection is
   [τ_min(Λ), τ_max(Λ)]. Milestone 1 reports τ_min(Λ) as a **sensitivity curve over Λ**, never as one point estimate.
2. **Finite-sample validity.**
   - τ_L(Λ) is a one-sided lower confidence bound for τ_min(Λ), and therefore for τ.
   - It is valid when the Clopper–Pearson boxes cover P and Q simultaneously (probability ≥ 1 − α by Bonferroni).
   - No bootstrap is used for the confirmatory bound: the bootstrap is unreliable at the boundary τ_min = 0 and for
     unseen states.
3. **Separation of claims.** Every report labels each statement as one of: mathematical result, simulation evidence,
   human-data observation, biological interpretation, or unverifiable assumption.
4. **No Bayesian claims.** No prior is placed on Λ in Milestone 1. Posterior concentration is not identification.
5. **Human data.**
   - Milestone 1 is a feasibility audit only.
   - No primary-vs-metastasis bound is computed on unpaired patients. τ is defined per patient, so pooling unpaired
     primaries and metastases changes the estimand.
