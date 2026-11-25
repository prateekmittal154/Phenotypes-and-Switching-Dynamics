# Independent critical review (Milestone 1, Task 8)

An attempt to break the Milestone 1 results. Each row states the attack, the outcome, and where it is documented.

| # | Attack | Outcome |
|---|---|---|
| 1 | "P2 is already known." | **Upheld in substance.** It is maximal coupling (Lindvall†) plus the TV projection onto a density-ratio band (DeRobertis–Hartigan / Huber–Strassen family†), with the partial-identification logic of the marginal sensitivity model (Tan†; Dorn 2024). No retrieved paper states it for the founder→metastasis estimand, but that is an application gap, not mathematical novelty. G1 fails. |
| 2 | "The closed form is wrong." | **Rejected.** An independent LP that uses neither lemma agrees to 1.2e-12 on 1,000 random instances. The limits Λ = 1 → TV and Λ = ∞ → SM hold to 4e-16. |
| 3 | "The bound is vacuous in practice." | **Largely upheld.** Under a correctly specified Λ in a realistic perturbation family, τ_min > 0 in 4% of instances at Λ = 10 and in 0% at Λ = 100, although τ_true ≈ 0.13. The bound is informative mainly for absent or extremely rare primary states. |
| 4 | "Λ does all the work." | **Upheld.** S8: assuming Λ = 5 or 10 when Λ_true = 50 gives 100% false induction claims, and the data cannot reveal this (P4). Any application must report the whole τ_min(Λ) curve and the break-even Λ\*. |
| 5 | "Finite samples break it." | Upheld for the plug-in (56–100% false claims) and **rejected for the robust LCB** (≤ 0.2%), which is valid by construction. |
| 6 | "A single metastasis isn't a mean-field object." | **Upheld.** P8 and E6: with N_f ≤ 10 founders, neutral seeding produces "induction" 89–100% of the time. Applications need polyclonal or pooled metastases, or explicit founder-count modelling. This is a hidden assumption of the prompt's starting equation. |
| 7 | "Sampling the primary from one region is selection in disguise." | **Upheld.** P9 and E5: ρ = 3 gives 100% false claims at 20k cells. More cells do not fix it; more regions do. |
| 8 | "Multiple metastases give a Luria–Delbrück test." | **Rejected for founder-level induction** (P10; KS p ≥ 0.58). It works only for per-cell, non-heritable induction. |
| 9 | "Clone labels solve it." | **Partly.** P7 holds only if within-clone selection is bounded and labels are independent of expression. Power drops by 6× at 20% label error. inferCNV clones from the same RNA do not qualify. |
| 10 | "Discretisation drives the result." | **Upheld as a caveat.** Coarsening lowers τ_min. Splitting rare primary states inflates apparent novelty. States must be frozen before metastatic data are seen. |
| 11 | "Simulation scenarios were tuned." | S5 was re-parametrised after the quick run (D1, logged). The original setting was uninformative (τ_true = 0.009). The G5 coverage criterion is unaffected. No other scenario was changed. |
| 12 | "P at surgery ≠ P at dissemination." | **Unaddressed.** This is equivalent to an extra unbounded reweighting, like P9. It remains an unverifiable assumption for every cross-sectional design, including the mouse lineage data. |

## Net assessment
- The Milestone 1 results are **correct and useful as cautions and as a design calculus**.
- They are not a new theorem.
- The most decision-relevant outputs are the negative ones (attacks 3, 4, 6, 7 and 8). They change how a
  "pre-existing vs induced" claim from primary–metastasis single-cell comparisons should be read and designed.
