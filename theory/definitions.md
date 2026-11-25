# Definitions (Milestone 1)

**State space.** X is a finite set of malignant-cell states, |X| = d. Examples: clusters, discretised programmes, or
states defined *before* the metastatic data are inspected. Results that depend on the discretisation are flagged.

**Primary distribution P.** P ∈ Δ(X) is the state distribution of malignant cells in the primary tumour *at the time of
dissemination*, over the whole tumour. The sampled primary is P_obs; P_obs ≠ P under regional sampling (Proposition 9).

**Founder selection w.** w : supp P → (0, ∞) is the relative probability that a primary cell in state x disseminates,
survives, seeds and contributes descendants to the sampled metastasis. Dissemination, survival, colonisation and expected
clonal expansion multiply into one founder-level weight. These steps are individually unidentifiable from (P, Q), and
they act identically on the estimand. The selected founder law is π(x) = P(x) w(x) / ⟨P, w⟩.

**Transition kernel K.** K(y | x) is the state distribution of the sampled metastatic descendants of a founder in state
x. It aggregates every transition and every state-dependent growth that happens *after* founding: induction, reversible
plasticity, and post-transition expansion.

**Observation model (mean field).** Q(y) = Σ_x π(x) K(y | x). The joint founder–descendant law is
γ(x, y) = π(x) K(y | x). The model already separates selection (w) from transitions (K). It does not separate
dissemination from survival or growth. Mechanism U (hidden rare states, sampling, batch) enters through P_obs ≠ P,
finite samples, finite founders and measurement error. These are treated separately in Propositions 4–6 and 9.

**Estimand τ (lineage-level transition fraction).** τ = γ(X ≠ Y): the fraction of metastatic cells whose state differs
from the state of the founder they descend from. Its properties:
- τ = 0 is pure selection (Mechanism S).
- τ > 0 means at least that fraction is attributable to transitions (Mechanisms I and M).
- τ counts cells, not transition events. A single early induced founder that expands contributes all its descendants.
- τ is not "the number of switching events". An even number of reversible switches can return a lineage to its founder
  state, so reversible plasticity can be invisible to τ.

**Selection-bound model class M(Λ).** M(Λ) = {(w, K) : max_{supp P} w / min_{supp P} w ≤ Λ}, with Λ ∈ [1, ∞].

**Identified set.** T(Λ; P, Q) = {τ(w, K) : (w, K) ∈ M(Λ), πK = Q}. Its endpoints are
τ_min(Λ) = inf T and τ_max(Λ) = sup T.

**Reference bounds.**
- TV(P, Q) = ½‖P − Q‖₁ is the "difference = transitions" reading, i.e. Λ = 1.
- SM(P, Q) = Q(X \ supp P) is the "novel state = transitions" reading, i.e. Λ = ∞.

**Hockey-stick divergence.** E_γ(Q ‖ P) = Σ_x (Q(x) − γ P(x))₊ (Polyanskiy, Poor and Verdú 2010; differential
privacy). It appears in Proposition 2.

**Clone-resolved data.** A heritable clone label c ∈ {1..C} is measured independently of expression. This gives the joint
laws P(c, x) and Q(c, x).
