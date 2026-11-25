# Propositions (Milestone 1)

Every proposition carries one status label:
- **[new-form]**: the statement is in this form for this problem. The tools are classical.
- **[known]**: the result exists in the literature under the name given.
- **[corollary]**: an immediate consequence.

Proofs are in `proofs/`. Numerical checks are in `results/` (E1–E7).

**P1 Observational equivalence [known in spirit; explicit construction].** Let supp Q ⊆ supp P. Pure selection
(w = Q/P, K = I, τ = 0) and pure induction (w ≡ 1, K(y | x) = Q(y), τ = 1 − ⟨P, Q⟩) produce the same (P, Q).
- Pure induction fits **every** (P, Q), including those with novel states.
- Proof: `proofs/P1_P2_P3.md`. Instance: `counterexamples/C1.md`.

**P2 Sharp lower bound under Λ-bounded selection [new-form].** For Λ ∈ [1, ∞),

  τ_min(Λ) = min_{m ∈ [1/Λ, 1]} max{ Σ_x (Q − ΛmP)₊ , Σ_x (mP − Q)₊ }
           = min_{m} max{ E_{Λm}(Q ‖ P), m · E_{1/m}(P ‖ Q) }.

- The minimum is attained, so the bound is sharp.
- T(Λ) is the closed interval [τ_min(Λ), τ_max(Λ)].
- The function is computable by bisection in O(d log(1/ε)).
- Proof: maximal coupling, then the box-projection lemma, then the crossing argument. See `proofs/P1_P2_P3.md`.

**P3 Limits, monotonicity, vacuity [corollary].**
- τ_min(1) = TV(P, Q).
- τ_min is non-increasing in Λ.
- τ_min(Λ) → SM(P, Q) as Λ → ∞.
- τ_min(Λ) = 0 iff supp Q = supp P and max(Q/P) / min(Q/P) ≤ Λ on supp P.
- τ_max(Λ) ≥ 1 − ⟨P, Q⟩ for every Λ. The upper end is uninformative unless K is restricted.

**P4 Asymmetry [corollary].** Cross-sectional (P, Q) data can never exclude induction. They can exclude pure selection
only through an assumed Λ, and only if Λ < Λ*(P, Q) := max(Q/P) / min(Q/P).

**P5 Finite samples [new-form].**
- (a) The plug-in τ_min(Λ; P̂, Q̂) is **not** a valid lower bound. If a pure-selection state s has primary frequency p
  and metastatic frequency q_s, then P(τ̂_min ≥ Q̂(s)) ≥ (1 − p)^{n_P}, even though τ = 0.
- (b) τ_L is a valid level-(1 − α) lower confidence bound for τ under M(Λ). It is the minimum of τ_min(Λ; P′, Q′) over
  simultaneous Clopper–Pearson boxes. τ_L is the value of one linear programme with 4d + 1 variables.

**P6 Measurement error [corollary].** Suppose metastatic states are misclassified by a kernel M with
max_x (1 − M(x, x)) ≤ ε. Then |τ_min(P, QM) − τ_min(P, Q)| ≤ ε. τ_min is 1-Lipschitz in Q (TV) and
(3Λ − 1)/2-Lipschitz in P.

**P7 Clone-resolved sharp bound [new-form].** Assume:
- clone labels are exact and heritable;
- transitions preserve the clone;
- w(c, x) = v(c) r_c(x) with v unconstrained and max r_c / min r_c ≤ Λ_in;
- supp Q(c, ·) ⊆ supp P(c, ·).

Then τ_min^clone = Σ_c Q(c) τ_min(P(· | c), Q(· | c); Λ_in), and the bound is sharp. Corollaries:
- (a) Without clone labels, the same assumptions give only τ_min^pop ≤ SM(P, Q). This is vacuous whenever clones are
  state-pure, since clone selection can then realise any reweighting.
- (b) Simpson reversal: P = Q at the population level is compatible with τ_min^clone > 0.
- (c) If every clone is state-pure, τ_min^clone = 0. Shared genotype then adds nothing.
- (d) Clone-label error is not covered. Misassignment creates spurious within-clone shifts (E4-S5).

**P8 Founder bottleneck [new-form].** Let a single metastasis arise from N_f founders with stochastic clone sizes Z_i
having positive density on (0, ∞). Then for neutral, transition-free seeding (Λ = 1, τ = 0), every Q with supp Q ⊆
supp P and |supp Q| ≤ N_f has positive probability density. Consequences:
- No deterministic bound applied to *one* metastasis controls false induction claims.
- The mean-field model holds only for E[Q] across exchangeable metastases, or as N_f → ∞.

**P9 Regional sampling ≡ selection [corollary].** Suppose the primary is sampled with regional weights v(x) of ratio ≤ ρ,
so P_obs ∝ P v. Then any (w, K) ∈ M(Λ) for (P, Q) corresponds to (w/v, K) ∈ M(Λρ) for (P_obs, Q). Consequences:
- An analysis of P_obs must use Λρ.
- More cells from the same region do not reduce ρ.

**P10 Founder-level induction is invisible to multi-metastasis fluctuation [new-form].** Compare two mechanisms:
- each founder independently converts x → y with probability K(y | x), and the new state is heritable;
- pure selection with founder law π′ = πK.

They induce the **same joint law** of the compositions of any number of metastases from the same primary, for any N_f
and any state-dependent clone-size law. The Luria–Delbrück logic therefore does not transfer. All metastases sample
**one** shared primary, so there are no independent "cultures" in which pre-existing variants could jackpot. Only
**per-cell, post-founding** induction is detectable, through reduced between-metastasis variance.

**P11 Dynamics [known: Gunnarsson et al. 2023; analogue of Guan et al. 2025].** In the 2-state normalised
replicator–mutator ODE x′ = a(1 − x) − c x + Δb x(1 − x):
- (a) two snapshots (x0, x1) are fitted exactly both by pure selection (Δb = Δlogit / t) and by pure induction
  (a = −ln((1 − x1)/(1 − x0)) / t);
- (b) the vector field identifies (a, c, Δb), and for d states it identifies A = B + T up to A + cI;
- (c) at stationarity, a one-parameter family of (B, T) shares the same equilibrium.

**P12 Design law [new-form; elementary].** Suppose a metastatic state of frequency q is never observed among n_P
primary cells. Excluding pure selection under M(Λ) at level α′ requires the Clopper–Pearson upper bound for 0/n_P to be
below p*(q, Λ) = q / (Λ(1 − q) + q). That means n_P > ln(1/α′) / (−ln(1 − p*)) ≈ ln(1/α′)(Λ(1 − q) + q)/q. Sampling bias
ρ multiplies Λ (P9). Founder-linked lineage recording identifies τ directly with n ≈ ln(1/α)/τ metastatic cells, for any
Λ.
