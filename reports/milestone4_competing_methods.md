# Critical assessment of the manuscript against competing methods

The brief asks for the manuscript to be assessed critically against PASTRI, SCOUT and the closest
lineage-based design work. This is written adversarially: the question is not whether our contribution is
defensible but where a referee from each of those groups would attack it.

## 1. Estimands first, because most comparisons would be unfair

| Method | Estimand | Overlaps ours? |
|---|---|---|
| **PATH** \[Schiffman 2024\] | One phylogenetic autocorrelation over the tree, converted to transition rates under a **shared-rate** model | **Partly.** Its autocorrelation mixes our two components; its rate is a population average, not per clone. |
| **PASTRI** \[Yang 2026\] | Transition rates resolved by **developmental stage**, using pairs at varying phylogenetic distance | **No.** Heterogeneity in time, not across clones. |
| **SCOUT** \[Stuart 2025\] | Per-gene Ornstein–Uhlenbeck drift and reversion on lineage trees | **No.** Continuous and per gene; the OU reversion rate is the analogue of $k$ but is never decomposed across clones. |
| **Multitype branching inference** \[Gunnarsson 2023\] | Absolute switching and net growth rates from **time series** of type abundances | **No.** Needs timepoints; a fixed-tree likelihood cannot see growth at all. |
| **This work** | $V_\pi=\mathrm{Var}_c(\pi_c)$, $V_{\log k}=\mathrm{Var}_c(\log k_c)$, plus an identifiability window and a capacity law | — |

A head-to-head accuracy table across these would be meaningless. The only fair comparison we can run is against
a PATH-style autocorrelation, which we reimplemented from its definition.

## 2. Where a referee would attack, and the honest answer

**"PASTRI already uses pairs at varying phylogenetic distance — your Figure 3 is their Figure 2."**
Partly fair, and it is disclosed. PASTRI's purpose is to let rates vary with developmental stage; ours is to
decompose a single snapshot's covariance into a decaying and a non-decaying part and to identify the asymptote
with between-clone occupancy variance. The shared ingredient is the idea of resolving by distance; the estimand
and the conclusion differ. **Our defence is the decomposition and the window, not the use of distance.** If a
referee insists the distance-resolved curve is PASTRI's, we lose nothing essential: Propositions 2–4 and the
capacity law stand.

**"SCOUT's OU reversion rate is your $k$, done better, and continuously."**
Largely fair for the single-clone case. SCOUT is a more flexible model of the same physics. What it does not do
is ask whether reversion rates differ *between clones* or state when that question is answerable. A referee
could reasonably ask why we did not simply add a clone random effect to SCOUT. **The honest answer is that we
should have tried**; it is listed as the main unaddressed alternative in §4.

**"Your method's headline real-data result is that nothing is identifiable. Why is that a method?"**
This is the sharpest attack and it has force. Our reply is that the contribution is the *diagnosis plus the
remedy*: the window tells you when an estimate is trustworthy, the empirical geometry check shows published
recorders fail it, and the capacity law says what to build. A method that reliably reports "not identifiable"
where competitors report a number is useful, but only if the field accepts a negative deliverable. That is a
real risk to acceptance, and it drives the venue recommendation.

**"Your calibration is on your own generative model."**
Partly answered: one scenario uses a three-state process collapsed to two, which is out of family, and the
occupancy test is anti-conservative there. But the misspecification tested is narrow. A referee could ask for
non-Markov state dynamics, rate variation along branches, or realistic reconstruction error, none of which we
did.

**"You take the published trees as truth."**
Fully fair. We do not propagate lineage-reconstruction uncertainty. Given that our central empirical claim is
about *tree geometry* — that real trees are star-like enough to cost 0.7 decades — a referee is entitled to ask
how much of that star-likeness is reconstruction artefact rather than biology. **This is the single most
damaging unanswered question in the manuscript**, because reconstruction methods are known to produce
poorly-resolved deep structure, and it is listed first among next experiments.

## 3. The fair comparison we did run
A PATH-style depth-weighted autocorrelation on the same pancreatic cells gives 0.524, 0.090, 0.098 (mouse 1,
three sites) and 0.109, 0.004, 0.007 (mouse 2). It agrees with our permutation test about where phylogenetic
structure exists. It cannot say whether that structure is between-clone occupancy or within-clone memory; in
mouse 1's primary our decomposition attributes $V_\pi = 0.053$ of it. The relationship is **complementary**:
the autocorrelation detects, the decomposition attributes.

| Property | PATH-style statistic | This work |
|---|---|---|
| Identifiability window stated | no | yes, and verified against real geometry |
| Separates occupancy from kinetics | no | yes, inside the window |
| Interval coverage reported | no | yes, including a regime where it fails (20%) |
| Dependence between pairs handled | n/a | quantified (1–221× inflation) |
| Design guidance | none | capacity law and Pareto designs |
| Computation | seconds | minutes to hours |
| Behaviour when unidentifiable | returns a number | flags it |

## 4. What we did not do, and should say so in the paper
1. **Did not run PATH, PASTRI or SCOUT implementations.** Only a PATH-style reimplementation. A proper
   benchmark on common data with matched estimands is the main methodological gap.
2. **Did not add a clone random effect to SCOUT**, which is the most obvious alternative route to our estimand.
3. **Did not propagate tree uncertainty**, which bears directly on our central geometric claim.
4. **Did not test non-Markov or branch-heterogeneous dynamics** beyond one collapsed three-state process.

## 5. Net assessment
The decomposition, the window, the dependence correction and the capacity law are defensible and, as far as our
searches found, unclaimed in this combination. The weakest points are the absence of a real benchmark and the
treatment of reconstruction error. Neither is fatal to the mathematics, but both are fixable and should be
fixed before submission rather than argued around.
