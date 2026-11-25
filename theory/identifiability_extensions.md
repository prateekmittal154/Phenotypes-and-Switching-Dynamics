# Identifiability extensions: what each measurement buys (Proposition T5)

Each row states what is identified *in addition* to the row above, under the model of `branching_model.md`.
"Identified" means the parameter is a function of the observable law; finite-sample power is a separate question
and is reported in `../reports/milestone2_simulation.md`.

| Design | Newly identified | Not identified | Evidence |
|---|---|---|---|
| **Two state distributions only** (Milestone 1) | τ bounded under an assumed selection limit Λ | everything else; pure induction always fits | M1 propositions P1–P4 |
| **+ clone labels** | Var_c(F_c), the variance of *realised* clone fractions | the split of Var_c(F_c) into propensity, state memory and founder sampling | T4, T5.1, counterexample C6; simulation: naive statistic reads 0.040 when true B = 0 |
| **+ a resolved lineage tree** | A and B separately, **iff** 2k(L_max − L_min) is bounded away from 0 | the split when the decay is unresolved (k ≈ 0) | T2, T3, counterexample C7; real data: condition failed in 7 of 8 compartments |
| **+ compartment labels** | induction as a compartment main effect | per-clone induction, unless the same clones appear in two compartments | T2.1; real data: M2 adds ~0 AIC over M0 |
| **+ the same clones in two compartments** | the clone-by-environment interaction | its separation from compartment-specific technical effects, without a technical control | real data: within-clone shift confounded (figure 11) |
| **+ a second time point** | absolute rates a_c, b_c and growth rates | — | T5.4: a fixed-tree likelihood conditions on the tree, so growth never enters |
| **+ founder-resolving lineage marks** | B without the π(1 − π)/N_eff inflation | — | T4 |

## Practical design implications
1. **Estimate B in the primary tumour, not in metastases.** Metastases are oligoclonal, and simulation shows B̂
   attenuated 3–6 fold there.
2. **Maximise the spread of MRCA depths, not the number of cells.** The binding constraint is the identification
   condition in T3, which depends on the range of lineage distances the recorder resolves, not on cell count. In
   the real data, mouse 2's primary was informative with 1,619 cells while mouse 1's peritoneal metastasis was not
   with 3,691.
3. **Always report a technical control for compartment comparisons.** Held-out programmes and QC covariates shifted
   as much as the programme of interest.
4. **Treat a clone-model selection win as a bound, not a mechanism.** Counterexample C8.
