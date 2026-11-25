# Milestone 4: experimental design, with the capacity law derived and verified

## 1. The Fisher-information claim, restated precisely

Milestone 3 reported an optimum at $kL^\ast \approx 0.405$ without stating four things. They are fixed here
(`experiments/sampling_design/fisher_audit.py`).

| Specification | Milestone 3 | Milestone 4 |
|---|---|---|
| Parameter of interest | unstated | **log k**, k = a + b the total switching rate |
| Observation | "a pair" | one within-clone leaf pair: a single draw from a 4-cell multinomial |
| Meaning of L | unstated | **MRCA-to-each-leaf distance**; the path between the pair is 2L and ρ = e^{−2kL}. In total-pair-distance units d = 2L the optimum is k·d\* ≈ 0.86 |
| Occupancy π | held fixed | reported both ways; profiled out via the Schur complement |
| Pair independence | assumed | **false in a tree**, quantified in §3 |

**Corrected optimum.** Over k ∈ [0.05, 5] and π ∈ [0.1, 0.7]: with π known, kL\* = **0.43** (range 0.39–0.49);
with π estimated jointly, 0.43 (range 0.39–0.50), at a cost of at most **6%** of the information. The Milestone 3
value 0.405 came from sweeping only π ∈ {0.3, 0.5} and is superseded.

**The optimum is flat.** At least half the peak information is retained over **kL ∈ [0.10, 1.0]**, a ten-fold
band. This is a numerical result, not a closed form.

## 2. The capacity law, derived

Statements are labelled **[identity]** (follows algebraically from the stated model), **[approximation]**
(numerical, with the accuracy given), or **[empirical]** (measured from data).

**Model.** A clone grows for $D$ divisions; the recorder has $S$ usable sites, so the reconstructed tree
resolves $D_{\text{eff}} = \min(S, D)$ coalescent levels, each spanning $b = D / D_{\text{eff}}$ divisions.

**Step 1 [identity].** In a *balanced* resolved tree, a pair whose MRCA sits at resolved level $h$ has
$L = (D_{\text{eff}} - h)\,b$ for $h = 0, \dots, D_{\text{eff}}-1$. Hence the available distances are
$\{b, 2b, \dots, D_{\text{eff}} b\}$, so
$$L_{\min} = b = D/D_{\text{eff}}, \qquad L_{\max} = D .$$

**Step 2 [approximation].** A pair at distance $L$ retains at least half the peak Fisher information about
$\log k$ when $kL \in [0.10, 1.0]$ (§1; verified over two decades of $k$ and $\pi \in [0.1, 0.7]$, band
endpoints stable to about $\pm 0.01$).

**Step 3 [identity, given Steps 1–2].** A rate is covered if some available distance places it in the band:
$$k \in \left[\frac{0.10}{L_{\max}},\ \frac{1.0}{L_{\min}}\right] = \left[\frac{0.10}{D},\ \frac{D_{\text{eff}}}{D}\right].$$

**Step 4 [identity].** The width in decades is
$$W = \log_{10}\!\frac{D_{\text{eff}}/D}{0.10/D} = \log_{10}(10\,D_{\text{eff}}) = 1 + \log_{10} D_{\text{eff}} .$$
**$D$ cancels.** Recording for longer shifts the window downward but does not widen it; only resolving more
levels widens it. This is the content of the law, and it is an algebraic consequence of Steps 1–3, not a new
theorem.

**What is new here** is not the algebra but Step 2 (the information band, computed in this work) and the
consequence that the *ceiling* $D_{\text{eff}}/D$ is set by the finest resolvable branch rather than by tree
depth or cell number.

## 3. Full-tree information versus the naive pair sum
Pairs in a tree share ancestors, so summing per-pair information double-counts. Exact Fisher information of the
full-tree likelihood was computed by enumerating all $2^n$ leaf configurations for 8-leaf trees and by Monte
Carlo (validated against the exact values) for deeper trees. Results:
`reports/milestone4_full_tree_information.md`.

## 4. Verification against the real tree geometry — the idealisation is optimistic

The law assumes balanced trees. Real reconstructions are not balanced, so the covered band was recomputed from
the **empirical** distribution of within-clone pair distances in the published trees
(`experiments/sampling_design/capacity_verification.py`), taking $L_{\min}$ and $L_{\max}$ as the 10th and 90th
percentiles of that distribution.

| System | dataset | clones | median max depth | empirical $L_{\min}$ | empirical covered k | empirical decades | law predicts |
|---|---|---|---|---|---|---|---|
| macsGESTALT | mouse 1 | 12 | 4.0 | 1.21 | 0.035 – 0.83 | 1.37 | 1.60 |
| macsGESTALT | mouse 2 | 24 | 3.5 | 0.50 | 0.052 – 2.00 | 1.35 | 1.54 |
| Cassiopeia | m10k | 11 | 8.0 | 4.26 | 0.015 – 0.235 | 1.19 | 1.90 |
| Cassiopeia | m100k | 19 | 9.0 | 3.31 | 0.013 – 0.302 | 1.28 | 1.95 |
| Cassiopeia | m30k | 4 | 13.0 | 9.19 | 0.008 – 0.110 | 1.11 | 2.11 |
| Cassiopeia | m5k | 82 | 20.0 | 10.81 | 0.005 – 0.092 | 1.24 | 2.30 |

**Median over-prediction: 0.70 decades** (up to 1.06 for the deepest Cassiopeia dataset).

**Why.** Real trees are star-like: most within-clone pairs coalesce near the root, so close-relative pairs are
scarce. The empirical $L_{\min}$ is 3–11 resolved branches in Cassiopeia rather than the single branch the
balanced model assumes. Deep trees therefore do **not** deliver the fine-grained pairs needed for a high
ceiling.

**Two consequences, one of which corrects an earlier claim.**
1. **Realised coverage is remarkably uniform across recorders: 1.1–1.4 decades.** Cassiopeia's extra depth buys
   essentially nothing in width, because the topology cancels it.
2. **An earlier Milestone 4 statement that "Cassiopeia raises the covered ceiling from 0.40 to 0.86 per
   division" was based on the idealised law and is wrong on the empirical geometry.** Empirically Cassiopeia's
   ceilings are *lower* (0.09–0.30 per resolved branch) than macsGESTALT's (0.83–2.00). Note the units are not
   comparable across systems — a resolved branch means different amounts of real time in each — so this is a
   within-system statement about whether fitted rates fall inside the covered band, not a ranking of recorders.

Only **13.8%** of real clones cover k = 1 per resolved branch, and **11.8%** cover k ≥ 2.

## 5. Budgeted designs
At a fixed 2,048-cell budget, scoring information about the occupancy *variance* (which scales with clone count,
unlike information about the grand occupancy, which the budget fixes — deviation D21): every Pareto-optimal
design uses 64 clones × 32 cells, so breadth of clone sampling wins for both objectives at this budget, and the
frontier then trades recorder capacity (S = 16 maximises worst-case rate information; S = 2 maximises
occupancy-variance information because coarser trees leave cells less correlated). 5 of 100 designs are
Pareto-optimal.

**Robust (minimax) design** over k ∈ [0.02, 5] selects the highest-capacity design, because the worst case is
always the fastest rate and that is governed by the finest branch.

**A design that does not work.** Splitting the budget across two edit rates was tested and **loses** (0.11
against 0.23 worst-case): the slower-editing arm has strictly coarser resolution. Mixing rates is not a
substitute for capacity.

## 6. Recommendations
1. **Set capacity from the fastest rate you need:** $S \gtrsim k_{\max} D$ — and then, because real trees are
   star-like, add a margin of roughly **0.7 decades** over the idealised law (§4).
2. **Set duration from the slowest:** $D \gtrsim 0.10 / k_{\min}$.
3. **Spend the cell budget on clones**, not on depth per clone.
4. **Report the empirical covered band** $[0.10/L_{\max}, 1.0/L_{\min}]$ computed from the realised pair-distance
   distribution — not from the idealised formula — alongside any switching-rate estimate. On the two datasets
   analysed here, the fitted rates fall outside that band in both systems.

## 7. Limits
The band in Step 2 is numerical. The law assumes a constant edit rate per division and the two-state stationary
model. §4 removes the balanced-tree assumption but still treats resolved branches as a clock. Dropout and
saturation reduce $D_{\text{eff}}$ below $S$, so the law must be applied with realised resolved levels.
