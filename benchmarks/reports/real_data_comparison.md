# Real-data comparison: Simeonov and Quinn

All four methods were run on the same 116 clone trees, exported through the loaders the
canonical Study D analysis uses, so the trees, the clone assignments and the state
definitions are identical to those in the manuscript. Animals are kept separate throughout
and cells are never pooled across animals.

Continuous trait: the published EMT programme score (Simeonov), the cell-cycle score
(Quinn). Binary state: the same thresholded definitions the canonical analysis uses.
Branch lengths: unit edges, because the published trees carry no usable branch times —
the same convention SCOUT applies internally and PASTRI implies by construction.

There is no ground truth here. The question is not accuracy but agreement.

Source: `results/tables/TableB5_real_data.csv`, `TableB5b_cross_method_agreement.csv`,
`results/metrics/real_clone_level_all.csv`. Figure B5.

---

## 1. Coverage

| Animal | Clone trees | Cells | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|---|---|
| Simeonov M1 | 8 | 615 | 8 | 7 | 6 | 8 |
| Simeonov M2 | 21 | 1,706 | 21 | 21 | 15 | 21 |
| Quinn 5k | 64 | 5,999 | 64 | 64 | 64 | 64 |
| Quinn 10k | 9 | 801 | 9 | 9 | 9 | 9 |
| Quinn 30k | 1 | 132 | 1 | 1 | 1 | 1 |
| Quinn 100k | 13 | 1,141 | 13 | 13 | 13 | 13 |

Counts are clones returning a finite occupancy estimate. PASTRI declines most often — 6 of
21 Simeonov M2 clones and 2 of 8 M1 clones — and for the switching rate it returns nothing
for 12 of 21 M2 clones. Study D and SCOUT return an estimate for every clone.

The Quinn 30k animal contributes a single eligible clone, so no across-clone variance can
be formed from it. It is carried in the tables and excluded from every across-clone
statistic, exactly as in the manuscript.

## 2. Occupancy: the methods agree

Across-clone occupancy variance, scaled within each method so the four are comparable
(Figure B5A):

| Animal | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| Simeonov M1 | 0.111 | 0.117 | 0.089 | 0.097 |
| Simeonov M2 | 0.072 | 0.074 | 0.068 | 0.045 |
| Quinn 5k | 0.033 | 0.035 | 0.032 | 0.032 |
| Quinn 10k | 0.088 | 0.077 | 0.085 | 0.092 |
| Quinn 100k | 0.089 | 0.079 | 0.100 | 0.069 |

**Every method finds clone-associated occupancy structure in every eligible animal**, at
closely similar magnitude. The agreement is not merely in the summary: within each animal,
the rank correlation between any two methods' per-clone occupancy estimates has a median
of **0.854**, and Study D against PATH reaches **0.957** (Figure B5C).

This is the real-data counterpart of the simulation result. The manuscript's central
observational claim — that clone-associated occupancy is heterogeneous in these animals —
does not depend on the Study D estimator. Three other programs see the same clones as high
and low, in the same order.

## 3. Switching rates: the methods disagree

Median rate per animal and the spread across clones:

| Animal | Study D median | PATH median | PASTRI median | Study D IQR(log) | PATH IQR(log) | PASTRI IQR(log) |
|---|---|---|---|---|---|---|
| Simeonov M1 | 2.12 | 1.23 | 3.05 | 3.18 | 0.54 | 1.19 |
| Simeonov M2 | 19.49 | 1.87 | 7.62 | 2.72 | 0.47 | 1.34 |
| Quinn 5k | 1.39 | 1.06 | 5.34 | 3.01 | 0.63 | 0.80 |
| Quinn 10k | 1.45 | 0.84 | 4.75 | 3.13 | 0.56 | 0.84 |
| Quinn 100k | 17.75 | 1.19 | 5.58 | 2.63 | 0.42 | 0.51 |

Two things are visible and both match the simulations.

**The levels differ by up to an order of magnitude.** On Simeonov M2 the same cells give a
median rate of 1.87 under PATH, 7.62 under PASTRI and 19.49 under Study D.

**The spreads differ by three- to six-fold.** Study D reports an interquartile range of
2.6–3.2 log units across clones; PATH reports 0.42–0.63. Benchmark I showed which of these
is closer to the truth when the truth is known: neither. Study D over-disperses, PATH
compresses.

**Study D's own boundary flag fires here.** Median rates of 17.75 and 19.49 on the two
largest-spread animals sit against the estimator's upper bound, which is precisely the
signal the manuscript says should be read as "this tree does not resolve the rate" rather
than as a measurement of fast switching.

Cross-method rank agreement on the rate has a median of **0.308**, against 0.854 for
occupancy. The weakest pair, PATH against SCOUT, agrees at **0.067** — essentially not at
all.

## 4. Where the disagreements come from

Every method saw identical cells, identical states and identical trees, so the differences
are attributable to assumptions, and the simulation benchmarks identify which:

- **Branch-time convention.** None of these trees carries usable branch times. Benchmark
  III showed that this alone shifts rate estimates by about four-fold and leaves occupancy
  bit-identical. PASTRI's rates are expressed in topological depth by construction.
- **Which pairs each method uses.** PATH weights one-node-depth pairs, PASTRI
  depth-stratified adjacencies, Study D the full per-clone likelihood. Benchmark II showed
  the available informative pairs govern the bias for all four.
- **Shrinkage versus dispersion.** Benchmark II showed PATH and PASTRI compress clone
  rates toward a constant on uninformative trees. These real trees are uninformative in
  exactly that sense, which is why their spreads are narrow here.
- **Estimand.** SCOUT's alpha is an OU reversion rate on a continuous trait, not a binary
  switching rate, which is why its rate column is reported but never placed on the same
  axis as the others.

## 5. What this does and does not support

**Supported.** Clone-associated occupancy heterogeneity in these animals is a robust
observation, reproduced by four independent implementations with high rank agreement. The
manuscript's occupancy claims stand without reference to our estimator.

**Not supported.** Any per-clone switching-rate claim from these datasets. Four methods on
identical input disagree about both the level and the spread, and the one method that
reports a boundary flag raises it.

**Not evidence of replication.** Agreement between methods applied to the same cells is
agreement about one dataset, not independent biological replication. The five-animal
replication claim in the manuscript rests on the animals, not on the methods.
