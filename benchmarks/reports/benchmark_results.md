# Benchmark results

All numbers come from executed runs. Tables are in `results/tables/`, per-clone outputs in
`results/metrics/`, figures in `results/figures/`. The rules that govern which comparisons
are admissible were fixed in `protocols/estimand_matrix.md` before any run.

Replicate counts: Benchmark I, 200 replicates (8 mechanisms x 25) x 2 branch-length
conventions x 12 clones = 4,800 clone-level fits per method. Benchmark II, 75 replicates.
Benchmark IV, 20 replicates x 5 trees. Benchmark V, 116 real clone trees.

---

## 1. Occupancy is recovered by every method

On the true tree with true branch times, the correlation between each method's clone-level
occupancy estimate and the truth (Table B3, Figure B2A):

| Mechanism | Study D | PATH | PASTRI | SCOUT (theta) |
|---|---|---|---|---|
| occupancy varies (S3) | 0.910 | 0.910 | 0.889 | 0.872 |
| both vary (S4) | 0.850 | 0.830 | 0.806 | 0.418 |
| induction (S5) | 0.936 | 0.903 | 0.892 | 0.872 |
| bottleneck (S6) | 0.890 | 0.882 | 0.861 | 0.860 |
| unequal expansion (S7) | 0.897 | 0.887 | 0.838 | 0.817 |
| measurement error (S8) | 0.877 | 0.866 | 0.842 | 0.874 |
| stable states (S1) | 0.567 | 0.387 | 0.401 | **-0.070** |

The four methods agree closely, with Study D marginally ahead and the ordering stable
across mechanisms. Two exceptions matter. In **S1**, where states are nearly frozen
(k = 0.02), every method degrades and SCOUT's OU optimum carries no signal at all: with
almost no within-clone variation the OU model has nothing to estimate a reversion level
from. In **S4** SCOUT also drops to 0.418, because when the reversion rate varies across
clones the fitted optimum absorbs some of that variation.

The variance component actually reported in papers behaves the same way (Table B3b). For
S3, where the true `V_pi` is 0.0346, the estimates are Study D 0.0364, PATH 0.0395,
PASTRI 0.0392 — all within 15%. For S2, where the true `V_pi` is exactly 0, all three
return a small positive value (0.017–0.022) that is estimation noise, smallest for Study D.

**Conclusion for Q1: yes.** Clone-associated phenotypic occupancy is recovered under
realistic sampling by all four methods, and this does not depend on which method you use.

---

## 2. No method recovers the switching rate, and they fail in opposite directions

Median log bias in the inferred rate against the truth (Table B3, Figure B2B). Zero is
unbiased; +1 is a 2.7-fold overestimate.

| Mechanism | Study D | PATH | PASTRI | SCOUT (alpha) |
|---|---|---|---|---|
| rate varies (S2) | +0.454 | -0.522 | +2.501 | +0.134 |
| occupancy varies (S3) | +0.274 | -0.519 | +2.580 | +0.217 |
| both vary (S4) | +0.488 | -0.811 | +2.380 | +0.117 |
| induction (S5) | +1.100 | -0.472 | +2.788 | +0.221 |
| bottleneck (S6) | +0.275 | -0.441 | +2.576 | +0.198 |
| unequal expansion (S7) | +0.261 | -0.211 | +2.367 | +0.247 |
| measurement error (S8) | +0.782 | -0.521 | +2.699 | +0.362 |
| stable states (S1) | +2.902 | +2.629 | +4.852 | +0.985 |

Three independent implementations, three different directions: **PATH systematically
underestimates** by a factor of about 1.6–2.3, **PASTRI overestimates by roughly 11–16
fold**, and **Study D overestimates mildly**, by 1.3–1.6 fold in most mechanisms. SCOUT's
`alpha` is the closest to its own generating value, but it is a different estimand on a
different data type and is not evidence that the binary switching rate is recoverable.

The disagreement is larger than any bootstrap interval, so it is driven by modelling
assumptions, not sampling noise. That the three methods bracket the truth from opposite
sides is the strongest available evidence that **the data do not pin the rate down**: if
the rate were well identified, three correct implementations of related estimators would
not disagree by more than an order of magnitude on the same trees.

**S1 is the clearest case.** With k = 0.02 the true rate sits far below the window the
trees resolve, and every method is wrong by 2.6 to 4.9 log units — a 14 to 130-fold
overestimate. PATH and PASTRI also decline to return a value for 13% and 18% of clones
respectively, because clones fixed in one state admit no transition matrix.

### Spurious rate heterogeneity

In S3, S5, S6, S7 and S8 the true `V_log k` is exactly zero: every clone shares one rate.
What the methods report (Table B3b, median over replicates):

| Mechanism | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| S3 | 2.590 | 0.249 | 0.369 | 1.763 |
| S5 | 2.935 | 0.208 | 0.198 | 1.357 |
| S6 | 2.555 | 0.284 | 0.270 | 1.520 |
| S7 | 2.177 | 0.345 | 0.315 | 0.984 |
| S8 | 3.252 | 0.323 | 0.215 | 1.659 |

**Study D manufactures the most spurious rate heterogeneity of the three binary methods.**
This is a genuine negative result for our estimator and we report it as such. The reason is
visible in the complementary column: in S2, where the true `V_log k` is 1.056, Study D
reports 3.76 while PATH reports 0.30 and PASTRI 0.50. Study D over-disperses; PATH and
PASTRI compress toward a constant. Neither behaviour recovers the truth — PATH and PASTRI
look steadier only because they are shrinking, and a method that always reports "little
variation" will look stable whether or not variation exists.

The practical implication is the same for all three: **a reported spread of per-clone
switching rates is not evidence that clones differ in switching rate.**

**Conclusion for Q2:** when the rate lies outside the information the genealogy supports,
none of the four methods communicates that fact. Three of them return a point estimate with
no uncertainty at all. Study D returns an interval and flags boundary cases, which is the
difference in kind, not in accuracy.

---

## 3. Discarding branch times moves every rate estimate, and no occupancy estimate

Benchmark III re-analyses every Benchmark I replicate with all edges set to 1, changing
nothing else (Table B4, Figure B4).

| Method | median log(unit / true) for the rate | 95% CI | fraction bit-identical |
|---|---|---|---|
| Study D | +1.372 | [1.316, 1.465] | 0.040 |
| PATH | +1.557 | [1.546, 1.564] | 0.000 |
| SCOUT | +1.745 | [1.679, 1.830] | 0.000 |
| **PASTRI** | **0.000** | [0.000, 0.000] | **1.000** |

Occupancy, by contrast, barely moves: the median absolute change is 0.0026 for Study D and
exactly **0** for PATH and PASTRI, which return bit-identical occupancy for 96% and 100% of
clones.

Two things follow.

**First, the clock assumption costs a factor of about four on the rate** (e^1.37 to e^1.75),
for three methods built on quite different machinery, while leaving occupancy untouched.
This independently reproduces, in external software, the asymmetry the Study D manuscript
reports from its own estimator.

**Second, PASTRI's rate estimate is exactly invariant to branch times** — 505 of 505 clones
unchanged to the bit. This is not an error and is not scored as one. It follows from the
method's definition of depth, which counts nodes
(`lca_normalized_height`, computed from distance-to-root and distance-to-farthest-leaf in
`R/LCA_depths_calculation.R`) and never reads `edge.length`. The consequence is worth
stating plainly for anyone using the method: **PASTRI's per-cell-cycle rates are expressed
in units of topological depth, so they cannot be read as rates per unit biological time
unless divisions are synchronous.**

SCOUT shows the same issue from the other side. Its own shipped example is generated at
alpha = 3, and the installed engine recovers a median alpha of 0.309 on it — a ten-fold
scale error — because the shipped tree carries no branch lengths and SCOUT substitutes unit
edges by its own rule (`R/dropout_validation.R:470-472`).

**Conclusion for Q3:** the branch-time convention is the single largest controllable source
of error in the inferred rate, and it is invisible in the occupancy estimate.

---

## 4. Where methods decline to answer

Fraction of clones returning no finite estimate on the true tree (Table B4b):

| Mechanism | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| stable states (S1) | 0.000 | 0.130 | 0.180 | 0.000 |
| rate varies (S2) | 0.000 | 0.000 | 0.063 | 0.000 |
| occupancy varies (S3) | 0.000 | 0.000 | 0.067 | 0.000 |
| both vary (S4) | 0.000 | 0.000 | 0.090 | 0.000 |
| induction (S5) | 0.000 | 0.000 | 0.097 | 0.000 |
| bottleneck (S6) | 0.000 | 0.000 | 0.077 | 0.000 |
| unequal expansion (S7) | 0.000 | 0.000 | 0.030 | 0.000 |
| measurement error (S8) | 0.000 | 0.000 | 0.103 | 0.000 |

PATH's failures are confined to S1 and have a precise cause: a clone fixed in one state
makes `catMat()` return a single column, so no 2x2 transition matrix exists. PASTRI's 3–18%
are optimiser failures and cases where the implied switch probability reaches 1, which makes
the converted rate infinite. Study D and SCOUT always return something — but Study D's own
`k_at_bound` flag fires on 5.7–13.7% of clones, marking exactly the cases where its estimate
sits against the optimiser boundary and should not be read as a measurement.

The useful distinction is not how often a method fails but whether it tells you. Study D
and PATH signal their failures; PASTRI's infinities are visible only if you check; SCOUT
returns a finite alpha in every case including S1, where the true rate is far outside what
the tree can resolve.

---

## 5. Tree geometry predicts how wrong the rate will be

Benchmark II holds the truth fixed — **k = 0.3 for every clone in every geometry**, so the
true `V_log k` is exactly zero — and varies only the shape of the sampled genealogy. The
independent variable is the fraction of within-clone pairs whose `kL` falls in the
informative band [0.10, 1.0], which ranges 180-fold across the five geometries.

| Geometry | frac. in window | median MRCA depth |
|---|---|---|
| G_clustered | 0.149 | 5.0 |
| G_shallow | 0.055 | 6.9 |
| G_dispersed | 0.014 | 11.0 |
| G_starlike | 0.0035 | 10.6 |
| G_deep | 0.0008 | 12.5 |

Median log bias in the rate (Table B3c):

| Geometry | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| G_clustered (most informative) | +0.275 | **+0.149** | +1.984 | +0.196 |
| G_shallow | +0.217 | +0.054 | +2.101 | +0.084 |
| G_dispersed | +0.098 | -0.414 | +2.331 | +0.050 |
| G_starlike | +0.491 | -0.505 | +2.752 | +0.138 |
| G_deep (least informative) | +0.376 | -0.826 | +2.719 | +0.468 |

**The window predicts the bias, for every method** (Table B3d). The rank correlation
between the informative-pair fraction and the absolute log bias is negative throughout:
PASTRI -0.751, PATH -0.748, SCOUT -0.249, Study D -0.246 (75 replicates each). More
informative trees give smaller bias, as the theory says they should.

The practical consequence is the encouraging half of this benchmark: on the most
informative geometry, PATH's bias falls to +0.149 and Study D's to +0.275 — both close to
unbiased. **The methods are not broken; the genealogies usually are.**

### A caution about apparent stability

Spurious rate heterogeneity, where the truth is zero (Table B3c):

| Geometry | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| G_clustered | 1.516 | 0.493 | 0.701 | 0.260 |
| G_shallow | 1.435 | 0.368 | 0.575 | 0.447 |
| G_dispersed | 0.399 | 0.367 | 0.335 | 0.189 |
| G_starlike | 2.927 | 0.289 | 0.235 | 1.534 |
| G_deep | 2.710 | 0.208 | 0.213 | 1.646 |

The correlation with informativeness has **opposite signs** by method: Study D -0.534 and
SCOUT -0.609 (informative trees reduce spurious heterogeneity, as they should), but PATH
+0.386 and PASTRI +0.557 (informative trees *increase* it).

This resolves the apparent advantage PATH and PASTRI showed in Benchmark I. On
uninformative trees they collapse every clone onto nearly the same rate — low variance,
large bias. Only when the tree carries real information do their estimates begin to vary at
all. Their steadiness on poor data is shrinkage, not accuracy, and it would be read by a
user as "these clones share a switching rate" regardless of whether they do.

---

## 6. Reconstruction costs much less than the clock

Benchmark IV gives every method five trees over identical cells and states: the true
genealogy with true times, the same topology with unit edges, and three reconstructions
from simulated barcodes (greedy, NJ, UPGMA; mean recorder saturation 0.497).

**Occupancy is essentially untouched by reconstruction** (Table B4c):

| Tree | Study D | PATH | PASTRI | SCOUT |
|---|---|---|---|---|
| true | 0.910 | 0.912 | 0.896 | 0.901 |
| unit | 0.931 | 0.912 | 0.896 | 0.902 |
| greedy | 0.932 | 0.915 | 0.801 | 0.903 |
| NJ | 0.932 | 0.915 | 0.907 | 0.903 |
| UPGMA | 0.932 | 0.918 | 0.902 | 0.904 |

**The rate is not.** Movement away from the true-tree answer, median absolute log shift
(Table B4d):

| Tree | Study D | PATH | SCOUT | PASTRI |
|---|---|---|---|---|
| unit edges only (clock) | 1.395 | 1.563 | 1.431 | 0.000 |
| greedy (clock + topology) | 1.561 | 1.579 | 1.702 | 0.599 |
| NJ | 2.438 | 1.605 | 1.900 | 0.504 |
| UPGMA | 1.855 | 1.605 | 1.735 | 0.637 |

For PATH the unit-edge step alone accounts for 1.563 of the 1.579 total under greedy
reconstruction: **topology error adds 0.016**. For Study D it accounts for 1.395 of 1.561,
so topology adds about 0.17. The ordering the Study D manuscript reports from its own
estimator — clock first, topology a distant second — holds in external software.

The corresponding occupancy shifts are two orders of magnitude smaller: 0.0036–0.0044 for
Study D, 0.031–0.035 for PATH, 0.045–0.066 for PASTRI.

**Conclusion for Q4:** reconstruction uncertainty changes occupancy estimates negligibly
and rate estimates substantially, but almost all of the rate damage is done by the
branch-time convention the reconstruction forces, not by getting the topology wrong.

---

## 7. Held-out prediction: the one axis all four methods share

Estimand group D. Within each clone 30% of cells are withheld; each method is fitted on the
remaining 70% and its per-clone central tendency predicts the withheld cells. Skill is
`1 - MSE(method) / MSE(grand mean)`, so 0 means the method adds nothing over ignoring clone
identity. 50 replicates from scenarios S3 and S4 (Table B3e).

| Method | Target | Median skill | 95% CI | Median MAE |
|---|---|---|---|---|
| Study D | binary state | 0.157 | [0.122, 0.174] | 0.418 |
| PATH | binary state | 0.130 | [0.100, 0.163] | 0.415 |
| PASTRI | binary state | 0.116 | [0.081, 0.148] | 0.409 |
| SCOUT | continuous trait | 0.088 | [0.067, 0.120] | 1.047 |

**Every method carries real predictive signal**: all four skill scores are positive with
intervals clear of zero, so knowing which clone a cell belongs to genuinely helps predict
its state in cells that were never used for fitting.

**The ordering among them is suggestive, not established.** Study D's interval overlaps
PATH's substantially and PASTRI's considerably; at 50 replicates these three cannot be
separated. SCOUT is scored on a different target and its number is not comparable in level,
only in sign.

This is the fairest single comparison in the package, because the quantity scored is a
prediction about data rather than a parameter, and it says what the rest of the benchmark
says: on occupancy the methods are interchangeable.
