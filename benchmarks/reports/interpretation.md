# Interpretation: what Study D adds, and what it does not

The question this benchmark was built to answer is not whether our estimator wins. It is
which biological properties lineage-resolved data actually support, how robustly existing
methods recover them, and what our framework contributes to reading them.

---

## 1. Conclusions shared across methods

**Clone-associated phenotypic occupancy is real, recoverable, and method-independent.**
Across 200 ground-truth replicates and four independent implementations resting on
different mathematics — a two-state CTMC likelihood, a Moran's-I style phylogenetic
autocorrelation, a depth-stratified adjacency optimisation, and an Ornstein-Uhlenbeck
model on a continuous trait — the correlation between the estimated and true clone
occupancy is 0.81–0.94 in six of eight mechanisms, and the three binary methods recover
the across-clone variance `V_pi` to within 15% of truth.

This matters for the manuscript because it converts a claim about our estimator into a
claim about the data. Occupancy heterogeneity is not an artefact of the Study D likelihood;
three programs written by other groups, with different assumptions, see the same thing.

**Discarding branch times costs a factor of about four on the rate and almost nothing on
occupancy.** Replacing true branch lengths with unit edges shifts the inferred rate by
+1.37 (Study D), +1.56 (PATH) and +1.75 (SCOUT) log units, while the occupancy estimate
moves by a median of 0.0026 or less and is bit-identical for 96–100% of clones under PATH
and PASTRI. The asymmetry the Study D manuscript reports from its own estimator is
reproduced in external software.

---

## 2. Quantities only particular methods estimate

- **Study D alone** reports an identifiability window and profile intervals, and flags
  clones whose rate estimate sits at an optimiser boundary (5.7–13.7% of clones).
- **PATH alone** provides an analytical significance test for phylogenetic signal itself,
  through the z-score on `xcor`.
- **PASTRI alone** allows the rate to change with depth, which is the right structure for a
  staged developmental process and which no other method here offers.
- **SCOUT alone** models a continuous phenotype natively and selects between evolutionary
  models, rather than assuming one.

These are complements, not competitors. A study that wants a significance test on
heritability should use PATH's; one that suspects stage-dependent dynamics should use
PASTRI's structure; one with a graded phenotype should not binarise it for us.

---

## 3. Where methods disagree, and why

**Switching rates.** On identical trees with identical states, PATH underestimates the true
rate by 1.6–2.3x, PASTRI overestimates by 11–16x, and Study D overestimates by 1.3–1.6x.
Three correct implementations, run on the same data, bracket the truth from both sides and
differ from each other by more than an order of magnitude.

The disagreement is traceable to assumptions rather than to noise:

- **PASTRI's depth is topological.** `lca_normalized_height` counts nodes and never reads
  `edge.length`, so its per-cell-cycle rate is expressed in units of topological depth. Its
  estimate is bit-identical on true and unit branch lengths in 505 of 505 clones. Its bias
  also grows as more distal pairs enter the window (log bias 1.42 at `d <= 0.2` rising to
  2.83 at `d <= 1.0`), which is exactly the behaviour the authors themselves describe.
- **PATH's weight matrix uses one-node-depth pairs**, so it draws on the closest relatives
  only. Its rate estimate is nearly flat in the true rate across the range tested
  (Figure B3C): it is reporting something closer to a fixed property of tree shape than a
  measurement of kinetics.
- **Study D fits a per-clone likelihood with few informative pairs**, which produces an
  over-dispersed spread of clone-level rates rather than a compressed one.

**This is the central result of the benchmark.** If switching rates were identified by
these data, three independent estimators of the same quantity would not disagree this
much. The disagreement is the evidence.

---

## 4. Does Study D improve uncertainty calibration, or only diagnose identifiability?

Honestly: **mainly the latter, and the diagnosis is the more useful contribution.**

Study D does not estimate the rate better than the alternatives in the sense of variance.
It produces the **largest** spurious rate heterogeneity of the three binary methods: in
mechanisms where every clone shares one rate, it reports a `V_log k` of 2.2–3.3 against
PATH's 0.21–0.35 and PASTRI's 0.20–0.37. We report this plainly because it is a real
weakness.

But the comparison that looks favourable to PATH and PASTRI on that metric is misleading in
a specific way. In the mechanism where rates genuinely do vary (`V_log k` = 1.056), PATH
reports 0.30 and PASTRI 0.50 — they under-report real variation by two to four fold, while
Study D over-reports it. Neither recovers the truth. PATH and PASTRI look steadier because
they shrink toward a constant, and a method that always says "little variation" will look
stable whether or not variation exists.

What Study D adds is not a better point estimate. It is that:

1. it reports an interval rather than a bare number;
2. it states the window of rates the tree can resolve, so a user can see when their
   estimate lies outside it;
3. it flags its own boundary cases rather than returning a confident-looking value.

The three competitors return a point estimate with no uncertainty in every case, including
the stable-state mechanism where the true rate lies far outside what the trees resolve and
every method is wrong by 14–130 fold. None of them signals that.

**So the contribution is a diagnostic, not a precision gain.** That is a weaker claim than
"our method is more accurate", and it is the claim the evidence supports.

---

## 5. Do the design recommendations generalise beyond our implementation?

Yes, and this is the result with the clearest practical reach.

The Study D manuscript argues that the supply of informative close-relative pairs, not the
reconstruction algorithm, is what limits kinetic inference, and that the branch-time
convention is the dominant controllable error. Both claims are properties of the data, so
they should hold for any method. They do:

- the unit-edge penalty appears in three independent methods at a similar magnitude;
- PASTRI's own bias grows monotonically as the depth window admits more distal pairs,
  reproducing the information-window argument in software written by another group;
- SCOUT recovers alpha = 0.309 on its own shipped example generated at alpha = 3, because
  that tree has no branch lengths and the package substitutes unit edges.

A reader therefore does not have to adopt our estimator to act on our design
recommendations. **Harvest at more than one timepoint, calibrate the clock, and report the
rate window rather than a point estimate** are instructions that apply whichever of these
four tools is used.

---

## 6. What a biological researcher should do differently

1. **Report clone-associated occupancy, and trust it.** It is recoverable, robust to
   method choice, and insensitive to the branch-time convention.
2. **Do not report a per-clone switching rate from a single-timepoint lineage experiment
   without a resolvable window.** Three independent tools disagree by more than an order of
   magnitude on the same trees.
3. **Do not read a spread of per-clone rates as biological heterogeneity.** Every method
   here manufactures one when the truth is a single shared rate.
4. **Check what your tool does with branch lengths before interpreting its rate.** PASTRI's
   rates are in units of topological depth; SCOUT substitutes unit edges silently when a
   tree has none; PATH's Q is scaled by the mean patristic distance it was handed.
5. **If the phenotype is graded, use a method that models it as graded.** Binarising to fit
   a two-state model discards information that SCOUT can use.
