# Final decision

**Verdict: GO — integrate into the Cell Systems manuscript, with one claim corrected.**

The benchmark was commissioned to close the gap the manuscript itself names as its largest:
that PATH, SCOUT and PASTRI had not been run. They have now been installed from their
authors' repositories, reproduced on the authors' own examples, and run on identical data.
The gap is closed, the manuscript's central claims are strengthened rather than weakened,
and one claim about our own estimator must be corrected before submission.

---

## What was executed

| Benchmark | Scale | Status |
|---|---|---|
| Method reproduction | 3 methods on their own examples | done |
| I — ground truth, 8 mechanisms | 200 replicates x 2 conventions x 12 clones | done |
| II — tree geometry | 75 replicates, 5 geometries | done |
| III — branch-time calibration | every Benchmark I replicate, twice | done |
| IV — reconstruction | 20 replicates x 5 trees | done |
| V — real data | 116 clone trees, 6 animals | done |
| Group D — held-out prediction | 50 replicates, all four methods | done |

Roughly 4,800 clone-level fits per method on the ground-truth set alone. Every figure and
table is reproducible from `src/finalise.sh`; checksums in `checksums/BENCHMARK_SHA256.txt`.

---

## The three findings that justify integration

**1. Occupancy is a property of the data, not of our estimator.**
Four independent implementations resting on different mathematics recover clone occupancy
at correlation 0.81–0.94 with the truth, agree with each other on real data at a median
rank correlation of 0.854, and give the three binary methods `V_pi` within 15% of truth.
The manuscript's central observational claim no longer depends on the Study D likelihood.

**2. The non-identifiability of switching kinetics is demonstrated, not argued.**
On identical trees PATH underestimates the rate by 1.6–2.3x, PASTRI overestimates by
11–16x, and Study D overestimates by 1.3–1.6x. On real data their median rates differ by up
to an order of magnitude and their rank agreement falls to 0.308, with PATH against SCOUT
at 0.067. Three correct implementations of the same estimand do not disagree this much on
a quantity the data identify. This is a stronger argument than the coverage analysis the
manuscript currently relies on.

**3. The rate–time confounding result reproduces in external software.**
Discarding branch times shifts the rate by +1.37 (Study D), +1.56 (PATH) and +1.75 (SCOUT)
log units while leaving occupancy bit-identical for 96–100% of clones. PASTRI's rates are
**exactly** invariant — 505 of 505 clones unchanged — because its depth counts nodes.
SCOUT recovers alpha = 0.309 on its own shipped example generated at alpha = 3, for the
same reason. Proposition V4 is no longer a statement about our model.

Two further results tighten the manuscript's design message. Benchmark II shows the
informative-pair fraction predicts the bias for **all four** methods (Spearman −0.25 to
−0.75), so the identifiability window is not a Study D construct. Benchmark IV decomposes
reconstruction error and finds the clock accounts for 1.563 of PATH's 1.579 total shift,
with topology adding 0.016 — the manuscript's ordering, in someone else's code.

---

## The claim that must be corrected

**Study D is not a better rate estimator, and the manuscript must not imply it is.**

In mechanisms where every clone shares one rate, Study D reports a spurious `V_log k` of
2.2–3.3 against PATH's 0.21–0.35 and PASTRI's 0.20–0.37. It produces the **largest**
spurious rate heterogeneity of the three binary methods. On real data it reports an
interquartile spread of clone rates three to six times wider than PATH's.

The honest reading, supported by Benchmark II, is that neither behaviour is correct: Study
D over-disperses, PATH and PASTRI shrink toward a constant and under-report real variation
by two- to four-fold when it exists. Their apparent steadiness is shrinkage, and it reverses
on informative trees.

What Study D contributes is therefore a **diagnostic, not a precision gain**: it reports an
interval, states the window of rates the tree resolves, and raises a boundary flag on
5.7–13.7% of clones — including on the real animals, where its median rates of 17.75 and
19.49 sit against the bound. The three competitors return a confident point estimate in
every case, including where the true rate lies far outside what the trees can resolve.

`reports/integration_memo.md` §6 states exactly how to word this.

---

## Residual risk

| Risk | Assessment |
|---|---|
| A referee runs SCOUT's full pipeline and gets different behaviour | **Real.** We ran the OUwie engine, not `SCOUT()` with EM, lambda screening and hybrid selection. Scoped explicitly in every report; `runSCOUT.identifiability()` not run. |
| PASTRI at 48-cell clones is outside its design range | Mitigated but not eliminated. Clone size was raised from 24 to 48 after a pilot, the depth window is swept rather than fixed, and Benchmark II varies size and depth. PASTRI's published applications use far larger trees. |
| The ggtree substitution changed PASTRI's behaviour | Low. One data-frame conversion was replaced, from the same authors' `tidytree`; the estimator is unmodified and reproduces the shipped example. |
| Study D's authors wrote the benchmark | Mitigated by running every competitor from published code at documented defaults, reproducing each on its own example first, including two mechanisms outside Study D's likelihood family, and reporting the result that is unfavourable to Study D as a headline. |

Two bugs in **our own** code were found and fixed during the work, both logged: a root-
selection error that collapsed four Simeonov M1 clones to a single tip (B11), and a failure
to remove unary nodes that made SCOUT reject all 116 real trees (B12). Every real-data
result was regenerated after both fixes.

---

## Recommended next actions

1. Fold the benchmark into the manuscript per `reports/integration_memo.md`: one new main
   figure, a new Results subsection, a rewritten Discussion paragraph, a new Methods
   subsection and a new supplementary section.
2. Strike the "benchmark against the authors' software remains to be done" sentence from
   the Discussion and the matching entry from `submission/external_assets_needed.md`.
3. Rewrite Supplementary Section 5 item 2 to record what was run and what remains.
4. Correct any wording that implies Study D estimates rates more accurately.
5. Optional, if a reviewer asks: run `SCOUT()` end to end, PASTRI's permutation test, and
   the GESTALT-like recorder in Benchmark IV. All three are implemented and unexecuted; see
   `reports/deviations.md`, "Not done".

**No further compute is required before submission.**
