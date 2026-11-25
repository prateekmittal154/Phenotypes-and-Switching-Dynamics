# Evaluation rules

Fixed before results were inspected.

## What counts as a fair comparison

1. A method is scored only on a quantity it natively estimates
   (`protocols/estimand_matrix.md`).
2. Methods are compared on identical inputs: same cells, same states, same tree, same
   branch-length convention.
3. Where a conversion is needed to put two methods on one scale, the conversion is stated
   in `benchmark_design.md` and applied to both.
4. A method that declines to return a value is recorded as "no estimate", with the count
   and the reason. It is never scored as if it had returned zero or infinity.
5. A method that reports no uncertainty is marked "not applicable" for coverage. Absence of
   an interval is not evidence of a bad interval.

## What is reported when a method fails

The failure, its frequency with a Wilson interval, and the message. Failures are a property
of the method on that data, not noise to be dropped. `TableB4b_no_estimate_rates.csv`
carries these counts for every scenario and method.

## Boundary behaviour

An estimate at an optimiser bound is not a measurement. We record:
- Study D: the `k_at_bound` flag the estimator already returns.
- PASTRI: cases where the implied switch probability reaches 1, which makes the converted
  rate infinite.
- PATH: cases where `logm` fails and the optimiser fallback is used.

These are reported alongside the point estimates, never silently filtered out before taking
a median.

## Direction of the claims we allow ourselves

- "Method X recovers quantity Q" requires correlation with the truth **and** a bias
  interval that is informative at the replicate count used.
- "Method X and Y agree" requires a rank correlation within animals, never pooled across
  animals.
- "Method X is better than Y" is only written where the two estimate the same quantity on
  the same input, and where the difference is larger than the bootstrap interval.
- Agreement between methods on the same dataset is **not** independent replication, and is
  never described as such.

## Deliberate asymmetries we accept

- Study D is the method under test and its authors wrote this benchmark. To limit the
  advantage that gives: every competitor was run from its authors' own published code; each
  was first reproduced on its authors' own example; settings follow each package's
  documentation and published usage rather than anything tuned here; and the scenario set
  includes two mechanisms (S2, S5) that Study D's own time-homogeneous model does not
  generate, so it is tested outside its own likelihood family.
- PASTRI was designed for trees with hundreds to thousands of tips. Running it per clone on
  48 cells is a harder setting than its published applications. Benchmark II varies tree
  size and depth so that this is visible rather than hidden, and the depth sweep is reported
  in full.
- SCOUT models a continuous trait. It is never scored on the binary state.
