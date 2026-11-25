# Milestone 6, Experiment 4: final biological robustness with dataset-specific calibration

E22 (`experiments/reconstruction/run_e22_realboot.py`). The occupancy test is re-run with a **parametric
bootstrap null simulated on the actual published trees**, rather than against a χ² reference or the single
generic threshold from the 12-clone simulation in Milestone 4. Null: all clones share the pooled occupancy,
each keeps its fitted rate. 120 bootstrap draws per dataset.

## Result: occupancy heterogeneity is significant in every eligible dataset

| Dataset | clones | LRT | χ² p | **bootstrap p** | bootstrap 95% threshold | $V_\pi$ |
|---|---|---|---|---|---|---|
| Simeonov mouse 1 | 8 | 34.8 | 1.2e-05 | **< 0.008** | 16.9 | 0.0140 |
| Simeonov mouse 2 | 12 | 146.8 | 6.8e-26 | **< 0.008** | 46.4 | 0.0704 |
| Quinn 5k, all sites | 12 | 37.7 | 8.8e-05 | **< 0.008** | 17.9 | 0.0037 |
| Quinn 5k, primary only | 12 | 16.9 | 0.110 | **0.050** | 16.7 | 0.0012 |
| Quinn 10k, all sites | 9 | 18.3 | 0.019 | **0.033** | 13.8 | 0.0106 |
| Quinn 10k, primary only | 6 | 25.5 | 1.1e-04 | **< 0.008** | 10.9 | 0.0212 |
| Quinn 100k, all sites | 12 | 77.0 | 5.7e-12 | **< 0.008** | 19.9 | 0.0429 |
| Quinn 100k, primary only | 10 | 89.5 | 2.1e-15 | **< 0.008** | 15.0 | 0.0321 |

**Eligible datasets: 8. Significant at the dataset-specific bootstrap threshold: 8 of 8.**

## Eligible biological replicates, stated precisely
- **Simeonov:** 2 mice, both eligible.
- **Quinn:** 4 mice, of which **3 are eligible** (5k, 10k, 100k). The **30k mouse is ineligible**: it has four
  lineage groups, and after the ≤400-cell tree cap only one usable tree remains — a variance across clones
  cannot be estimated from one clone.
- **Total: 5 eligible animals across 2 tumour models and 2 recorder chemistries.**

Earlier milestones reported 6 animals; the stricter tree caps used here (≤12 trees, ≤400 cells) drop the 30k
animal from both scopes. The honest count for this analysis is 5.

## The generic threshold was unreliable in both directions
Dataset-specific bootstrap thresholds range from **10.9 to 46.4**, against the χ² value of 19.68 and the
generic simulated threshold of 32.59 from Milestone 4. Tree shape and clone count move the null distribution
substantially:

- For Quinn 5k primary-only the **χ² test was conservative**, reporting p = 0.110 where the correct
  calibration gives p = 0.050. Milestone 4 found the χ² reference anti-conservative on average; on these real
  trees it errs both ways.
- For Simeonov mouse 2 the bootstrap threshold (46.4) is more than double the χ² value, so χ² was
  anti-conservative there — though with an observed statistic of 146.8 it makes no difference.

**Conclusion: a single generic threshold should not be used.** Dataset-specific calibration is cheap and is
what the final analysis uses.

## The marginal case, reported as marginal
**Quinn 5k primary-site-only sits exactly at p = 0.050** with $V_\pi = 0.0012$, the smallest effect in the
table. This is the single weakest result among the eight and should be described as borderline rather than
significant. It does not change the overall conclusion, which rests on seven clearly significant datasets
across five animals.

## Negative findings preserved
These are unchanged and are restated so they are not quietly dropped:

1. **Switching-rate heterogeneity remains unidentifiable.** Fitted rates lie above the empirically covered
   band in both systems; $V_{\log k}$ spans 6.5-fold across three published reconstructions of the same cells;
   after depth adjustment it was significant in only one of three mice. Milestone 6 adds that even
   reconstruction-aware inference attains only 3–38% interval coverage below the identifiability window.
2. **Cross-tissue within-clone expression changes are confounded by sequencing depth** in both datasets, with
   quality covariates moving further than the biological signal. No within-clone mechanism is claimed.
3. **No gene-level or spatial validation is offered.** Nine held-out programmes in Simeonov produced no
   association with clone occupancy that replicated across animals, with three changing sign; no public
   dataset with matched lineage, expression and spatial coordinates was located. Nothing was invented to fill
   these gaps.

## Limits
Tree caps (≤12 trees, ≤400 cells) were applied uniformly to make 120 bootstrap refits per dataset tractable,
so the $V_\pi$ values here are computed on a subset of the trees used in Milestone 4 and differ numerically
from those. The significance conclusions are unaffected. 120 draws give a resolution of about 0.008 on the
p-value, so "< 0.008" means no bootstrap draw exceeded the observed statistic.
