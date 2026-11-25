# Milestone 6, Experiment 3: sampling-design optimisation at fixed budget

E21 (`experiments/reconstruction/run_e21_sampling.py`). Five strategies, each spending the same **96 cells**,
across 2 growth models × 2 recorders × 2 true rates × 5 replicates = 40 observations per strategy. Inference
uses the best variant from E20 (ensemble over three reconstructions with ultrametric branch times).

**Scale warning.** The run took 3.2 hours for 5 replicates; a 30-replicate version would take ~19 hours and
was not run. The close-relative-pair measurements are precise, but the coverage comparisons are underpowered
and their Wilson intervals overlap. This is stated rather than glossed.

## Strategies
| Strategy | Practical? | Description |
|---|---|---|
| random | yes | uniform from the final live population — what is usually done |
| balanced | yes | equal numbers per founder clone |
| twotime | yes | half the budget at an earlier generation, half at the end |
| recorder_bc | yes | **prospective**: pre-screen 4× the budget by barcode similarity, then preferentially take cells that look closely related. Uses only recorded information, so an experimenter could actually do it |
| oracle | **no** | picks sibling-rich sets using the true genealogy — an unattainable upper benchmark |

## Results

| Strategy | close-relative pairs | median log bias | interval coverage [95% Wilson] | log interval width |
|---|---|---|---|---|
| random | 0.017 | +0.554 | 0.725 [0.572, 0.839] | 3.87 |
| balanced | 0.018 | +0.812 | **0.575 [0.422, 0.715]** | 4.13 |
| recorder_bc | 0.021 | +0.554 | 0.750 [0.598, 0.858] | 3.87 |
| **twotime** | **0.052** | **+0.468** | **0.825 [0.680, 0.913]** | 3.87 |
| oracle (unattainable) | 0.115 | +0.210 | 0.850 [0.709, 0.929] | 3.61 |

### What is solid
**Two-timepoint sampling triples the supply of informative pairs** (0.052 against 0.017), consistently across
both growth models and both recorders (0.051–0.064 against 0.014–0.023). This is the one effect large enough
to be unambiguous at this sample size.

**The oracle is only about twice as good again** (0.115). So a practical, zero-cost change recovers roughly
half of what perfect genealogical foresight would buy.

**The mechanism is sampling fraction, not timing.** Sampling 48 cells from the ~150-cell population present at
generation 7 captures a far larger share of that population than 96 cells from the ~580 present at generation
10. Close relatives are abundant when you sample a large fraction of a small population. This is the same
mechanism that makes published trees star-like (Milestone 5), seen from the other side, and it generalises:
**dense sampling of a sub-population works as well as early sampling.**

**The clever strategy disappoints.** Recorder-informed prospective selection raises the close-pair fraction
only from 0.017 to 0.021. Barcode similarity is too noisy a proxy for recent shared ancestry to beat simply
harvesting at a time when the population is small.

**Balanced clone sampling is the worst for rate inference** (coverage 0.575, widest intervals). Spreading the
budget evenly across clones deliberately avoids close relatives. This is the occupancy–kinetics tension made
concrete: balanced sampling is what you want for estimating occupancy variance, and it is actively harmful for
estimating rates.

### What is not solid
The coverage ordering (twotime 0.825, oracle 0.850, random 0.725) is suggestive but the intervals overlap at
n = 40. The bias ordering points the same way (+0.47, +0.21, +0.55). A properly powered repeat is the obvious
next computation, and it is cheap to specify: 30 replicates, which needs roughly a day of wall time.

## Recommendation
At fixed cost, in priority order:
1. **Harvest at two or more timepoints**, or equivalently sample densely within small sub-populations. Free,
   and it triples the informative-pair supply.
2. **Do not balance across clones if kinetics is the target** — and do balance if occupancy is. The two
   objectives want opposite designs and an experiment should choose explicitly.
3. **Do not expect barcode-guided prospective selection to help much** on its own.
4. Combine with an independent clock calibration (Experiment 2), without which the rate remains biased
   whatever the sampling.
