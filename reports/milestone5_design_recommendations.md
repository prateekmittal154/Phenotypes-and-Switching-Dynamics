# Milestone 5: experimental design, revised by the reconstruction results

Milestone 4 derived a recorder-capacity law and found it over-predicts real coverage by 0.70 decades.
Milestone 5 identifies *why*, and the answer changes the recommendation.

## 1. What the hierarchy experiment changes

| Milestone 4 belief | Milestone 5 finding |
|---|---|
| Real trees are star-like, costing 0.70 decades; cause unknown, possibly reconstruction artefact | **Sampling, not reconstruction.** The *true* sampled genealogy already has only 1.6–1.9% close-relative pairs; reconstructions give 1.2–3.2%. No algorithm can recover pairs the sample does not contain. |
| Topology error is the main inference risk | **Clock miscalibration is larger.** Treating one resolved branch as one time unit biases $\hat k$ by 1.83×; topology error adds only 1.04–1.10× more for NJ, UPGMA or an ensemble. |
| More recorder capacity widens the window | **Still true, and it also doubles topology recovery** (0.56–0.60 of true clades at 25 sites against 0.19–0.26 at 9 sites, where the recorder is 97% saturated). |

## 2. The revised design hierarchy, in priority order

**1. Fix the clock before anything else.** The largest single error source is converting resolved branches to
time. Two routes: calibrate the edit rate independently (a timed in-vitro arm with the same construct), or
record timestamps so that elapsed time is observed rather than inferred. A time-resolved recorder is valuable
here not because it widens the identifiability window — it does not — but because it removes the dominant
bias.

**2. Sample densely within clades, not broadly across the tumour.** This is the sharpest new recommendation
and it is counter-intuitive. The informative pairs for switching kinetics are close relatives, and sampling
*n* cells uniformly from a large expanding population yields almost none: the expected fraction of pairs
coalescing recently falls as the population grows. Deliberate dense sampling of a few microregions, or
micro-dissection of small clonal patches, concentrates the budget on the pairs that carry rate information.

**This conflicts with the design that maximises occupancy precision**, which wants many clones sampled broadly
(Milestone 4, §5). The two objectives want opposite sampling schemes, and an experiment must choose or split
its budget deliberately.

**3. Raise recorder capacity to avoid saturation.** At 9 sites the simulated recorder is 97% saturated within
ten generations, so recent divisions leave no mark and topology recovery halves. Capacity should be chosen so
that saturation stays well below 1 over the intended recording window: with per-site edit probability $\mu$ per
division and $g$ generations, the expected unedited fraction is $(1-\mu)^g$, so sites should satisfy
$S \gtrsim k_{\max} D$ (Milestone 4) **and** $\mu$ should be low enough that $(1-\mu)^{D} $ is not negligible.
These pull in opposite directions and are the real design tension in recorder chemistry.

**4. Deposit more than one reconstruction.** Costless, and it lets readers see the 6.5-fold spread in rate
estimates that a single tree hides. The Quinn deposit does this for one mouse and it is the only reason the
real-data ensemble analysis was possible.

## 3. Cost framing
We can state the direction but not a credible price. What the experiments support:

- **Capacity** is the cheapest lever for topology: doubling resolved sites roughly doubled clade recovery.
- **Dense within-clade sampling** costs nothing extra in sequencing — it is a change in *which* cells are
  taken, not how many — and is the only lever that addresses the missing close pairs.
- **A timed calibration arm** removes a 1.83× bias for the cost of one additional in-vitro experiment, with no
  change to the animal work.

We do **not** quantify the number of cells needed for a target precision. The honest obstacle is that the
naive pair-based sample-size calculation overstates information by 1–221× (Milestone 4 §3), and we have not
built a full-tree power calculation. Any figure we quoted would be unreliable, so none is given.

## 4. What would make switching rates measurable
Collecting the previous milestones' constraints, a design that could identify clone-specific switching needs
all four:

1. resolved levels spanning the plausible rate range, with a margin of ~0.7 decades over the idealised law;
2. **dense within-clade sampling** so that close-relative pairs exist at all;
3. an independent clock calibration, or timestamps;
4. enough clones (≥ 20) for a variance across clones, which pulls against (2) at fixed budget.

No published dataset we examined satisfies (1) and (2) simultaneously, which is a sufficient explanation for
the uniform failure to identify rates — without invoking any deficiency of the inference methods.
