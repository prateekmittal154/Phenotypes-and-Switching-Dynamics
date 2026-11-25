# Milestone 5: tree reconstruction and its effect on inference

Two complementary experiments: synthetic, where the true genealogy is known (E17), and real, where three
*published* reconstructions of the same alleles are available (E18).

## 1. Real data: three published reconstructions of the same cells

The Quinn deposit contains three reconstructions per lineage group — the authors' primary
`hybrid_priors.alleleThresh`, plus `nj` and `nj_priors`. These are defensible alternatives produced by the
original authors from identical allele data, so they give a reconstruction ensemble without any modelling
choice of ours.

**Availability, stated precisely.** All three variants are deposited **only for the 5k mouse** (82 hybrid, 81
nj, 81 nj_priors trees). The 30k, 10k and 100k mice have the hybrid tree only (4, 11 and 19 trees). The
ensemble analysis therefore covers one animal — a limit of the deposit, not of the method. Twenty-five clones
are present in all three reconstructions and were analysed.

### Geometry differs enormously between equally defensible trees

| Reconstruction | median height | mean MRCA depth / height | fraction of close pairs | $L_{\min}$ | $L_{\max}$ |
|---|---|---|---|---|---|
| hybrid (authors' primary) | **15** | 0.196 | 0.001 | 10 | 13 |
| nj | **28** | 0.220 | 0.005 | 15 | 25 |
| nj_priors | **35** | 0.244 | 0.010 | 21 | 32 |

The same cells and the same alleles give trees whose heights differ **more than twofold**. Any quantity measured
in units of "resolved branches" — including the covered rate band of the Milestone 4 capacity law — inherits
that ambiguity directly.

**But all three agree that close-relative pairs are almost absent**: 0.1–1.0% of pairs coalesce in the top 30%
of the tree. The star-likeness that Milestone 4 identified is therefore **not an artefact of one algorithm**.
Every available reconstruction of these data lacks the short-distance pairs that switching-rate estimation
needs.

### What survives reconstruction choice, and what does not

| Quantity | hybrid | nj | nj_priors | ensemble | spread |
|---|---|---|---|---|---|
| median $\hat\pi_c$ | 0.343 | 0.340 | 0.338 | 0.320 | ±1% |
| **occupancy variance $V_\pi$** | **0.0092** | **0.0083** | **0.0084** | **0.0094** | **±6%** |
| median $\hat k_c$ | 0.663 | 0.664 | 0.942 | 0.785 | ±20% |
| **rate variance $V_{\log k}$** | **0.634** | **0.837** | **1.808** | **0.277** | **6.5-fold** |

**Clone-associated occupancy is robust to reconstruction**: the variance moves by 6% across trees whose heights
differ twofold. This is the strongest evidence yet that the occupancy finding is a property of the data rather
than of an analysis choice, and it satisfies Milestone 5 success criterion 5.

**Clone-associated switching rate is not robust**: $V_{\log k}$ ranges 0.28–1.81 across equally defensible
trees — a 6.5-fold spread that dwarfs the biological effect it purports to measure. Notably the **ensemble**,
which integrates over the three trees instead of committing to one, gives the *smallest* value (0.277),
consistent with much of the apparent rate heterogeneity being reconstruction noise that averaging removes.

This is a direct answer to the Milestone 4 question. The conclusion that switching rates are unidentifiable in
these data **survives, and is strengthened**: not only do the fitted rates lie outside the covered band, the
rate-heterogeneity estimate itself is dominated by which tree one happens to use.

## 2. Synthetic: how much star-likeness does reconstruction create?

E17 compares the true sampled genealogy against reconstructions of simulated recordings of that same
genealogy; results in `milestone5_information_loss.md`. The headline geometric comparison is reported there,
along with the decomposition of estimator bias into branch-length miscalibration and topology error.

## 3. Methodological consequences
1. **Report the reconstruction, and preferably more than one.** A switching-rate estimate from a single tree is
   not reproducible across equally defensible trees of the same data.
2. **Integrate over trees where the deposit allows it.** The ensemble estimate is both more conservative and
   less dependent on an arbitrary choice.
3. **Occupancy may be reported from a single tree.** Its stability across reconstructions is good enough that
   the choice does not matter at the precision available.
4. **Do not convert "resolved branches" into time** without a calibration: the same data give heights of 15, 28
   and 35 depending on the algorithm.

## 4. Limits
The real-data ensemble rests on one animal, because the other three deposits contain a single tree each. The
three variants are not independent — they share the allele data and two of them share the NJ framework — so
their spread understates full reconstruction uncertainty rather than overstating it.
