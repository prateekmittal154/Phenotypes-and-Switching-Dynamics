# Information across the observation hierarchy (Propositions V1–V3)

Notation. $\theta=(\pi,k)$ are the switching parameters. $G_0$ is the true genealogy of all cells, $G_1$ the
genealogy induced on the sampled cells, $B$ the recorded allele matrix, $\hat G = f(B)$ a reconstructed tree
from a fixed algorithm $f$. $X_S$ denotes the observed phenotypic states of the sampled cells, $X_{\text{all}}$
the states of every cell.

## The key structural fact
**The recorder writes marks that depend on the genealogy, not on $\theta$.** Editing happens at cell divisions
and is independent of phenotype; likewise sampling and reconstruction do not consult $\theta$. Formally,

$$\theta \;\longrightarrow\; (G_0, X_{\text{all}}) \;\longrightarrow\; (G_1, X_S) \;\longrightarrow\; (B, X_S)
\;\longrightarrow\; (\hat G, X_S)$$

is a Markov chain in which every arrow after the first is a **parameter-independent** stochastic or
deterministic map.

## Proposition V1 (information-processing chain)
For mutual information with any prior on $\theta$,
$$I\big(\theta; X_{\text{all}}, G_0\big) \;\ge\; I\big(\theta; X_S, G_1\big) \;\ge\; I\big(\theta; X_S, B\big)
\;\ge\; I\big(\theta; X_S, \hat G\big).$$

*Proof.* Each step is a $\theta$-independent channel applied to the conditioning variable, so the data-processing
inequality applies to each link of the Markov chain. Sampling cells maps $(G_0,X_{\text{all}})$ to
$(G_1,X_S)$ without reference to $\theta$; the recorder maps $G_1$ (together with the unsampled structure) to
$B$ without reference to $\theta$; and $\hat G = f(B)$ is a deterministic function of $B$. ∎

**Three consequences that matter for practice.**

1. **A single reconstructed tree can never beat the barcodes.** $I(\theta; X_S,\hat G) \le I(\theta; X_S,B)$,
   with equality only if $\hat G$ is a sufficient statistic of $B$ for $\theta$, which it is not in general.
   This is the formal justification for integrating over trees rather than committing to one, and it bounds
   what approach C can deliver.
2. **No reconstruction method can recover information the recorder never wrote.** The bound
   $I(\theta; X_S,B) \le I(\theta; X_S,G_1)$ is a property of the recorder and the sampling, not of the
   algorithm. If barcodes saturate before the divisions of interest, the gap is permanent.
3. **The chain gives an audit.** Measuring the loss at each link separates *biological sampling* loss
   ($G_0 \to G_1$), *recording* loss ($G_1 \to B$) and *reconstruction* loss ($B \to \hat G$), which is exactly
   the decomposition Milestone 5 set out to produce.

## Proposition V2 (what the inequality does **not** say)
V1 is a statement about mutual information for a fixed prior, or equivalently about expected information. Three
things do not follow.

- **Fisher information at a point is not monotone in the same way for a misspecified model.** When the analyst
  conditions on a wrong tree, the likelihood is misspecified; the curvature of that wrong likelihood can exceed
  the curvature of the correct one. A narrower interval from a wrong tree is overconfidence, not information.
- **Estimator stability is not information.** An estimator computed on $\hat G$ may have smaller variance than
  one computed on $G_1$ while being badly biased. Reporting a reduction in variance as "recovered information"
  is therefore invalid; bias and coverage must be reported together. This is why the Milestone 5 experiments
  report paired log-ratios and coverage, not variance alone.
- **The inequality is about the experiments, not about any particular estimator.** A poor estimator on $G_1$ can
  easily be worse than a good estimator on $\hat G$; V1 bounds only what is achievable.

## Proposition V3 (saturation gives a hard ceiling)
Let the recorder have $S$ sites, each editable once, with the per-division edit probability $\mu$. After $g$
generations the expected number of still-unedited sites is $S(1-\mu)^g$, so the expected number of divisions
after which no further marks can accrue is finite and the recorder becomes uninformative about coalescences
more recent than
$$g^\ast \;\approx\; \frac{\log S}{-\log(1-\mu)} .$$
Beyond $g^\ast$, $B$ carries no information distinguishing recent ancestry, so by V1 neither does any
reconstruction. **Fast switching is exactly what recent coalescences would reveal**, which is why saturation
caps the detectable rate rather than merely adding noise.

In the simulation this is directly observable: the low-capacity setting reaches 97% saturation within ten
generations, and its reconstructions cannot resolve the recent structure, whereas the high-capacity setting
reaches 73%.

## How these are used
- V1 motivates comparing four estimators at four levels of the hierarchy (experiment E17) and gives the
  direction each comparison must respect; a violation would indicate an implementation error.
- V2 is the reason the headline metric is paired bias plus interval coverage, not estimator variance.
- V3 connects the hierarchy to the Milestone 4 capacity law: the law's ceiling and the saturation ceiling are
  the same constraint seen from two directions.

## Proposition V4 (rate–time confounding, and why a single scale factor cannot fix it)

**Statement.** For a two-state CTMC on a fixed tree, multiplying every branch length by $c>0$ and dividing $k$
by $c$ leaves the likelihood unchanged:
$$\mathcal{L}(\pi, k;\ \{L_e\}) \;=\; \mathcal{L}(\pi,\ k/c;\ \{cL_e\}).$$

*Proof.* The likelihood depends on branch lengths only through the transition matrices $P(L_e)$, whose entries
depend on $(k, L_e)$ only via the product $kL_e$ (Lemma T1 and the two-state exponential). Scaling every $L_e$
by $c$ and $k$ by $1/c$ leaves every product invariant. ∎

**Consequence 1 — unit edges do not estimate a rate, they estimate a rate in unknown units.** An analyst who
sets $L_e \equiv 1$ is implicitly asserting that one resolved branch equals one unit of time. If the true mean
is $\bar{L}$ divisions per resolved branch, the estimate is inflated by approximately $\bar{L}$. This is the
mechanism behind the 1.83-fold bias measured in Milestone 5, and it explains why the bias is common to every
reconstruction: all of them share the convention, not the error.

**Consequence 2 — a single global scale is not a valid correction.** Setting $c = D / h$, with $D$ the number
of generations and $h$ the tree height, makes only the *deepest* root-to-leaf path equal $D$. In an unbalanced
tree most leaves are shallower, so their implied paths are too short and the rate is inflated further.
Experiment E20 confirms this: the naive global scale *worsened* both bias (+1.85 to +2.11 in log units) and
coverage (0.484 to 0.328).

**Consequence 3 — the correct calibration is ultrametric.** All sampled cells are the same age, so every
root-to-leaf path must equal $D$. Assigning node times
$$t(v) \;=\; D\,\frac{\mathrm{depth}(v)}{\mathrm{depth}(v) + \mathrm{maxdepth\ below}(v)}$$
makes the tree ultrametric and uses only information the experiment already has — the elapsed number of
generations. This reduced the bias from +1.85 to +1.37 in log units where the naive scale made it worse.

**Consequence 4 — the residual is irreducible without a clock.** Even ultrametric calibration leaves bias,
because $D$ itself is usually known only approximately and the division count varies between lineages.
Treating $D$ as uncertain and integrating over it (variant V6) widens intervals and improves coverage from
0.48 to 0.73, but does not reach nominal: by Proposition V1 the missing information cannot be manufactured.
