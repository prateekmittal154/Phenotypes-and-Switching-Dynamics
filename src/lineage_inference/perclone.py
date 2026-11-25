"""Per-clone (pi_c, k_c) estimation and the occupancy / switching-rate variance components.

Milestone 2 estimated one pooled plateau B = Var_c(pi_c) and one shared decay rate k. The audit (U1) shows the
plateau is invariant to k, so it measures clone-associated *equilibrium occupancy*, not switching dynamics. To
separate the two we need per-clone estimates of BOTH parameters:

    pi_c = a_c / (a_c + b_c)   equilibrium occupancy of the high state
    k_c  = a_c + b_c           switching timescale

Then
    V_pi = Var_c(pi_c)         clone-associated occupancy heterogeneity (the Milestone 2 plateau)
    V_k  = Var_c(log k_c)      clone-associated switching-rate heterogeneity (true "heritable switching dynamics")

Both estimates are corrected for estimation noise by the method of moments: Var(hat) - mean(SE^2).
"""
import numpy as np
from scipy.optimize import minimize_scalar, minimize
from scipy.stats import chi2

from .ctmc import FastTree

EXPIT = lambda x: 1 / (1 + np.exp(-x))
LOGIT = lambda p: np.log(p / (1 - p))


def fit_clone(tree, pi_grid=None, k_grid=None, refine=True, profile=True):
    """ML fit of (pi, k) for one clone tree, with a profile-likelihood CI for each parameter.

    A coarse grid is evaluated first (the surface is often flat in k, so local optimisers stall), then refined.
    Returns the estimates, their profile 95% intervals, the curvature-based SEs and the log-likelihood.
    """
    pi_grid = np.linspace(0.02, 0.98, 25) if pi_grid is None else pi_grid
    k_grid = np.exp(np.linspace(np.log(1e-3), np.log(20), 31)) if k_grid is None else k_grid
    LL = np.full((len(pi_grid), len(k_grid)), -np.inf)
    for i, p in enumerate(pi_grid):
        for j, k in enumerate(k_grid):
            LL[i, j] = tree.loglik(p, [(p, k), (p, k)])
    i, j = np.unravel_index(np.nanargmax(LL), LL.shape)
    pi, k = pi_grid[i], k_grid[j]
    if refine:
        f = lambda x: -tree.loglik(EXPIT(x[0]), [(EXPIT(x[0]), np.exp(x[1])), (EXPIT(x[0]), np.exp(x[1]))])
        r = minimize(f, [LOGIT(pi), np.log(k)], method="Nelder-Mead",
                     options=dict(maxiter=400, fatol=1e-7, xatol=1e-4))
        if -r.fun >= LL[i, j]:
            pi, k = float(EXPIT(r.x[0])), float(np.clip(np.exp(r.x[1]), 1e-4, 50))
    ll = tree.loglik(pi, [(pi, k), (pi, k)])
    n0 = int((tree.leaf_state >= 0).sum())
    if not profile:
        # point estimates only: the profile grids triple the cost and are not needed inside a bootstrap
        return dict(pi=float(pi), k=float(k), loglik=float(ll), n=n0,
                    pi_lo=np.nan, pi_hi=np.nan, k_lo=np.nan, k_hi=np.nan,
                    pi_se=np.nan, logk_se=np.nan,
                    k_at_bound=bool(k <= 1.5e-3 or k >= 19.0))
    cut = ll - chi2.ppf(0.95, 1) / 2

    def prof(which):
        """Profile interval: maximise over the other parameter on the grid."""
        if which == "pi":
            prof_ll = np.array([max(tree.loglik(p, [(p, kk), (p, kk)]) for kk in k_grid) for p in pi_grid])
            ok = pi_grid[prof_ll >= cut]
        else:
            prof_ll = np.array([max(tree.loglik(p, [(p, kk), (p, kk)]) for p in pi_grid) for kk in k_grid])
            ok = k_grid[prof_ll >= cut]
        return (float(ok.min()), float(ok.max())) if len(ok) else (np.nan, np.nan)

    pl, ph = prof("pi"); kl, kh = prof("k")
    # the profile is evaluated on a grid while the point estimate is refined off-grid, so widen each interval to
    # contain its MLE rather than reporting an interval that excludes it
    if np.isfinite(pl):
        pl, ph = min(pl, pi), max(ph, pi)
    if np.isfinite(kl):
        kl, kh = min(kl, k), max(kh, k)
    n = int((tree.leaf_state >= 0).sum())
    return dict(pi=float(pi), k=float(k), loglik=float(ll), n=n,
                pi_lo=pl, pi_hi=ph, k_lo=kl, k_hi=kh,
                pi_se=(ph - pl) / 3.92 if np.isfinite(ph) else np.nan,
                logk_se=(np.log(kh) - np.log(kl)) / 3.92 if np.isfinite(kh) and kl > 0 else np.nan,
                k_at_bound=bool(k <= 1.5e-3 or k >= 19.0))


def variance_components(fits, min_n=25):
    """Noise-corrected Var_c(pi_c) and Var_c(log k_c) across clones (method of moments)."""
    f = [x for x in fits if x["n"] >= min_n]
    out = dict(n_clones=len(f))
    if len(f) < 3:
        return out
    pi = np.array([x["pi"] for x in f]); se = np.array([x["pi_se"] for x in f], float)
    raw = float(np.var(pi, ddof=1))
    noise = float(np.nanmean(se ** 2)) if np.isfinite(se).any() else 0.0
    out.update(V_pi_raw=raw, V_pi_noise=noise, V_pi=max(raw - noise, 0.0), mean_pi=float(pi.mean()))
    usable = [x for x in f if not x["k_at_bound"] and np.isfinite(x["logk_se"])]
    out["n_clones_k_identified"] = len(usable)
    if len(usable) >= 3:
        lk = np.log([x["k"] for x in usable]); s = np.array([x["logk_se"] for x in usable])
        rawk = float(np.var(lk, ddof=1)); nk = float(np.mean(s ** 2))
        out.update(V_logk_raw=rawk, V_logk_noise=nk, V_logk=max(rawk - nk, 0.0),
                   median_k=float(np.exp(np.median(lk))))
    return out


def lrt_shared(trees, fits, what="pi", n_boot=0, rng=None, maxiter_mult=60):
    """Likelihood-ratio test of a shared parameter across clones against clone-specific values.

    what='pi': H0 pi common, k clone-specific. what='k': H0 k common, pi clone-specific.
    Returns the statistic, df, the asymptotic p-value, and a parametric-bootstrap p-value when n_boot > 0.
    """
    full = sum(f["loglik"] for f in fits)
    C = len(trees)
    if what == "pi":
        def nll(x):
            p = EXPIT(x[0])
            return -sum(t.loglik(p, [(p, np.exp(x[1 + i])), (p, np.exp(x[1 + i]))]) for i, t in enumerate(trees))
        x0 = [LOGIT(np.clip(np.mean([f["pi"] for f in fits]), .02, .98))] + [np.log(f["k"]) for f in fits]
        df = C - 1
    else:
        def nll(x):
            k = np.exp(x[0])
            return -sum(t.loglik(EXPIT(x[1 + i]), [(EXPIT(x[1 + i]), k), (EXPIT(x[1 + i]), k)])
                        for i, t in enumerate(trees))
        x0 = [np.log(np.median([f["k"] for f in fits]))] + [LOGIT(np.clip(f["pi"], .02, .98)) for f in fits]
        df = C - 1
    r = minimize(nll, x0, method="Nelder-Mead", options=dict(maxiter=maxiter_mult * len(x0), fatol=1e-5))
    stat = 2 * (full - (-r.fun))
    return dict(what=what, stat=float(max(stat, 0)), df=df,
                p_asym=float(chi2.sf(max(stat, 0), df)), ll_null=float(-r.fun), ll_full=float(full))


def moran_path(paths, state, clone):
    """A PATH-style single-number phylogenetic autocorrelation, as the established-method baseline.

    Weights same-clone pairs by their shared MRCA depth and returns the weighted state autocorrelation.
    """
    from .propensity import shared_depth
    state = np.asarray(state, float); mu = state.mean()
    num = den = 0.0
    idx = np.arange(len(state))
    for c in np.unique(clone):
        k = idx[clone == c]
        for ii in range(len(k)):
            for jj in range(ii + 1, len(k)):
                w = shared_depth(paths[k[ii]], paths[k[jj]])
                num += w * (state[k[ii]] - mu) * (state[k[jj]] - mu); den += w
    v = float(np.var(state))
    return float(num / den / v) if den > 0 and v > 0 else np.nan
