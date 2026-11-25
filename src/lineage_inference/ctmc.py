"""Two-state CTMC on lineage trees: likelihood (Felsenstein pruning), MRCA-depth covariance, and estimators.

State 0 = low, 1 = high. The generator is
    Q = [[-a, a], [c, -c]],
so the stationary fraction of state 1 is pi1 = a / (a + c) and the total rate is k = a + c. The transition matrix
over an edge of length L is
    P(L) = [[c + a E, a - a E], [c - c E, a + c E]] / (a + c),  with E = exp(-k L).

The parameters (pi1, k) are used throughout, because the covariance decomposition (Proposition T2) is natural in them.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp


def trans(pi1, k, L):
    """P(L)[i, j] for a 2-state CTMC, vectorised over L. Returns an array of shape (len(L), 2, 2)."""
    L = np.atleast_1d(np.asarray(L, float))
    E = np.exp(-k * L)
    P = np.empty((len(L), 2, 2))
    P[:, 0, 0] = 1 - pi1 * (1 - E); P[:, 0, 1] = pi1 * (1 - E)
    P[:, 1, 0] = (1 - pi1) * (1 - E); P[:, 1, 1] = 1 - (1 - pi1) * (1 - E)
    return P


class Tree:
    """A rooted tree held as arrays, with children listed before parents (post-order).

    parent[i] is the index of i's parent, or -1 for the root. elen[i] is the edge length above i. leaf_state[i] is
    0/1 for a leaf, or -1 for an internal node. group[i] identifies the terminal environment of a leaf, used to let
    the terminal edge carry organ-specific rates.
    """

    def __init__(self, parent, elen, leaf_state, group=None):
        self.parent = np.asarray(parent, int); self.elen = np.asarray(elen, float)
        self.leaf_state = np.asarray(leaf_state, int)
        self.group = np.zeros(len(parent), int) if group is None else np.asarray(group, int)
        self.n = len(parent)
        self.order = np.argsort(-self._depth())          # deepest first = post-order
        self.children = [[] for _ in range(self.n)]
        for i, p in enumerate(self.parent):
            if p >= 0:
                self.children[p].append(i)
        self.root = int(np.where(self.parent < 0)[0][0])

    def _depth(self):
        d = np.zeros(self.n, int)
        for i in range(self.n):
            j, c = i, 0
            while self.parent[j] >= 0:
                j = self.parent[j]; c += 1
            d[i] = c
        return d

    def loglik(self, pi1_root, params):
        """params[g] = (pi1, k) for edges whose child has group g. Returns the log-likelihood of the leaf states."""
        L = np.zeros((self.n, 2))                     # log partial likelihoods
        for i in self.order:
            if self.leaf_state[i] >= 0:
                L[i] = [-np.inf, -np.inf]; L[i, self.leaf_state[i]] = 0.0
            for ch in self.children[i]:
                pi1, k = params[self.group[ch]]
                P = trans(pi1, k, self.elen[ch])[0]
                with np.errstate(divide="ignore"):
                    L[i] += logsumexp(np.log(P) + L[ch][None, :], axis=1)
        pri = np.log([1 - pi1_root, pi1_root])
        return float(logsumexp(pri + L[self.root]))


def fit_trees(trees, n_groups=1, share_root=True, bounds_k=(1e-3, 20.0), x0=None):
    """Maximum likelihood over a list of Tree objects. Parameters: pi1_root, then (pi1, k) per group."""
    def unpack(x):
        pr = 1 / (1 + np.exp(-x[0]))
        pars = []
        for g in range(n_groups):
            pars.append((1 / (1 + np.exp(-x[1 + 2 * g])), np.exp(x[2 + 2 * g])))
        return pr, pars

    def nll(x):
        pr, pars = unpack(x)
        if any(not (bounds_k[0] <= k <= bounds_k[1]) for _, k in pars):
            return 1e9
        return -sum(t.loglik(pr, pars) for t in trees)

    best = None
    for s in ([x0] if x0 is not None else [np.array([0.0] + [0.0, 0.0] * n_groups),
                                           np.array([-1.0] + [-1.0, 1.0] * n_groups),
                                           np.array([1.0] + [1.0, -1.0] * n_groups)]):
        r = minimize(nll, s, method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-6, fatol=1e-8))
        if best is None or r.fun < best.fun:
            best = r
    pr, pars = unpack(best.x)
    return dict(loglik=-best.fun, pi1_root=pr, params=pars, npar=1 + 2 * n_groups, x=best.x)


# ---------------------------------------------------------------- covariance by MRCA depth

def mrca_depth_pairs(paths, states, max_pairs=400_000, rng=None):
    """For all (or a random subset of) leaf pairs, return the MRCA depth and the product of centred states.

    paths[i] is the tuple of visible ancestor nodes from the root to leaf i. The MRCA depth is the length of the
    shared prefix, which is the number of shared edits.
    """
    n = len(paths)
    idx = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(idx) > max_pairs:
        rng = rng or np.random.default_rng(0)
        idx = [idx[t] for t in rng.choice(len(idx), max_pairs, replace=False)]
    m = np.mean(states)
    d, pr = [], []
    for i, j in idx:
        a, b = paths[i], paths[j]
        s = 0
        for x, y in zip(a, b):
            if x != y:
                break
            s += 1
        d.append(s); pr.append((states[i] - m) * (states[j] - m))
    return np.array(d), np.array(pr)


def cov_curve(depth, prod, min_n=30):
    """Mean centred product per MRCA depth, with the pair count."""
    out = []
    for h in np.unique(depth):
        m = depth == h
        if m.sum() >= min_n:
            out.append((int(h), float(prod[m].mean()), int(m.sum()),
                        float(prod[m].std(ddof=1) / np.sqrt(m.sum()))))
    return np.array(out) if out else np.zeros((0, 4))


def fit_plateau(curve, w=True, kbounds=(0.02, 8.0)):
    """Weighted least-squares fit of C(h) = A exp(-k (D - h)) + B, with A >= 0 and B >= 0.

    h is the MRCA edit depth, and D = max(h), so (D - h) increases with lineage distance between the pair.
    By Proposition T2, A = E_c[pi_c (1 - pi_c)] * exp(-2 k0 L0) (state memory, decays with distance) and
    B = Var_c(pi_c) (heritable propensity, does not decay). Both are variances, hence the non-negativity
    constraints: without them the fit is unidentified when the decay is shallow.
    """
    from scipy.optimize import least_squares
    if len(curve) < 3:
        return None
    h, c, n = curve[:, 0], curve[:, 1], curve[:, 2]
    D = h.max(); x = D - h
    wt = np.sqrt(n) if w else np.ones_like(n)

    def res(p):
        A, k, B = p
        return (A * np.exp(-k * x) + B - c) * wt

    lo = [0.0, kbounds[0], 0.0]; hi = [2.0, kbounds[1], 1.0]
    best = None
    for s in ((max(c.max() - min(c.min(), 0), 1e-3), 0.5, max(c.min(), 0.0)),
              (max(c.max(), 1e-3), 2.0, 0.0),
              (1e-3, 0.1, max(c.mean(), 0.0))):
        s = np.clip(s, lo, hi)
        r = least_squares(res, s, bounds=(lo, hi), max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    A, k, B = best.x
    return dict(A=float(A), k=float(k), B=float(B), cost=float(best.cost),
                frac_plateau=float(B / (A + B)) if (A + B) > 1e-12 else np.nan)


def two_step_plateau(curve_within_clone, curve_within_site, kbounds=(0.02, 8.0)):
    """Two-step estimator of the heritable-propensity plateau B (Proposition T2, estimator form).

    Step 1. Centring each clone at its own mean removes the between-clone term, so the within-clone curve is
            C_w(h) = A exp(-k (D - h)) + d, where d is the (negative) finite-clone centring bias. Fitting this
            identifies the decay rate k. Without it, a flat site-centred curve is explained equally well by a
            plateau (B > 0, A = 0) and by perfectly heritable state (A > 0, k -> 0): see counterexample C7.
    Step 2. With k fixed, the site-centred curve C(h) = A exp(-k (D - h)) + B is linear in (A, B), so B follows
            from a weighted non-negative least-squares fit.

    Returns k, A, B and the plateau fraction, or None when either curve has fewer than three depth classes.
    """
    from scipy.optimize import least_squares, nnls
    if curve_within_clone is None or curve_within_site is None:
        return None
    if len(curve_within_clone) < 3 or len(curve_within_site) < 3:
        return None
    hw, cw, nw = curve_within_clone[:, 0], curve_within_clone[:, 1], curve_within_clone[:, 2]
    D = max(hw.max(), curve_within_site[:, 0].max())
    xw = D - hw; ww = np.sqrt(nw)

    def res(p):
        A, k, d = p
        return (A * np.exp(-k * xw) + d - cw) * ww

    best = None
    for s in ((max(cw.max() - cw.min(), 1e-3), 0.5, cw.min()), (max(cw.max(), 1e-3), 2.0, 0.0)):
        r = least_squares(res, s, bounds=([0.0, kbounds[0], -1.0], [2.0, kbounds[1], 1.0]), max_nfev=20000)
        if best is None or r.cost < best.cost:
            best = r
    k = float(best.x[1])
    hs, cs, ns = curve_within_site[:, 0], curve_within_site[:, 1], curve_within_site[:, 2]
    ws = np.sqrt(ns)
    M = np.column_stack([np.exp(-k * (D - hs)), np.ones_like(hs)]) * ws[:, None]
    coef, _ = nnls(M, cs * ws)
    A, B = float(coef[0]), float(coef[1])
    return dict(k=k, A=A, B=B, k_step1_A=float(best.x[0]), bias=float(best.x[2]),
                frac_plateau=float(B / (A + B)) if (A + B) > 1e-12 else np.nan)


class FastTree:
    """Level-vectorised version of Tree for large real trees (all edge lengths equal 1).

    Nodes are grouped by depth and processed deepest-first with numpy, which avoids the per-node Python loop.
    group[i] selects the parameter set used on the edge above node i.
    """

    def __init__(self, parent, leaf_state, group=None):
        self.parent = np.asarray(parent, int)
        self.leaf_state = np.asarray(leaf_state, int)
        self.n = len(self.parent)
        self.group = np.zeros(self.n, int) if group is None else np.asarray(group, int)
        d = np.zeros(self.n, int)
        for i in range(self.n):
            j, c = i, 0
            while self.parent[j] >= 0:
                j = self.parent[j]; c += 1
            d[i] = c
        self.depth = d
        self.levels = [np.where(d == h)[0] for h in range(d.max(), -1, -1)]
        self.root = int(np.where(self.parent < 0)[0][0])
        self.is_leaf = ~np.isin(np.arange(self.n), self.parent[self.parent >= 0])

    def loglik(self, pi1_root, params):
        L = np.zeros((self.n, 2))
        leaf = self.leaf_state >= 0
        L[leaf] = -np.inf
        L[np.where(leaf)[0], self.leaf_state[leaf]] = 0.0
        logP = {}
        for g, (pi1, k) in enumerate(params):
            with np.errstate(divide="ignore"):
                logP[g] = np.log(trans(pi1, k, 1.0)[0])
        for lv in self.levels:
            nodes = lv[self.parent[lv] >= 0]
            if len(nodes) == 0:
                continue
            msg = np.empty((len(nodes), 2))
            for g in np.unique(self.group[nodes]):
                m = self.group[nodes] == g
                msg[m] = logsumexp(logP[g][None, :, :] + L[nodes[m]][:, None, :], axis=2)
            np.add.at(L, self.parent[nodes], msg)
        pri = np.log([max(1 - pi1_root, 1e-300), max(pi1_root, 1e-300)])
        return float(logsumexp(pri + L[self.root]))


def fit_fast(trees, n_groups=1, maxiter=800, starts=None):
    """ML fit over FastTree objects. Parameters: logit(pi1_root), then logit(pi1), log(k) per group."""
    def unpack(x):
        pr = 1 / (1 + np.exp(-x[0]))
        return pr, [(1 / (1 + np.exp(-x[1 + 2 * g])), float(np.clip(np.exp(x[2 + 2 * g]), 1e-3, 20)))
                    for g in range(n_groups)]

    def nll(x):
        pr, pars = unpack(x)
        return -sum(t.loglik(pr, pars) for t in trees)

    starts = starts or [np.array([0.0] + [0.0, -0.7] * n_groups), np.array([-0.5] + [-0.5, 0.5] * n_groups)]
    best = None
    for s in starts:
        r = minimize(nll, s, method="Nelder-Mead", options=dict(maxiter=maxiter, fatol=1e-6, xatol=1e-4))
        if best is None or r.fun < best.fun:
            best = r
    pr, pars = unpack(best.x)
    npar = 1 + 2 * n_groups
    return dict(loglik=-best.fun, pi1_root=pr, params=pars, npar=npar, aic=2 * npar + 2 * best.fun, x=best.x)


class VarTree:
    """Level-vectorised tree likelihood with arbitrary per-edge lengths.

    Needed to separate two distinct reconstruction errors: getting the topology wrong, and getting the branch
    lengths wrong. FastTree assumes every edge has length 1, which is what an analyst implicitly assumes when
    treating edit depth as a clock.
    """

    def __init__(self, parent, leaf_state, elen=None):
        self.parent = np.asarray(parent, int)
        self.leaf_state = np.asarray(leaf_state, int)
        self.n = len(self.parent)
        self.elen = np.ones(self.n) if elen is None else np.asarray(elen, float)
        d = np.zeros(self.n, int)
        for i in range(self.n):
            j, c = i, 0
            while self.parent[j] >= 0:
                j = self.parent[j]; c += 1
            d[i] = c
        self.levels = [np.where(d == h)[0] for h in range(d.max(), -1, -1)]
        self.root = int(np.where(self.parent < 0)[0][0])
        self.ulen = np.unique(np.round(self.elen, 6))

    def loglik(self, pi1_root, pi1, k):
        L = np.zeros((self.n, 2))
        leaf = self.leaf_state >= 0
        L[leaf] = -np.inf
        L[np.where(leaf)[0], self.leaf_state[leaf]] = 0.0
        with np.errstate(divide="ignore"):
            logP = {float(u): np.log(trans(pi1, k, max(u, 1e-9))[0]) for u in self.ulen}
        for lv in self.levels:
            nodes = lv[self.parent[lv] >= 0]
            if len(nodes) == 0:
                continue
            msg = np.empty((len(nodes), 2))
            rl = np.round(self.elen[nodes], 6)
            for u in np.unique(rl):
                m = rl == u
                msg[m] = logsumexp(logP[float(u)][None, :, :] + L[nodes[m]][:, None, :], axis=2)
            np.add.at(L, self.parent[nodes], msg)
        pri = np.log([max(1 - pi1_root, 1e-300), max(pi1_root, 1e-300)])
        return float(logsumexp(pri + L[self.root]))


def fit_var(tree, pi_grid=None, k_grid=None):
    """ML (pi, k) on a VarTree by grid search then local refinement; returns estimates and the bound flag."""
    pi_grid = np.linspace(0.04, 0.96, 24) if pi_grid is None else pi_grid
    k_grid = np.exp(np.linspace(np.log(2e-3), np.log(20), 30)) if k_grid is None else k_grid
    best = (-np.inf, np.nan, np.nan)
    for p in pi_grid:
        for k in k_grid:
            ll = tree.loglik(p, p, k)
            if ll > best[0]:
                best = (ll, p, k)
    ll, pi, k = best
    f = lambda x: -tree.loglik(1 / (1 + np.exp(-x[0])), 1 / (1 + np.exp(-x[0])), np.exp(x[1]))
    r = minimize(f, [np.log(pi / (1 - pi)), np.log(k)], method="Nelder-Mead",
                 options=dict(maxiter=300, fatol=1e-7, xatol=1e-4))
    if -r.fun >= ll:
        pi = float(1 / (1 + np.exp(-r.x[0]))); k = float(np.clip(np.exp(r.x[1]), 1e-3, 50)); ll = -r.fun
    return dict(pi=pi, k=k, loglik=ll, k_at_bound=bool(k <= 3e-3 or k >= 19.0))
