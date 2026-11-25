"""Heritable-propensity decomposition from lineage paths, shared by the synthetic and real-data analyses.

Input per environment (one site of one animal):
    paths[i]   tuple of shared-ancestor identifiers, root to leaf (edited nodes, or shared indel prefix)
    state[i]   0/1
    clone[i]   clone label
Output: the two MRCA-depth covariance curves (clone-centred and site-centred) and the two-step (k, A, B) fit,
with a clone-bootstrap CI for B. See theory/branching_model.md, Proposition T2.
"""
import numpy as np
from .ctmc import two_step_plateau


def shared_depth(a, b):
    s = 0
    for x, y in zip(a, b):
        if x != y:
            return s
        s += 1
    return s


def pair_table(paths, state, clone, max_pairs_per_clone=60_000, min_clone=20, rng=None):
    """All within-clone leaf pairs: MRCA depth, the two states, and the clone."""
    rng = rng or np.random.default_rng(0)
    out_d, out_x, out_y, out_c = [], [], [], []
    for c in np.unique(clone):
        idx = np.where(clone == c)[0]
        if len(idx) < min_clone:
            continue
        n = len(idx)
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        if len(pairs) > max_pairs_per_clone:
            pairs = [pairs[t] for t in rng.choice(len(pairs), max_pairs_per_clone, replace=False)]
        for i, j in pairs:
            out_d.append(shared_depth(paths[idx[i]], paths[idx[j]]))
            out_x.append(state[idx[i]]); out_y.append(state[idx[j]]); out_c.append(c)
    return (np.array(out_d), np.array(out_x, float), np.array(out_y, float), np.array(out_c))


def curves(d, x, y, c, state, clone, min_pairs=40):
    """Site-centred and clone-centred covariance curves by MRCA depth."""
    site_mu = float(np.mean(state))
    cmu = {k: float(np.mean(state[clone == k])) for k in np.unique(clone)}
    pr_site = (x - site_mu) * (y - site_mu)
    mu_c = np.array([cmu[k] for k in c])
    pr_clone = (x - mu_c) * (y - mu_c)
    out = []
    for pr in (pr_clone, pr_site):
        rows = []
        for h in np.unique(d):
            m = d == h
            if m.sum() >= min_pairs:
                rows.append((h, pr[m].mean(), m.sum()))
        out.append(np.array(rows) if rows else None)
    return out[0], out[1]


def decompose(paths, state, clone, n_boot=200, rng=None, min_clone=20):
    """Point estimate plus a clone-bootstrap CI for the heritable-propensity plateau B."""
    rng = rng or np.random.default_rng(0)
    state = np.asarray(state, float); clone = np.asarray(clone)
    d, x, y, c = pair_table(paths, state, clone, min_clone=min_clone, rng=rng)
    if len(d) == 0:
        return None
    cw, cs = curves(d, x, y, c, state, clone)
    fit = two_step_plateau(cw, cs)
    if fit is None:
        return None
    clones = np.unique(c)
    boot = []
    for _ in range(n_boot):
        pick = rng.choice(clones, len(clones), replace=True)
        m = np.concatenate([np.where(c == k)[0] for k in pick])
        keep = np.isin(clone, pick)
        cw_b, cs_b = curves(d[m], x[m], y[m], c[m], state[keep], clone[keep])
        f = two_step_plateau(cw_b, cs_b)
        if f:
            boot.append((f["B"], f["A"], f["k"], f["frac_plateau"]))
    boot = np.array(boot) if boot else np.zeros((0, 4))
    res = dict(fit, n_clones=len(clones), n_pairs=len(d), n_cells=len(state), mean_state=float(state.mean()),
               clone_frac_var=float(np.var([np.mean(state[clone == k]) for k in clones], ddof=1))
               if len(clones) > 1 else np.nan)
    if len(boot) >= 20:
        res.update(B_lo=float(np.percentile(boot[:, 0], 2.5)), B_hi=float(np.percentile(boot[:, 0], 97.5)),
                   A_lo=float(np.percentile(boot[:, 1], 2.5)), A_hi=float(np.percentile(boot[:, 1], 97.5)),
                   k_lo=float(np.percentile(boot[:, 2], 2.5)), k_hi=float(np.percentile(boot[:, 2], 97.5)),
                   frac_lo=float(np.nanpercentile(boot[:, 3], 2.5)),
                   frac_hi=float(np.nanpercentile(boot[:, 3], 97.5)), n_boot=len(boot))
    return res


def curve_frame(paths, state, clone, min_clone=20, rng=None):
    """The two curves as a tidy array for plotting: columns depth, clone_centred, site_centred, n_pairs."""
    state = np.asarray(state, float); clone = np.asarray(clone)
    d, x, y, c = pair_table(paths, state, clone, min_clone=min_clone, rng=rng or np.random.default_rng(0))
    cw, cs = curves(d, x, y, c, state, clone)
    if cw is None or cs is None:
        return None
    m = {int(h): (v, n) for h, v, n in cw}
    rows = [(int(h), m.get(int(h), (np.nan, 0))[0], v, n) for h, v, n in cs]
    return np.array(rows)
