"""Nested 2-state CTMC models M0-M4 on real within-clone lineage trees.

Each clone contributes one tree. Edges are grouped: group 0 = internal edges and terminal edges into primary-tumour
leaves; group 1 = terminal edges into metastatic leaves. A clone's root state is drawn from that clone's stationary
fraction.

  M0  one stationary fraction pi and one rate k                              2 params
  M1  clone-specific pi_c, shared k                                      1 + C params
  M2  shared pi, plus a logit offset on metastatic terminal edges            3 params
  M4  clone-specific pi_c plus the metastatic offset                     2 + C params

M3 (differential growth and founder bottlenecks) is not expressible as a fixed-tree CTMC likelihood: the tree is
conditioned on, so growth differences do not enter. It is assessed by simulation instead.
"""
import numpy as np
from scipy.optimize import minimize
from ..lineage_inference.ctmc import FastTree

LOGIT = lambda p: np.log(p / (1 - p))
EXPIT = lambda x: 1 / (1 + np.exp(-x))


def build(parent, leaf_state, is_met_leaf):
    group = np.where(np.asarray(is_met_leaf), 1, 0)
    return FastTree(parent, leaf_state, group)


def _ll(trees, pis, k, met_off):
    tot = 0.0
    for t, pi in zip(trees, pis):
        p1 = float(np.clip(pi, 1e-4, 1 - 1e-4))
        p2 = float(np.clip(EXPIT(LOGIT(p1) + met_off), 1e-4, 1 - 1e-4))
        tot += t.loglik(p1, [(p1, k), (p2, k)])
    return tot


def fit(trees, model, maxiter=1500):
    C = len(trees)
    if model == "M0":
        x0 = [0.0, -0.7]
        f = lambda x: -_ll(trees, [EXPIT(x[0])] * C, np.exp(x[1]), 0.0)
        np_ = 2
    elif model == "M1":
        x0 = [-0.7] + [0.0] * C
        f = lambda x: -_ll(trees, EXPIT(np.asarray(x[1:])), np.exp(x[0]), 0.0)
        np_ = 1 + C
    elif model == "M2":
        x0 = [0.0, -0.7, 0.0]
        f = lambda x: -_ll(trees, [EXPIT(x[0])] * C, np.exp(x[1]), x[2])
        np_ = 3
    elif model == "M4":
        x0 = [-0.7, 0.0] + [0.0] * C
        f = lambda x: -_ll(trees, EXPIT(np.asarray(x[2:])), np.exp(x[0]), x[1])
        np_ = 2 + C
    else:
        raise ValueError(model)
    best = None
    for jig in (0.0, 0.4):
        s = np.array(x0) + jig * np.linspace(-1, 1, len(x0))
        r = minimize(f, s, method="Nelder-Mead",
                     options=dict(maxiter=maxiter * max(1, len(x0) // 3), fatol=1e-6, xatol=1e-4))
        if best is None or r.fun < best.fun:
            best = r
    ll = -best.fun
    return dict(model=model, loglik=ll, npar=np_, aic=2 * np_ - 2 * ll, x=best.x.tolist())
