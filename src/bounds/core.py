"""Sharp bounds on the induced-transition fraction tau under Lambda-bounded founder selection.

Model: primary P, founder weights w (max/min <= Lam on supp P), pi = P*w/<P,w>, kernel K, metastasis Q = pi K.
tau = sum_x pi(x) (1 - K(x,x)) = P(founder state != metastatic state).
See theory/propositions.md for statements and theory/proofs/ for proofs.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.stats import beta as _beta


def tv(p, q):
    return 0.5 * np.abs(np.asarray(p, float) - np.asarray(q, float)).sum()


def support_mismatch(P, Q):
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    return Q[P <= 0].sum()


def _AB(P, Q, Lam, m):
    A = np.clip(Q - Lam * m * P, 0, None).sum()
    B = np.clip(m * P - Q, 0, None).sum()
    return A, B


def tau_min(P, Q, Lam, tol=1e-13):
    """Closed form (Proposition 2): min over m in [1/Lam, 1] of max(A(m), B(m)).

    A is non-increasing and B non-decreasing in m, so the minimum is at their crossing or at an endpoint.
    """
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    P, Q = P / P.sum(), Q / Q.sum()
    if np.isinf(Lam):
        return support_mismatch(P, Q)
    lo, hi = 1.0 / Lam, 1.0
    A_lo, B_lo = _AB(P, Q, Lam, lo)
    if A_lo <= B_lo:
        return B_lo
    A_hi, B_hi = _AB(P, Q, Lam, hi)
    if A_hi >= B_hi:
        return A_hi
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        A, B = _AB(P, Q, Lam, mid)
        if A > B:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return max(_AB(P, Q, Lam, 0.5 * (lo + hi)))


def tau_min_lp(P, Q, Lam, maximise=False):
    """Independent LP over (pi, m, gamma): does not use the box lemma or the coupling lemma.

    Variables: pi (d), m (1), gamma (d*d). Constraints: m P <= pi <= Lam m P, sum pi = 1,
    row sums of gamma = pi, column sums = Q. Objective: off-diagonal mass of gamma.
    """
    P, Q = np.asarray(P, float), np.asarray(Q, float)
    d = len(P)
    nv = d + 1 + d * d
    c = np.zeros(nv)
    off = np.ones((d, d)) - np.eye(d)
    c[d + 1:] = off.ravel() * (-1 if maximise else 1)
    A_ub, b_ub = [], []
    for x in range(d):
        r = np.zeros(nv); r[x] = -1; r[d] = P[x]; A_ub.append(r); b_ub.append(0)       # m P - pi <= 0
        r = np.zeros(nv); r[x] = 1; r[d] = -Lam * P[x]; A_ub.append(r); b_ub.append(0)  # pi - Lam m P <= 0
    A_eq, b_eq = [], []
    r = np.zeros(nv); r[:d] = 1; A_eq.append(r); b_eq.append(1)
    for x in range(d):
        r = np.zeros(nv); r[x] = -1
        r[d + 1 + x * d: d + 1 + (x + 1) * d] = 1; A_eq.append(r); b_eq.append(0)
    for y in range(d):
        r = np.zeros(nv); r[d + 1 + y::d][:d] = 1; A_eq.append(r); b_eq.append(Q[y])
    res = linprog(c, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq,
                  bounds=[(0, None)] * nv, method="highs")
    if res.status != 0:
        return np.nan
    return -res.fun if maximise else res.fun


def tau_max_lp(P, Q, Lam):
    return tau_min_lp(P, Q, Lam, maximise=True)


def max_selected_freq(p, Lam):
    """Largest metastatic frequency a primary state of frequency p can reach with no transitions."""
    return Lam * p / (Lam * p + 1 - p)


def p_star(q, Lam):
    """Pure selection explains metastatic frequency q of a state iff its primary frequency is >= p_star."""
    return q / (Lam * (1 - q) + q)


def cp_box(counts, alpha):
    """Simultaneous (Bonferroni over states) two-sided Clopper-Pearson box for multinomial proportions."""
    k = np.asarray(counts, float); n = k.sum(); d = len(k); a = alpha / d
    lo = np.where(k > 0, _beta.ppf(a / 2, k, n - k + 1), 0.0)
    hi = np.where(k < n, _beta.ppf(1 - a / 2, k + 1, n - k), 1.0)
    return np.nan_to_num(lo), np.nan_to_num(hi, nan=1.0)


def tau_lower_robust(kP, kQ, Lam, alpha=0.05, eps_meas=0.0):
    """Finite-sample lower confidence bound for tau (Proposition 5).

    min over P' in box_P, Q' in box_Q (each at level alpha/2) of tau_min(P', Q', Lam), minus a measurement
    allowance eps_meas (organ-specific misclassification, Proposition 6). Solved as one LP with u = m P':
    vars [pi(d), u(d), m, Q'(d), s(d)]; s >= Q' - pi; minimise sum s = TV(pi, Q').
    """
    LP_, UP_ = cp_box(kP, alpha / 2)
    LQ, UQ = cp_box(kQ, alpha / 2)
    d = len(LP_)
    iP, iU, iM, iQ, iS = 0, d, 2 * d, 2 * d + 1, 3 * d + 1
    nv = 4 * d + 1
    c = np.zeros(nv); c[iS:iS + d] = 1
    A_ub, b_ub = [], []
    def row():
        return np.zeros(nv)
    for x in range(d):
        r = row(); r[iM] = LP_[x]; r[iU + x] = -1; A_ub.append(r); b_ub.append(0)     # m L <= u
        r = row(); r[iU + x] = 1; r[iM] = -UP_[x]; A_ub.append(r); b_ub.append(0)     # u <= m U
        r = row(); r[iU + x] = 1; r[iP + x] = -1; A_ub.append(r); b_ub.append(0)      # u <= pi
        r = row(); r[iP + x] = 1; r[iU + x] = -Lam; A_ub.append(r); b_ub.append(0)    # pi <= Lam u
        r = row(); r[iQ + x] = 1; r[iP + x] = -1; r[iS + x] = -1; A_ub.append(r); b_ub.append(0)  # Q'-pi <= s
    A_eq, b_eq = [], []
    r = row(); r[iU:iU + d] = 1; r[iM] = -1; A_eq.append(r); b_eq.append(0)   # sum u = m
    r = row(); r[iP:iP + d] = 1; A_eq.append(r); b_eq.append(1)
    r = row(); r[iQ:iQ + d] = 1; A_eq.append(r); b_eq.append(1)
    bounds = [(0, None)] * (2 * d) + [(0, None)] + list(zip(LQ, UQ)) + [(0, None)] * d
    res = linprog(c, A_ub=np.array(A_ub), b_ub=b_ub, A_eq=np.array(A_eq), b_eq=b_eq,
                  bounds=bounds, method="highs")
    if res.status != 0:
        return np.nan
    return max(0.0, res.fun - eps_meas)


def tau_min_clone(Pc, Qc, Lam_in):
    """Clone-resolved sharp bound (Proposition 7). Pc, Qc: (C, d) joint clone x state distributions.

    Clone-level selection is unconstrained, within-clone state selection is Lam_in-bounded, and transitions keep
    the clone. tau_min = sum_c Q(c) tau_min(P(.|c), Q(.|c); Lam_in) when supp Q(c) is within supp P(c).
    """
    Pc, Qc = np.asarray(Pc, float), np.asarray(Qc, float)
    tot = 0.0
    for c in range(Pc.shape[0]):
        qc = Qc[c].sum()
        if qc <= 0:
            continue
        if Pc[c].sum() <= 0:
            return np.nan  # metastatic clone unseen in primary: not identified
        tot += qc * tau_min(Pc[c], Qc[c], Lam_in)
    return tot


def tau_lower_clone_robust(kPc, kQc, Lam_in, alpha=0.05):
    """Clone-resolved lower confidence bound: per-clone robust bounds at level alpha/(2C).

    These are combined with the least favourable clone weights Q(c) inside a Clopper-Pearson box (level alpha/2).
    """
    kPc, kQc = np.asarray(kPc, float), np.asarray(kQc, float)
    C = kPc.shape[0]
    t = np.zeros(C)
    for c in range(C):
        if kQc[c].sum() == 0:
            t[c] = 0.0
        elif kPc[c].sum() == 0:
            t[c] = 0.0  # unseen primary clone: no constraint
        else:
            t[c] = tau_lower_robust(kPc[c], kQc[c], Lam_in, alpha=alpha / (2 * C))
    lo, hi = cp_box(kQc.sum(1), alpha / 2)
    # minimise sum_c qc t_c over the box with sum qc = 1: greedy, filling the smallest t_c first
    q = lo.copy(); rem = 1 - q.sum()
    for c in np.argsort(t):
        add = min(hi[c] - q[c], rem); q[c] += add; rem -= add
    return float((q * t).sum())
