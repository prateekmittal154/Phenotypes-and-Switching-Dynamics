"""Ground-truth scenarios S1-S8 (protocols/simulation_protocol.md). Each returns the exact law plus truth."""
import numpy as np
from scipy.linalg import expm

D = 6
BASE = np.array([0.40, 0.25, 0.15, 0.10, 0.10, 0.0])  # state 5 = the 'metastatic' state s


def _norm(v):
    v = np.asarray(v, float)
    return v / v.sum()


def joint(P, w, K):
    """Exact founder x metastasis law gamma, Q, tau and Lambda_true (on supp P)."""
    pi = _norm(P * w)
    G = pi[:, None] * K
    sp = P > 0
    return dict(P=P, Q=G.sum(0), pi=pi, tau=1 - np.trace(G), lam=w[sp].max() / w[sp].min(), G=G)


def with_rare(p):
    P = BASE.copy(); P[5] = p
    P[:5] *= (1 - p) / P[:5].sum()
    return P


def s1(p=0.005, ratio=50.0):
    P = with_rare(p); w = np.ones(D); w[5] = ratio
    return joint(P, w, np.eye(D))


def s2(h=0.3):
    P = with_rare(0.0); K = (1 - h) * np.eye(D); K[:, 5] += h
    return joint(P, np.ones(D), K)


def s3(ratio=10.0, h=0.3):
    P = with_rare(0.0); P[4] = 0.01; P = _norm(P)
    w = np.ones(D); w[4] = ratio
    K = np.eye(D); K[0, 0] = 1 - h; K[0, 1] = h
    return joint(P, w, K)


def s4(p, q=0.3):
    ratio = q * (1 - p) / ((1 - q) * p)
    return s1(p=p, ratio=ratio)


def induction_twin(sc):
    """Same (P, Q), with w = 1 and K(y|x) = Q(y)."""
    P, Q = sc["P"], sc["Q"]
    return joint(P, np.ones(len(P)), np.tile(Q, (len(P), 1)))


def s5(ratio_in=5.0, clone_sel=(1.0, 1000.0, 30.0), h=0.2, seed=5):
    """3 clones; within-clone selection ratio <= ratio_in; unconstrained clone selection.

    Clone 1 has induction h into its rarest primary state (deviation D1: parameters fixed after the quick run).
    """
    rng = np.random.default_rng(seed)
    C = 3
    Pc_cond = rng.dirichlet(np.ones(D) * 1.0, size=C)
    Pc = np.array([0.6, 0.3, 0.1])[:, None] * Pc_cond
    Gs, Qc = [], np.zeros((C, D))
    tau = 0.0
    r_all = []
    for c in range(C):
        r = rng.uniform(1, ratio_in, D); r[0] = 1; r[1] = ratio_in  # attains the bound exactly
        r_all.append(r)
        w = clone_sel[c] * r
        K = np.eye(D)
        if c == 1:
            tgt = int(np.argmin(Pc_cond[1]))
            K = (1 - h) * np.eye(D); K[:, tgt] += h  # induction into the rarest state of clone 1
        Gs.append((Pc[c] * w)[:, None] * K)
    Z = sum(G.sum() for G in Gs)
    for c in range(C):
        G = Gs[c] / Z
        Qc[c] = G.sum(0); tau += G.sum() - np.trace(G)
    return dict(Pc=Pc, Qc=Qc, tau=tau, P=Pc.sum(0), Q=Qc.sum(0))


def s6(t=8.0, a=0.05, c=0.10, b=(1.0, 1.3), w_dis=(1.0, 5.0), P=(0.95, 0.05)):
    """Reversible 2-state switching during expansion. Joint (founder, descendant) = pi(x) E(x,y), E = exp(t(B+T))."""
    P = np.array(P); T = np.array([[-a, a], [c, -c]]); B = np.diag(b)
    E = expm(t * (B + T))
    w = np.array(w_dis) * E.sum(1)          # effective founder weight = dissemination x expected output
    K = E / E.sum(1, keepdims=True)
    return joint(P, w, K)


def misclassify(Q, eps):
    d = len(Q); M = (1 - eps) * np.eye(d) + eps * (np.ones((d, d)) - np.eye(d)) / (d - 1)
    return Q @ M


def sample(dist, n, rng):
    return rng.multinomial(n, _norm(dist))


def sample_clones(Pc, n, rng, eps=0.0):
    """Sample n cells from a clone x state law; each clone label is replaced at random with probability eps."""
    C, d = Pc.shape
    k = rng.multinomial(n, _norm(Pc.ravel())).reshape(C, d)
    if eps <= 0:
        return k
    out = np.zeros_like(k)
    for c in range(C):
        for x in range(d):
            if k[c, x] == 0:
                continue
            flip = rng.binomial(k[c, x], eps); out[c, x] += k[c, x] - flip
            if flip:
                others = [o for o in range(C) if o != c]
                out[others, x] += rng.multinomial(flip, np.ones(C - 1) / (C - 1))
    return out
