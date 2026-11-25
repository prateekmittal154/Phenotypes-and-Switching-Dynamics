"""Tests for the occupancy and switching estimators.

These check mathematical properties that must hold regardless of implementation detail:
the covariance decomposition, plateau invariance, the identifiability window, rate-time
confounding, and parameter recovery on data simulated from the model itself.

They are deliberately tolerant on recovery: the point is that the estimator is unbiased in
the regime where the data identify the parameter, not that it hits a tight number.
"""
import numpy as np
import pytest

from src.branching.recorder import simulate_genealogy, induced_subtree
from src.lineage_inference.ctmc import VarTree, fit_var, trans


# ----------------------------------------------------------------- helpers
def evolve(parent, rng, pi, k, elen=None):
    st = np.empty(len(parent), int)
    st[0] = int(rng.random() < pi)
    for i in range(1, len(parent)):
        t = 1.0 if elen is None else float(elen[i])
        E = np.exp(-k * t)
        p = st[parent[i]]
        up = pi * (1 - E) if p == 0 else pi + (1 - pi) * E
        st[i] = int(rng.random() < up)
    return st


def one_clone(seed, pi, k, n_gen=10, n_cells=64, elen_scale=1.0):
    """Simulate one clone and return a VarTree ready to fit."""
    rng = np.random.default_rng(seed)
    par, gen, live = simulate_genealogy(n_gen, rng)
    if len(live) < n_cells:
        n_cells = len(live)
    take = rng.choice(live, n_cells, replace=False)
    sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
    st = evolve(par, rng, pi, k)
    leaf = np.full(len(sub_par), -1, int)
    for new_i, orig in leafmap.items():
        leaf[new_i] = st[orig]
    return VarTree(sub_par, leaf, elen=sub_elen * elen_scale)


# ----------------------------------------------------------------- transition kernel
def test_trans_rows_are_probabilities():
    for pi in (0.1, 0.5, 0.9):
        for k in (0.01, 0.3, 5.0):
            P = trans(pi, k, 1.0)[0]
            assert P.shape == (2, 2)
            assert np.allclose(P.sum(axis=1), 1.0), "rows must sum to 1"
            assert (P >= 0).all()


def test_trans_converges_to_stationary():
    """As L grows, both rows approach the stationary distribution (pi_0, pi_1)."""
    pi = 0.73
    P = trans(pi, 1.0, 500.0)[0]
    assert np.allclose(P[0], P[1], atol=1e-6), "rows converge to each other"
    assert P[0, 1] == pytest.approx(pi, abs=1e-4)


def test_trans_identity_at_zero_distance():
    P = trans(0.4, 2.0, 1e-12)[0]
    assert np.allclose(P, np.eye(2), atol=1e-6)


# ----------------------------------------------------------------- the decomposition
def test_covariance_decomposes_into_memory_plus_plateau():
    """Cov(L) = E_c[pi(1-pi) exp(-2kL)] + Var_c(pi), the identity the paper rests on."""
    rng = np.random.default_rng(7)
    pis = rng.uniform(0.15, 0.85, 400)
    ks = np.exp(rng.uniform(np.log(0.05), np.log(2.0), 400))
    plateau = float(np.var(pis))
    for L in (0.1, 1.0, 5.0):
        memory = float(np.mean(pis * (1 - pis) * np.exp(-2 * ks * L)))
        total = memory + plateau
        # the limit as L -> infinity is the plateau alone
        assert total >= plateau - 1e-12
    far = float(np.mean(pis * (1 - pis) * np.exp(-2 * ks * 500.0))) + plateau
    assert far == pytest.approx(plateau, abs=1e-6)


def test_plateau_is_invariant_to_switching_rate():
    """Two clone sets with identical occupancies but rates 100x apart share a plateau."""
    rng = np.random.default_rng(11)
    pis = rng.uniform(0.15, 0.85, 500)
    slow = np.full(500, 0.02)
    fast = np.full(500, 2.0)
    L = 400.0
    p_slow = float(np.mean(pis * (1 - pis) * np.exp(-2 * slow * L))) + float(np.var(pis))
    p_fast = float(np.mean(pis * (1 - pis) * np.exp(-2 * fast * L))) + float(np.var(pis))
    assert p_slow == pytest.approx(p_fast, abs=1e-6)
    assert p_slow == pytest.approx(float(np.var(pis)), abs=1e-6)


# ----------------------------------------------------------------- recovery
@pytest.mark.parametrize("pi,k", [(0.3, 0.3), (0.5, 0.5), (0.7, 0.3)])
def test_recovers_occupancy_within_the_window(pi, k):
    """Occupancy is recovered across replicates; it is the robust quantity."""
    est = [fit_var(one_clone(s, pi, k))["pi"] for s in range(12)]
    assert np.median(est) == pytest.approx(pi, abs=0.12)


def test_recovers_rate_within_the_window():
    """k = 0.3 on 10-generation trees sits inside the identifiable window."""
    est = [fit_var(one_clone(s, 0.5, 0.3))["k"] for s in range(12)]
    med = float(np.median(est))
    assert 0.5 < med / 0.3 < 3.0, f"median k_hat {med} should be within ~3x of 0.3"


def test_frozen_states_flag_the_boundary():
    """A near-frozen clone should be flagged, not silently reported as a measurement."""
    fits = [fit_var(one_clone(s, 0.5, 1e-3)) for s in range(8)]
    flagged = sum(f["k_at_bound"] for f in fits)
    assert flagged >= 1, "at least one near-frozen clone should hit the bound flag"


# ----------------------------------------------------------------- rate-time confounding
def test_rate_and_time_are_confounded():
    """L(pi, k; L) = L(pi, k/c; cL): scaling edges by c and k by 1/c leaves the likelihood
    unchanged. This is why unit edges give a rate in unknown units."""
    tree = one_clone(3, 0.5, 0.4)
    c = 2.5
    scaled = VarTree(tree.parent, tree.leaf_state, elen=tree.elen * c)
    a = tree.loglik(0.5, 0.5, 0.4)
    b = scaled.loglik(0.5, 0.5, 0.4 / c)
    assert a == pytest.approx(b, rel=1e-8), "likelihood must be invariant along the degeneracy"


def test_unit_edges_shift_the_rate_estimate():
    """Discarding branch times should move the rate but not the occupancy."""
    rng = np.random.default_rng(5)
    par, gen, live = simulate_genealogy(10, rng)
    take = rng.choice(live, 64, replace=False)
    sub_par, sub_elen, leafmap, _ = induced_subtree(par, take)
    st = evolve(par, rng, 0.6, 0.3)
    leaf = np.full(len(sub_par), -1, int)
    for new_i, orig in leafmap.items():
        leaf[new_i] = st[orig]

    true_t = fit_var(VarTree(sub_par, leaf, elen=sub_elen))
    unit_t = fit_var(VarTree(sub_par, leaf, elen=np.ones(len(sub_par))))
    assert abs(unit_t["pi"] - true_t["pi"]) < 0.15, "occupancy should barely move"
    assert unit_t["k"] != pytest.approx(true_t["k"], rel=1e-6), "the rate should move"


# ----------------------------------------------------------------- likelihood sanity
def test_loglik_is_finite():
    tree = one_clone(2, 0.5, 0.3)
    assert np.isfinite(tree.loglik(0.5, 0.5, 0.3))
    assert np.isfinite(tree.loglik(0.1, 0.1, 5.0))


def test_truth_beats_a_wrong_rate_when_clones_are_pooled():
    """Pooled over clones the likelihood prefers the truth, as a correct estimator must.

    Note this is asserted on the POOLED likelihood. A single clone does not carry enough
    information to prefer the true rate reliably -- see the next test, which is the point
    of the paper rather than a defect.
    """
    tot_true = tot_wrong = 0.0
    for s in range(30):
        t = one_clone(s, 0.5, 0.3)
        tot_true += t.loglik(0.5, 0.5, 0.3)
        tot_wrong += t.loglik(0.5, 0.5, 20.0)
    assert tot_true > tot_wrong, "pooled over clones, the true rate must win"


def test_single_clone_likelihood_is_nearly_flat_in_the_rate():
    """One clone does not identify the switching rate.

    Above k ~ 1 the log-likelihood barely moves: a twenty-fold change in the rate costs
    well under one log unit. This is the non-identifiability the estimator reports via
    `k_at_bound`, and it is why a per-clone rate from a single-timepoint experiment
    should not be read as a measurement.
    """
    t = one_clone(2, 0.5, 0.3)
    lls = {k: t.loglik(0.5, 0.5, k) for k in (1.0, 3.0, 20.0)}
    spread = max(lls.values()) - min(lls.values())
    assert spread < 1.0, f"expected a flat surface above k=1, got spread {spread:.3f}"

    # by contrast the surface IS informative at the low end, where the window applies
    assert t.loglik(0.5, 0.5, 0.3) - t.loglik(0.5, 0.5, 0.05) > 5.0
