"""Milestone 1 experiments E1-E7 (protocols/simulation_protocol.md). Writes results/tables, results/metrics.

Usage: python -I experiments/run_all.py [--quick]
"""
import json, os, sys, time
import numpy as np, pandas as pd
from scipy.optimize import least_squares
from scipy.stats import ks_2samp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.bounds.core import (tau_min, tau_min_lp, tau_max_lp, tv, support_mismatch, tau_lower_robust,
                             tau_min_clone, tau_lower_clone_robust, p_star, cp_box)
from src.simulators import scenarios as S

QUICK = "--quick" in sys.argv
R = 60 if QUICK else 400
SEED = 20261123
TAB, MET = os.path.join(ROOT, "results/tables"), os.path.join(ROOT, "results/metrics")
os.makedirs(TAB, exist_ok=True); os.makedirs(MET, exist_ok=True)
metrics = {}


def save(df, name):
    df.to_csv(os.path.join(TAB, name), index=False)


# ---------------- E1: closed form vs LP, interval, timing ----------------
def e1():
    rng = np.random.default_rng(SEED)
    rows = []
    for i in range(1000):
        d = int(rng.integers(3, 13))
        P, Q = rng.dirichlet(np.ones(d) * .5), rng.dirichlet(np.ones(d) * .5)
        if rng.random() < .25:
            P[rng.integers(d)] = 0; P /= P.sum()
        lam = float(np.exp(rng.uniform(0, np.log(1000))))
        a, b = tau_min(P, Q, lam), tau_min_lp(P, Q, lam)
        rows.append(dict(i=i, d=d, lam=lam, closed=a, lp=b, diff=abs(a - b), tv=tv(P, Q), sm=support_mismatch(P, Q),
                         tmax=tau_max_lp(P, Q, lam) if i < 200 else np.nan,
                         t1=tau_min(P, Q, 1.0), tinf=tau_min(P, Q, np.inf), mono_viol=tau_min(P, Q, lam * 2) - a))
    df = pd.DataFrame(rows); save(df, "e1_closed_vs_lp.csv")
    tim = {}
    for d in (100, 1000):
        P, Q = rng.dirichlet(np.ones(d) * .5), rng.dirichlet(np.ones(d) * .5)
        t0 = time.time(); tau_min(P, Q, 10.0); tim[f"closed_d{d}_s"] = time.time() - t0
        kP, kQ = rng.multinomial(50 * d, P), rng.multinomial(50 * d, Q)
        t0 = time.time(); tau_lower_robust(kP, kQ, 10.0); tim[f"robustLP_d{d}_s"] = time.time() - t0
    metrics["E1"] = dict(n=len(df), max_abs_diff=float(df["diff"].max()),
                         lam1_equals_tv=float((df.t1 - df.tv).abs().max()),
                         laminf_equals_sm=float((df.tinf - df.sm).abs().max()),
                         max_monotonicity_violation=float(df.mono_viol.max()),
                         tmax_median=float(df.tmax.median()), tmax_min=float(df.tmax.min()), **tim)


# ---------------- E2: nontriviality (G3) ----------------
def e2():
    rng = np.random.default_rng(SEED + 2)
    rows = []
    for i in range(1000):
        d = int(rng.integers(3, 13))
        P, Q = rng.dirichlet(np.ones(d) * .5), rng.dirichlet(np.ones(d) * .5)
        if rng.random() < .25:
            P[rng.integers(d)] = 0; P /= P.sum()
        T, SM = tv(P, Q), support_mismatch(P, Q)
        for lam in (2, 5, 10, 100):
            t = tau_min(P, Q, lam)
            rows.append(dict(i=i, lam=lam, tau_min=t, tv=T, sm=SM, gap_low=t - SM, gap_high=T - t))
    df = pd.DataFrame(rows); save(df, "e2_nontriviality.csv")
    g = df.assign(both=(df.gap_low >= .05) & (df.gap_high >= .05)).groupby("lam").both.mean()
    metrics["E2"] = {f"frac_strictly_between_lam{k}": float(v) for k, v in g.items()}
    # A second, biologically shaped family: Q is a modest perturbation of P (frozen family = the Dirichlet one above).
    rows = []
    for i in range(1000):
        d = 8; P = rng.dirichlet(np.ones(d) * 2)
        w = np.exp(rng.normal(0, 1.0, d)); K = 0.85 * np.eye(d) + 0.15 * rng.dirichlet(np.ones(d), d)
        Q = (P * w / (P * w).sum()) @ K
        for lam in (2, 5, 10, 100):
            t = tau_min(P, Q, lam)
            G = (P * w / (P * w).sum())[:, None] * K
            rows.append(dict(i=i, lam=lam, tau_min=t, tv=tv(P, Q), sm=support_mismatch(P, Q), tau_true=1 - np.trace(G),
                             lam_true=w.max() / w.min()))
    d2 = pd.DataFrame(rows); save(d2, "e2_nontriviality_perturbation_family.csv")
    metrics["E2_perturbation_family"] = {
        str(k): dict(median_ratio_to_tv=float((v.tau_min / v.tv).median()),
                     median_ratio_to_tau_true=float((v.tau_min / v.tau_true).median()),
                     frac_positive=float((v.tau_min > 1e-9).mean()),
                     frac_lam_true_le_lam=float((v.lam_true <= k).mean())) for k, v in d2.groupby("lam")}


# ---------------- E3: counterexamples, numerically ----------------
def e3():
    out = {}
    P = np.array([0.90, 0.09, 0.01]); Q = np.array([0.40, 0.30, 0.30])
    sel = S.joint(P, Q / P, np.eye(3)); ind = S.joint(P, np.ones(3), np.tile(Q, (3, 1)))
    out["C1"] = dict(P=P.tolist(), Q=Q.tolist(), Q_sel=sel["Q"].tolist(), Q_ind=ind["Q"].tolist(),
                     tau_sel=sel["tau"], tau_ind=ind["tau"], lam_sel=sel["lam"], tv=tv(P, Q),
                     tau_min_curve={str(l): tau_min(P, Q, l) for l in (1, 2, 5, 10, 30, 67.5, 100)},
                     tau_max_lp_lam100=tau_max_lp(P, Q, 100))
    # C2 Simpson: population P = Q but within-clone induction
    Pc = np.array([[0.4, 0.1], [0.1, 0.4]]); Qc = np.array([[0.1, 0.4], [0.4, 0.1]])
    out["C2"] = dict(pop_P=Pc.sum(0).tolist(), pop_Q=Qc.sum(0).tolist(),
                     pop_tau_min_lam1=tau_min(Pc.sum(0), Qc.sum(0), 1.0),
                     clone_tau_min_lam2=tau_min_clone(Pc, Qc, 2.0), clone_tau_min_lam4=tau_min_clone(Pc, Qc, 4.0),
                     clone_tau_min_lam16=tau_min_clone(Pc, Qc, 16.0))
    # C3 clone-pure states: clone data add nothing
    Pc = np.diag([0.7, 0.2, 0.1]); Qc = np.diag([0.1, 0.3, 0.6])
    out["C3"] = dict(clone_tau_min_lam1=tau_min_clone(Pc, Qc, 1.0), pop_tau_min_lam60=tau_min(Pc.sum(0), Qc.sum(0), 60))
    # C4 two snapshots, 2-state
    x0, x1, t = 0.05, 0.30, 1.0
    db = (np.log(x1 / (1 - x1)) - np.log(x0 / (1 - x0))) / t
    a = -np.log((1 - x1) / (1 - x0)) / t
    out["C4"] = dict(x0=x0, x1=x1, pure_selection_delta_b=db, pure_induction_a=a,
                     x1_sel=float(_traj(x0, 0, 0, db, [t])[0]), x1_ind=float(_traj(x0, a, 0, 0, [t])[0]))
    # C5 equilibrium: same stationary fraction from selection-dominated and switching-dominated generators
    out["C5"] = _equilibrium_pair()
    metrics["E3"] = out


def _rhs(x, a, c, db):
    return a * (1 - x) - c * x + db * x * (1 - x)


def _traj(x0, a, c, db, ts, n=4000):
    """RK4 integration of the normalised 2-state replicator-mutator ODE."""
    ts = np.asarray(ts, float); out = []; x = x0; tcur = 0.0
    for tt in ts:
        h = (tt - tcur) / n if tt > tcur else 0
        for _ in range(n if h else 0):
            k1 = _rhs(x, a, c, db); k2 = _rhs(x + h * k1 / 2, a, c, db)
            k3 = _rhs(x + h * k2 / 2, a, c, db); k4 = _rhs(x + h * k3, a, c, db)
            x += h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        tcur = tt; out.append(x)
    return np.array(out)


def _equilibrium_pair():
    # stationary x* solves a - (a + c - db) x - db x^2 = 0. Both parameter sets below give x* = 0.3.
    xs = 0.3
    A = dict(a=0.01, db=0.5); A["c"] = (A["a"] * (1 - xs) + A["db"] * xs * (1 - xs)) / xs
    B = dict(a=0.3, db=0.0); B["c"] = B["a"] * (1 - xs) / xs
    return dict(selection_dominated=A, switching_dominated=B,
                xstar_A=float(_traj(xs, A["a"], A["c"], A["db"], [50])[0]),
                xstar_B=float(_traj(xs, B["a"], B["c"], B["db"], [50])[0]))


# ---------------- E4: scenarios S1-S8 ----------------
def _eval(sc, lam, nP, nQ, rng, eps_meas=0.0, Qobs=None, reps=R):
    Q = sc["Q"] if Qobs is None else Qobs
    rec = []
    for _ in range(reps):
        kP, kQ = S.sample(sc["P"], nP, rng), S.sample(Q, nQ, rng)
        Ph, Qh = kP / nP, kQ / nQ
        rec.append(dict(plug=tau_min(Ph, Qh, lam), robust=tau_lower_robust(kP, kQ, lam, eps_meas=eps_meas),
                        tv=tv(Ph, Qh), sm=support_mismatch(Ph, Qh)))
    return pd.DataFrame(rec)


def _summ(df, tau, name, **kw):
    r = dict(scenario=name, tau_true=tau, **kw)
    for e in ("plug", "robust", "tv", "sm"):
        r[f"{e}_mean"] = df[e].mean(); r[f"{e}_pos"] = (df[e] > 1e-9).mean()
        r[f"{e}_cover"] = (df[e] <= tau + 1e-9).mean()
    return r


def e4():
    rng = np.random.default_rng(SEED + 4)
    rows = []
    scen = dict(S1=(S.s1(), 50.0), S2=(S.s2(), 50.0), S3=(S.s3(), 10.0), S6=(S.s6(), None))
    for n in (200, 1000, 2000, 5000, 20000):
        for name, (sc, lam) in scen.items():
            lam = sc["lam"] if lam is None else lam
            rows.append(_summ(_eval(sc, lam, n, n, rng), sc["tau"], name, n=n, lam=lam, lam_true=sc["lam"],
                              tau_min_pop=tau_min(sc["P"], sc["Q"], lam)))
        for p in (1e-4, 5e-4, 2e-3):
            sc = S.s4(p)
            rows.append(_summ(_eval(sc, sc["lam"], n, n, rng), sc["tau"], f"S4_p{p:g}", n=n, lam=sc["lam"],
                              lam_true=sc["lam"], tau_min_pop=tau_min(sc["P"], sc["Q"], sc["lam"]),
                              p_unseen=(1 - p) ** n))
    # S7: pure selection + organ-specific misclassification
    sc = S.s1()
    for eb in (0.02, 0.05, 0.10):
        Qo = S.misclassify(sc["Q"], eb)
        for allow in (0.0, eb):
            rows.append(_summ(_eval(sc, 50.0, 2000, 2000, rng, eps_meas=allow, Qobs=Qo), 0.0, "S7", n=2000,
                              lam=50.0, eps_b=eb, allowance=allow, tau_min_pop=tau_min(sc["P"], Qo, 50.0)))
    # S8: misspecified Lambda
    for lam in (5.0, 10.0, 20.0, 50.0):
        rows.append(_summ(_eval(sc, lam, 2000, 2000, rng), 0.0, "S8", n=2000, lam=lam, lam_true=50.0,
                          tau_min_pop=tau_min(sc["P"], sc["Q"], lam)))
    df = pd.DataFrame(rows); save(df, "e4_scenarios.csv")
    # equivalence twins: same (P,Q), different tau
    tw = []
    for name, sc, lam in (("S1", S.s1(), 50.0), ("S3", S.s3(), 10.0), ("S4_p5e-4", S.s4(5e-4), None)):
        lam = sc["lam"] if lam is None else lam
        t = S.induction_twin(sc)
        a = _eval(sc, lam, 2000, 2000, np.random.default_rng(1), reps=200)
        b = _eval(t, lam, 2000, 2000, np.random.default_rng(1), reps=200)
        tw.append(dict(scenario=name, tau_true=sc["tau"], tau_twin=t["tau"], maxdiff_PQ=float(
            np.abs(sc["Q"] - t["Q"]).max()), robust_identical=bool(np.allclose(a.robust, b.robust)),
            robust_mean=a.robust.mean(), cover_twin=(b.robust <= t["tau"] + 1e-9).mean()))
    save(pd.DataFrame(tw), "e4_equivalence_twins.csv")
    # S5 clones
    c5 = []
    for eps in (0.0, 0.05, 0.2):
        sc = S.s5()
        for _ in range(R // 2):
            kPc = S.sample_clones(sc["Pc"], 3000, rng, eps); kQc = S.sample_clones(sc["Qc"], 3000, rng, eps)
            c5.append(dict(eps=eps, tau_true=sc["tau"], tau_min_clone_pop=tau_min_clone(sc["Pc"], sc["Qc"], 5.0),
                           clone_robust=tau_lower_clone_robust(kPc, kQc, 5.0),
                           pop_robust_lam5=tau_lower_robust(kPc.sum(0), kQc.sum(0), 5.0),
                           pop_robust_lam5000=tau_lower_robust(kPc.sum(0), kQc.sum(0), 5000.0)))
    c5 = pd.DataFrame(c5); save(c5, "e4_s5_clones.csv")
    # S5 null twin: same clones, no induction -> clone bound should be ~0, population lam5 bound misattributes
    sc0 = S.s5(h=0.0); c50 = []
    for _ in range(R // 2):
        kPc = S.sample_clones(sc0["Pc"], 3000, rng); kQc = S.sample_clones(sc0["Qc"], 3000, rng)
        c50.append(dict(clone_robust=tau_lower_clone_robust(kPc, kQc, 5.0),
                        pop_robust_lam5=tau_lower_robust(kPc.sum(0), kQc.sum(0), 5.0)))
    c50 = pd.DataFrame(c50); save(c50, "e4_s5_null.csv")
    metrics["E4_S5"] = dict(tau_true=float(sc["tau"]), clone_tau_min_population=float(c5.tau_min_clone_pop.iloc[0]),
                            **{f"eps{e}_{k}": float(v) for e, g in c5.groupby("eps")
                               for k, v in dict(clone_cover=(g.clone_robust <= g.tau_true + 1e-9).mean(),
                                                clone_mean=g.clone_robust.mean(),
                                                pop5_mean=g.pop_robust_lam5.mean(),
                                                pop5000_mean=g.pop_robust_lam5000.mean()).items()},
                            null_clone_pos=float((c50.clone_robust > 1e-9).mean()),
                            null_pop_lam5_pos=float((c50.pop_robust_lam5 > 1e-9).mean()),
                            s5_null_tau=float(sc0["tau"]))


# ---------------- E5: design calculus ----------------
def e5():
    rows = []
    for q in (0.05, 0.1, 0.3):
        for lam in (2, 10, 100, 1000, 10000):
            for rho in (1, 10):
                ps = p_star(q, lam * rho)
                # smallest n with CP upper bound for k = 0 (level alpha/(2d) per state, d = 6) below p*
                a = 0.05 / 2 / 6
                n_cp = int(np.ceil(np.log(a / 2) / np.log(1 - ps)))
                rows.append(dict(q=q, lam=lam, rho=rho, p_star=ps, n_required_cells=n_cp,
                                 n_lineage_recording=int(np.ceil(np.log(0.05) / np.log(1 - q)))))
    df = pd.DataFrame(rows); save(df, "e5_design_required_n.csv")
    # simulation check of the analytic law: pure induction into absent state, q = 0.3
    rng = np.random.default_rng(SEED + 5); chk = []
    sc = S.s2(0.3)
    for lam in (10, 100, 1000):
        for n in (100, 300, 1000, 3000, 10000, 30000, 100000):
            pos = np.mean([tau_lower_robust(S.sample(sc["P"], n, rng), S.sample(sc["Q"], 5000, rng), lam) > 1e-9
                           for _ in range(60 if not QUICK else 20)])
            chk.append(dict(lam=lam, nP=n, power=pos))
    save(pd.DataFrame(chk), "e5_design_power_check.csv")
    # regional sampling bias: the primary sample comes from a region where state s is depleted rho-fold.
    # True Lambda = 50 (S1, pure selection) but the observed primary looks as if p were rho times smaller.
    rb = []
    for rho in (1, 3, 10):
        sc = S.s1(); Pobs = sc["P"].copy(); Pobs[5] /= rho; Pobs /= Pobs.sum()
        for n in (2000, 20000):
            f = np.mean([tau_lower_robust(S.sample(Pobs, n, rng), S.sample(sc["Q"], n, rng), 50.0) > 1e-9
                         for _ in range(100 if not QUICK else 20)])
            rb.append(dict(rho=rho, n=n, false_induction=f))
    save(pd.DataFrame(rb), "e5_regional_bias.csv")


# ---------------- E6: multi-metastasis / founder bottleneck ----------------
def _met(pi_f, nf, rng, mode, h=0.3, ncell=2000):
    """One metastasis: nf founders drawn from pi_f; Yule clone sizes ~ Exp(1).

    mode 'heritable' is pure selection; 'founder_ind' converts each founder to state s with prob h (heritable);
    'cell_ind' converts each descendant cell independently with prob h.
    """
    d = len(pi_f)
    fs = rng.choice(d, nf, p=pi_f)
    if mode == "founder_ind":
        fs = np.where(rng.random(nf) < h, d - 1, fs)
    z = rng.exponential(1.0, nf); comp = np.zeros(d)
    np.add.at(comp, fs, z); comp /= comp.sum()
    if mode == "cell_ind":
        comp = (1 - h) * comp; comp[-1] += h
    return rng.multinomial(ncell, comp)


def e6():
    rng = np.random.default_rng(SEED + 6)
    P_sel = S.s1(p=0.005, ratio=50)["P"]
    pi_sel = S.s1(p=0.005, ratio=50)["pi"]       # selection: s pre-exists, founder prob ~0.2
    h_match = pi_sel[5]                           # induction matched to the same founder-level s probability
    pi_ind_base = pi_sel.copy(); pi_ind_base[5] = 0; pi_ind_base /= pi_ind_base.sum()
    rows = []
    M = 8
    for nf in (1, 3, 10, 100):
        for rep in range(R // 4):
            fr = {}
            for mode in ("heritable", "founder_ind", "cell_ind"):
                if mode == "heritable":
                    mets = [_met(pi_sel, nf, rng, mode) for _ in range(M)]
                else:
                    mets = [_met(pi_ind_base, nf, rng, mode, h=h_match) for _ in range(M)]
                f = np.array([m[5] / m.sum() for m in mets]); fr[mode] = f
                rows.append(dict(nf=nf, rep=rep, mode=mode, mean_s=f.mean(), var_s=f.var(ddof=1)))
    df = pd.DataFrame(rows); save(df, "e6_fluctuation.csv")
    eq = []
    for nf, g in df.groupby("nf"):
        a, b, c = (g[g["mode"] == m].var_s.values for m in ("heritable", "founder_ind", "cell_ind"))
        eq.append(dict(nf=nf, ks_p_sel_vs_founderind=ks_2samp(a, b).pvalue, ks_p_sel_vs_cellind=ks_2samp(a, c).pvalue,
                       var_sel=np.median(a), var_founder=np.median(b), var_cell=np.median(c)))
    save(pd.DataFrame(eq), "e6_fluctuation_tests.csv")
    # single metastasis plug-in bound under neutral, transition-free seeding (tau = 0, Lambda = 1)
    sm = []
    for nf in (1, 3, 10, 100, 1000):
        for lam in (1.0, 10.0, 100.0):
            pos = []
            for _ in range(R // 4):
                kQ = _met(P_sel, nf, rng, "heritable"); kP = S.sample(P_sel, 5000, rng)
                pos.append(tau_lower_robust(kP, kQ, lam) > 1e-9)
            sm.append(dict(nf=nf, lam=lam, false_induction=np.mean(pos)))
    save(pd.DataFrame(sm), "e6_single_met_neutral.csv")


# ---------------- E7: dynamics identifiability ----------------
def e7():
    rng = np.random.default_rng(SEED + 7)
    truth = dict(a=0.08, c=0.02, db=0.6); x0 = 0.02
    rows = []
    for design, ts in (("2 snapshots", [1.5]), ("4 snapshots", [0.5, 1.5, 3, 5]), ("8 snapshots", np.linspace(.5, 6, 8))):
        xt = _traj(x0, truth["a"], truth["c"], truth["db"], ts, n=400)
        for ncell in (2000, 20000, None):
            for rep in range(30 if ncell else 1):
                obs = xt if ncell is None else rng.binomial(ncell, xt) / ncell
                x0o = x0 if ncell is None else rng.binomial(ncell, x0) / ncell
                best = None
                for init in ((0.1, 0.1, 0.1), (0.5, 0.01, 0.01), (0.01, 0.01, 1.0)):
                    f = least_squares(lambda th: _traj(x0o, *np.abs(th), ts, n=200) - obs, init, bounds=(0, 5))
                    if best is None or f.cost < best.cost:
                        best = f
                a, c, db = np.abs(best.x)
                rows.append(dict(design=design, ncell=ncell or 0, rep=rep, a=a, c=c, db=db, cost=best.cost))
    df = pd.DataFrame(rows); save(df, "e7_dynamics_fits.csv")


if __name__ == "__main__":
    for f in (e1, e2, e3, e4, e5, e6, e7):
        t0 = time.time(); f(); print(f.__name__, f"{time.time() - t0:.1f}s", flush=True)
    with open(os.path.join(MET, "milestone1_metrics.json"), "w") as fh:
        json.dump(metrics, fh, indent=1, default=float)
