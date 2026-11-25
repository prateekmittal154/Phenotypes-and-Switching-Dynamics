"""Milestone 4 figures. Usage: python figures_m4.py <tables_dir> <out_fig_dir>"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
MK = ["o", "s", "^", "D", "v"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e6", "grid.linewidth": .6, "lines.linewidth": 2, "lines.markersize": 6})
T, F = sys.argv[1], sys.argv[2]
os.makedirs(F, exist_ok=True)
have = lambda p: os.path.exists(os.path.join(T, p))
rd = lambda p: pd.read_csv(os.path.join(T, p))


def save(fig, n):
    fig.tight_layout(); fig.savefig(os.path.join(F, n + ".png"), dpi=200); fig.savefig(os.path.join(F, n + ".pdf"))
    plt.close(fig)


# Fig 22: audited Fisher information — optimum location and flatness
if have("e12_fisher_full.csv"):
    d = rd("e12_fisher_full.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.3))
    g = d[(d.k == 0.5)]
    for j, pi in enumerate(sorted(g.pi.unique())):
        s = g[g.pi == pi].sort_values("kL")
        axs[0].plot(s.kL, s.I_logk_pi_joint / s.I_logk_pi_joint.max(), color=C[j % 5], label=f"π = {pi}")
    axs[0].axvspan(0.10, 1.0, color="#999", alpha=.15)
    axs[0].axvline(0.43, color="#333", ls="--", lw=1.2)
    axs[0].annotate("k·L* ≈ 0.43", xy=(0.47, .45), fontsize=7.5)
    axs[0].annotate("≥50% of peak\nover a 10-fold band", xy=(1.15, .62), fontsize=7)
    axs[0].set_xscale("log"); axs[0].set_xlabel("k · L"); axs[0].set_ylabel("information, scaled to peak")
    axs[0].set_title("Optimum depends mildly on occupancy, and is flat", fontsize=9); axs[0].legend(fontsize=7)
    if have("e12_optimum_audit.csv"):
        o = rd("e12_optimum_audit.csv")
        axs[1].scatter(o.kL_opt_pi_known, o.kL_opt_pi_joint, c=[C[0]] * len(o), s=28)
        lim = [0.36, 0.52]
        axs[1].plot(lim, lim, color="#888", ls="--", lw=1)
        axs[1].set_xlim(lim); axs[1].set_ylim(lim)
        axs[1].set_xlabel("k·L* with occupancy known"); axs[1].set_ylabel("k·L* with occupancy estimated")
        axs[1].set_title("Estimating occupancy costs ≤6% of the information", fontsize=9)
    save(fig, "m4_fig22_fisher_audit")

# Fig 23: full-tree vs naive pair-sum information
if have("e13_tree_vs_pair_mc.csv"):
    d = rd("e13_tree_vs_pair_mc.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.3))
    for j, k in enumerate(sorted(d.k.unique())[:5]):
        g = d[(d.k == k) & (d.pi == 0.5)].sort_values("depth")
        if len(g):
            axs[0].plot(g.depth, g.inflation, color=C[j % 5], marker=MK[j % 5], label=f"k = {k}")
    axs[0].axhline(1, color="#333", lw=1)
    axs[0].set_yscale("log"); axs[0].set_xlabel("tree depth"); axs[0].set_ylabel("naive pair-sum / full-tree")
    axs[0].set_title("Treating pairs as independent overstates information", fontsize=9)
    axs[0].legend(fontsize=7)
    v = d.dropna(subset=["I_exact_logk"])
    if len(v):
        axs[1].scatter(v.I_exact_logk, v.I_tree_logk, color=C[0], s=28)
        lim = [min(v.I_exact_logk.min(), v.I_tree_logk.min()), max(v.I_exact_logk.max(), v.I_tree_logk.max())]
        axs[1].plot(lim, lim, color="#888", ls="--", lw=1)
        axs[1].set_xlabel("exact (full enumeration)"); axs[1].set_ylabel("Monte Carlo estimate")
        axs[1].set_title("Monte Carlo validated against exact enumeration", fontsize=9)
    save(fig, "m4_fig23_full_tree_information")

# Fig 24: calibration
if have("e14_calibration_summary.csv"):
    s = rd("e14_calibration_summary.csv")
    piv = s.pivot(index="scen", columns="test", values="rate")
    lo = s.pivot(index="scen", columns="test", values="lo"); hi = s.pivot(index="scen", columns="test", values="hi")
    order = [x for x in ["H0_pi", "H0_k", "H0_k_frozen", "MIS", "H1_pi", "H1_k"] if x in piv.index]
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.4))
    x = np.arange(len(order)); w = .38
    for j, (t, lab) in enumerate([("chi2_pi", "shared-occupancy test"), ("chi2_k", "shared-rate test")]):
        if t not in piv:
            continue
        v = piv.loc[order, t].values
        err = [v - lo.loc[order, t].values, hi.loc[order, t].values - v]
        axs[0].bar(x + (j - .5) * w, v, w, color=C[j], yerr=err, capsize=2, ecolor="#444",
                   hatch=["", "//"][j], edgecolor="white", label=lab)
    axs[0].axhline(.05, color="#333", ls=":", lw=1.2)
    axs[0].set_xticks(x, order, rotation=18, ha="right"); axs[0].set_ylabel("rejection rate")
    axs[0].set_title("χ² reference: rejection rates (95% Wilson)", fontsize=9); axs[0].legend(fontsize=7)
    if "cover_V_pi" in piv:
        v = piv.loc[order, "cover_V_pi"].values
        err = [v - lo.loc[order, "cover_V_pi"].values, hi.loc[order, "cover_V_pi"].values - v]
        axs[1].bar(x, v, .55, color=C[2], yerr=err, capsize=2, ecolor="#444", edgecolor="white")
        axs[1].axhline(.95, color="#333", ls=":", lw=1.2)
        axs[1].set_xticks(x, order, rotation=18, ha="right")
        axs[1].set_ylabel("coverage of the V$_\\pi$ interval"); axs[1].set_ylim(0, 1.05)
        axs[1].set_title("Interval coverage collapses in the frozen regime", fontsize=9)
    save(fig, "m4_fig24_calibration")

# Fig 25: recorder capacity design law
if have("e15b_capacity_designs.csv"):
    d = rd("e15b_capacity_designs.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.4))
    g = d.groupby(["S", "D"]).decades_covered.first().reset_index()
    for j, D in enumerate(sorted(g.D.unique())):
        s = g[g.D == D].sort_values("S")
        axs[0].plot(s.S, s.decades_covered, color=C[j % 5], marker=MK[j % 5], label=f"{int(D)} divisions")
    axs[0].set_xscale("log", base=2); axs[0].set_xlabel("recording sites S")
    axs[0].set_ylabel("decades of switching rate covered")
    axs[0].set_title("Capacity law: W = 1 + log$_{10}$(resolved levels)", fontsize=9); axs[0].legend(fontsize=7)
    if have("e16_capacity_vs_prediction.csv"):
        e = rd("e16_capacity_vs_prediction.csv").sort_values("k_lo")
        y = np.arange(len(e))
        for i, (_, row) in enumerate(e.iterrows()):
            # idealised law band
            lo_i = 0.10 / max(row.depth_max, 1e-9); hi_i = 1.0
            axs[1].plot([lo_i, hi_i], [i + .17, i + .17], lw=5, color="#bbb", solid_capstyle="butt",
                        label="idealised law" if i == 0 else None)
            axs[1].plot([row.k_lo, row.k_hi], [i - .17, i - .17], lw=5,
                        color=C[0] if row.system.startswith("macs") else C[1], solid_capstyle="butt",
                        label=("macsGESTALT" if row.system.startswith("macs") else "Cassiopeia")
                        if row.system[:4] not in [x[:4] for x in e.system[:i]] else None)
            axs[1].text(row.k_hi * 1.15, i - .17, f"{row.decades:.2f} dec", va="center", fontsize=6.5)
        axs[1].set_yticks(y, [f"{a}\n{b}" for a, b in zip(e.system, e.dataset)], fontsize=6.5)
        axs[1].set_xscale("log"); axs[1].set_xlim(0.002, 40)
        axs[1].set_xlabel("switching rate k covered (per resolved branch)")
        axs[1].axvspan(1.0, 40, color="#999", alpha=.13)
        axs[1].annotate("rates fitted\nin real trees", xy=(2.2, 2.2), fontsize=6.5)
        axs[1].legend(fontsize=6, loc="lower right")
        axs[1].set_title("Grey: idealised law. Colour: realised geometry", fontsize=9)
    save(fig, "m4_fig25_capacity_law")

# Fig 26: Quinn independent system
if have("m4_quinn_variance_components.csv"):
    d = rd("m4_quinn_variance_components.csv")
    fig, axs = plt.subplots(1, 3, figsize=(12.6, 3.4))
    x = np.arange(len(d))
    lab = [f"{m}\n({int(n)} cl.)" for m, n in zip(d.mouse, d.n_clones)]
    axs[0].bar(x, d.V_pi, .55, color=C[0], edgecolor="white")
    for xi, p1 in enumerate(d.lrt_pi_p):
        axs[0].text(xi, d.V_pi.iloc[xi] + .0015, "*" if p1 < .05 else "ns", ha="center", fontsize=7)
    axs[0].set_xticks(x, lab, fontsize=7.5); axs[0].set_ylabel("V$_\\pi$ (occupancy)")
    axs[0].set_title("Occupancy: heterogeneous in all four mice", fontsize=9)
    pr = rd("m4_quinn_primary_only.csv") if have("m4_quinn_primary_only.csv") else None
    if pr is not None:
        mp = dict(zip(pr.mouse, pr.V_pi))
        axs[0].plot(x, [mp.get(m, np.nan) for m in d.mouse], "k_", ms=16, mew=2,
                    label="primary-site only,\ndepth-adjusted")
        axs[0].legend(fontsize=6.5)
    axs[1].bar(x, d.V_logk, .55, color=C[1], hatch="//", edgecolor="white")
    for xi, p2 in enumerate(d.lrt_k_p):
        axs[1].text(xi, d.V_logk.iloc[xi] + .02, "*" if p2 < .05 else "ns", ha="center", fontsize=7)
    axs[1].set_xticks(x, lab, fontsize=7.5); axs[1].set_ylabel("V$_{\\log k}$ (switching)")
    axs[1].set_title("Switching: significant in 2 of 4 before adjustment", fontsize=9)
    if have("m4_quinn_technical_control.csv"):
        t = rd("m4_quinn_technical_control.csv")
        p = t.pivot_table(index=["mouse", "site"], columns="feature", values="within_clone_shift_sd")
        cols = [c for c in ["cc_score", "n_genes", "lib", "mt_frac"] if c in p.columns]
        y = np.arange(len(p))
        axs[2].axvline(0, color="#444", lw=1)
        for j, c in enumerate(cols):
            axs[2].scatter(p[c], y, color=C[j % 5], marker=MK[j % 5], label=c, zorder=3)
        axs[2].set_yticks(y, [f"{a} {b}" for a, b in p.index], fontsize=6.5)
        axs[2].set_xlabel("within-clone shift (primary SD)"); axs[2].legend(fontsize=6.5)
        axs[2].set_title("Within-clone shifts track sequencing depth", fontsize=9)
    save(fig, "m4_fig26_quinn")

print("written:", sorted(x for x in os.listdir(F) if x.startswith("m4_") and x.endswith(".png")))
