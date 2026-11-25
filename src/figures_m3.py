"""Milestone 3 figures. Usage: python figures_m3.py <tables_dir> <out_fig_dir>"""
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


# Fig 15: the audit — the plateau is blind to switching rate
expit = lambda x: 1 / (1 + np.exp(-x)); logit = lambda p: np.log(p / (1 - p))
fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3))
L = np.linspace(0, 6, 300)
pi = np.array([0.2, 0.8]); B = pi.var(); A0 = (pi * (1 - pi)).mean()
for j, k in enumerate([0.02, 0.3, 3.0]):
    axs[0].plot(L, A0 * np.exp(-2 * k * L) + B, color=C[j], marker="", ls=["-", "--", ":"][j],
                label=f"k = {k}")
axs[0].axhline(B, color="#333", lw=1.2)
axs[0].annotate(f"plateau = Var$_c(\\pi_c)$ = {B:.3f}\nidentical for every k",
                xy=(3.2, B + .006), fontsize=7.5)
axs[0].set_xlabel("lineage distance L"); axs[0].set_ylabel("within-clone state covariance")
axs[0].set_title("U1: the plateau measures occupancy, not switching", fontsize=9); axs[0].legend(fontsize=7)
sh = np.linspace(0, 3, 100)
v = [expit(logit(pi) + s).var() for s in sh]
v0 = [expit(logit(np.array([0.5, 0.5])) + s).var() for s in sh]
axs[1].plot(sh, v, color=C[0], label="clones differ ($\\pi_c$ = 0.2, 0.8)")
axs[1].plot(sh, v0, color=C[2], ls="--", label="clones identical ($\\pi_c$ = 0.5)")
axs[1].set_xlabel("common induction shift (logit units)"); axs[1].set_ylabel("plateau Var$_c(\\pi_c)$")
axs[1].set_title("U3: induction rescales a plateau but cannot create one", fontsize=9); axs[1].legend(fontsize=7)
save(fig, "m3_fig15_audit")

# Fig 16: Fisher information and the optimal depth
if have("e11_fisher_information.csv"):
    d = rd("e11_fisher_information.csv"); d = d[d.pi == 0.5]
    fig, axs = plt.subplots(1, 2, figsize=(8.8, 3.3))
    for j, k in enumerate(sorted(d.k.unique())):
        g = d[d.k == k].sort_values("L")
        axs[0].plot(g.L, g.I_logk, color=C[j], label=f"k = {k}")
        axs[1].plot(g.kL, g.I_logk / g.I_logk.max(), color=C[j], marker="", label=f"k = {k}")
    axs[0].set_xscale("log"); axs[0].set_xlabel("lineage distance L between the pair")
    axs[0].set_ylabel("Fisher information about log k"); axs[0].legend(fontsize=7)
    axs[0].set_title("Information peaks at an interior distance", fontsize=9)
    axs[1].axvline(0.405, color="#333", ls="--", lw=1.2)
    axs[1].annotate("k·L* ≈ 0.405", xy=(0.46, .5), fontsize=8)
    axs[1].set_xscale("log"); axs[1].set_xlabel("k · L (lineage distance in switching times)")
    axs[1].set_ylabel("information, scaled to its maximum")
    axs[1].set_title("The optimum collapses onto one constant", fontsize=9)
    save(fig, "m3_fig16_optimal_depth")

# Fig 17: identifiability windows of real recorders
if have("m3_depth_spectra_summary.csv"):
    d = rd("m3_depth_spectra_summary.csv")
    g = d.groupby(["system", "dataset"])[["k_lo_identifiable", "k_hi_identifiable", "depth_max"]].median()
    g = g.sort_values("k_lo_identifiable")
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    y = np.arange(len(g))
    for i, ((sysname, ds), r) in enumerate(g.iterrows()):
        col = C[0] if sysname == "macsGESTALT" else C[1]
        ax.plot([r.k_lo_identifiable, r.k_hi_identifiable], [i, i], color=col, lw=7, solid_capstyle="butt",
                label=sysname if sysname not in ax.get_legend_handles_labels()[1] else None)
        ax.text(r.k_hi_identifiable * 1.08, i, f"max depth {int(r.depth_max)}", va="center", fontsize=7)
    ax.set_yticks(y, [f"{s}\n{d_}" for s, d_ in g.index]); ax.set_xscale("log")
    ax.set_xlabel("switching rate k identifiable from these trees (1/branch length)")
    ax.set_xlim(0.02, 2.0)
    ax.set_title("Neither recorder spans a quarter of one decade", fontsize=9); ax.legend(fontsize=7)
    save(fig, "m3_fig17_identifiability_windows")

# Fig 18: synthetic separation of occupancy from switching rate
if have("e10_focus.csv") or have("e10_occupancy_vs_rate.csv"):
    d = rd("e10_focus.csv") if have("e10_focus.csv") else rd("e10_occupancy_vs_rate.csv")
    fig, axs = plt.subplots(1, 3, figsize=(11.4, 3.3), sharex=True)
    order = [m for m in ["OCC", "RATE", "FROZE", "BOTH"] if m in set(d.mech)]
    x = np.arange(len(order)); w = .38
    for j, (col, tcol, name) in enumerate([("V_pi", "true_V_pi", "occupancy  V$_\\pi$"),
                                           ("V_logk", "true_V_logk", "switching  V$_{\\log k}$")]):
        if col not in d.columns:
            continue
        est = [d[d.mech == m][col].median() for m in order]
        tru = [d[d.mech == m][tcol].median() for m in order]
        axs[j].bar(x - w / 2, tru, w, color="#999", label="truth")
        axs[j].bar(x + w / 2, est, w, color=C[j], label="estimate", edgecolor="white")
        axs[j].set_xticks(x, order); axs[j].set_title(name, fontsize=9); axs[j].legend(fontsize=7)
    axs[0].set_ylabel("variance component")
    if {"lrt_pi_p", "lrt_k_p"} <= set(d.columns):
        rej = d.assign(pi=d.lrt_pi_p < .05, k=d.lrt_k_p < .05).groupby("mech")[["pi", "k"]].mean().reindex(order)
        axs[2].bar(x - w / 2, rej.pi, w, color=C[0], label="rejects shared $\\pi$")
        axs[2].bar(x + w / 2, rej.k, w, color=C[1], hatch="//", edgecolor="white", label="rejects shared k")
        truth_k = [d[d.mech == m].true_V_logk.median() > 1e-9 for m in order]  # identical rates give ~1e-33, not exactly 0
        for xi, t in zip(x, truth_k):
            axs[2].text(xi, 1.04, "k varies" if t else "k shared", ha="center", fontsize=6.5, color="#555")
        axs[2].axhline(.05, color="#888", ls=":", lw=1)
        axs[2].set_ylim(0, 1.15); axs[2].set_xticks(x, order)
        axs[2].set_ylabel("rejection rate"); axs[2].legend(fontsize=7, loc="center right")
        axs[2].set_title("likelihood-ratio tests stay calibrated", fontsize=9)
    fig.suptitle("OCC and FROZE share an occupancy distribution and differ only in timescale; "
                 "FROZE sits below the identifiability window", fontsize=8)
    save(fig, "m3_fig18_separation")

# Fig 19: real per-clone occupancy and rate
if have("m3_real_perclone.csv"):
    d = rd("m3_real_perclone.csv")
    d = d[d.get("scope").isna()] if "scope" in d.columns else d
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for j, m in enumerate((1, 2)):
        g = d[d.mouse == m].sort_values("pi")
        if not len(g):
            continue
        y = np.arange(len(g))
        lo = np.minimum(g.pi_lo, g.pi); hi = np.maximum(g.pi_hi, g.pi)   # see note in perclone.fit_clone
        axs[j].errorbar(g.pi, y, xerr=[np.clip(g.pi - lo, 0, None), np.clip(hi - g.pi, 0, None)],
                        fmt="o", color=C[0], ms=4, lw=1.2, capsize=2)
        axs[j].set_yticks(y, [f"clone {int(c)} (n={int(n)})" for c, n in zip(g.clone, g.n)], fontsize=6.5)
        axs[j].set_xlabel("estimated occupancy $\\pi_c$ (95% profile CI)")
        nb = int((~g.k_at_bound).sum())
        axs[j].set_title(f"mouse {m}: {len(g)} clones; rate off-bound in {nb}", fontsize=9)
        axs[j].set_xlim(0, 1)
    save(fig, "m3_fig19_real_perclone")

# Fig 20: held-out site prediction ladder
if have("m3_real_heldout_ladder.csv"):
    d = rd("m3_real_heldout_ladder.csv")
    cols = [c for c in d.columns if c.startswith("err_")]
    nice = {"err_p_grand": "grand mean", "err_p_pt_mean": "primary clone mean",
            "err_p_clone_mean": "clone mean (all training sites)", "err_p_tree_pi": "tree-based clone occupancy"}
    fig, ax = plt.subplots(figsize=(6.8, 3.4))
    x = np.arange(len(cols)); w = .38
    for j, m in enumerate(sorted(d.mouse.unique())):
        g = d[d.mouse == m]
        ax.bar(x + (j - .5) * w, [g[c].mean() for c in cols], w, color=C[j], hatch=["", "//"][j],
               edgecolor="white", label=f"mouse {int(m)} (all clones)")
    for j, m in enumerate(sorted(d.mouse.unique())):
        g = d[(d.mouse == m) & (d.n < d[d.mouse == m].n.max())]
        if len(g):
            ax.plot(x + (j - .5) * w, [g[c].mean() for c in cols], "k_", ms=14, mew=2,
                    label="excluding the dominant clone" if j == 0 else None)
    ax.set_xticks(x, [nice.get(c, c) for c in cols], rotation=18, ha="right")
    ax.set_ylabel("mean |error|, held-out peritoneal site"); ax.legend(fontsize=7)
    ax.set_title("Predicting an unseen metastatic site", fontsize=9)
    save(fig, "m3_fig20_heldout_ladder")

print("written:", sorted(x for x in os.listdir(F) if x.startswith("m3_") and x.endswith(".png")))

# Fig 21: equal-budget designs
if have("e11_equal_budget.csv"):
    e = rd("e11_equal_budget.csv")
    g = e.groupby(["design", "C", "cells_per_clone", "depth"]).median(numeric_only=True).reset_index()
    g = g.sort_values("C")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), sharex=True)
    x = np.arange(len(g))
    lab = [f"{int(r.C)} clones\n{int(r.cells_per_clone)} cells, depth {int(r.depth)}" for _, r in g.iterrows()]
    for ax, (est, tru, name) in zip(axs, [("V_pi", "true_V_pi", "occupancy  V$_\\pi$"),
                                          ("V_logk", "true_V_logk", "switching  V$_{\\log k}$")]):
        if est not in g.columns:
            continue
        ax.bar(x - .2, g[tru], .38, color="#999", label="truth")
        ax.bar(x + .2, g[est], .38, color=C[0], label="recovered", edgecolor="white")
        ax.set_xticks(x, lab, fontsize=6.5); ax.set_title(name, fontsize=9); ax.legend(fontsize=7)
    axs[0].set_ylabel("variance component")
    fig.suptitle("Same 1,280-cell budget spent four ways: the two components want opposite designs", fontsize=8.5)
    save(fig, "m3_fig21_equal_budget")
    print("fig21 written")
