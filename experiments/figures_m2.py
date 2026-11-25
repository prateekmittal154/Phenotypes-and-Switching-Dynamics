"""Milestone 2 figures. Usage: python figures_m2.py <tables_dir> <out_fig_dir>"""
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
RT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ST = os.path.join(RT, "results/tables")
os.makedirs(F, exist_ok=True)


def save(fig, n):
    fig.tight_layout(); fig.savefig(os.path.join(F, n + ".png"), dpi=200); fig.savefig(os.path.join(F, n + ".pdf"))
    plt.close(fig)


def have(p):
    return os.path.exists(p)


# --- Fig 8: mechanism signatures in simulation (synthetic ground truth)
p = os.path.join(ST, "e9_curves.csv")
if have(p):
    cu = pd.read_csv(p); cu = cu[cu.organ == 0]
    lab = {"A_stable": "A stable-state selection", "B_switching": "B heritable propensity",
           "C_induction": "C organ induction", "D_bottleneck": "D founder bottleneck", "E_coupled": "E coupled"}
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3), sharey=True)
    for j, (mech, l) in enumerate(lab.items()):
        g = cu[cu.mech == mech].sort_values("depth")
        if not len(g):
            continue
        axs[0].plot(g.depth, g.clone_centred, color=C[j % 5], marker=MK[j % 5], ls=["-", "--", ":", "-.", "-"][j], label=l)
        axs[1].plot(g.depth, g.site_centred, color=C[j % 5], marker=MK[j % 5], ls=["-", "--", ":", "-.", "-"][j], label=l)
    axs[0].set_title("clone-centred: plateau removed by construction", fontsize=9)
    axs[1].set_title("site-centred: plateau = heritable propensity", fontsize=9)
    for a in axs:
        a.set_xlabel("MRCA depth (deeper = closer relatives)"); a.axhline(0, color="#888", lw=1)
    axs[0].set_ylabel("state covariance"); axs[1].legend(fontsize=7)
    save(fig, "m2_fig8_mechanism_signatures")

# --- Fig 9: recovery of the plateau B against ground truth
p = os.path.join(ST, "e9_mechanism_recovery.csv")
if have(p):
    d = pd.read_csv(p); d = d[d.organ == 0]
    fig, axs = plt.subplots(1, 2, figsize=(8.6, 3.3))
    order = ["A_stable", "B_switching", "C_induction", "D_bottleneck", "E_coupled"]
    order = [o for o in order if o in set(d.mech)]
    x = np.arange(len(order))
    med = [d[d.mech == o].B.median() for o in order]
    lo = [d[d.mech == o].B.quantile(.1) for o in order]
    hi = [d[d.mech == o].B.quantile(.9) for o in order]
    tru = [d[d.mech == o].truth_B.median() for o in order]
    nai = [d[d.mech == o].clone_frac_var.median() for o in order]
    axs[0].bar(x - .2, nai, .35, color=C[1], hatch="//", edgecolor="white", label="naive: variance of clone fractions")
    axs[0].bar(x + .2, med, .35, color=C[0], yerr=[np.array(med) - lo, np.array(hi) - np.array(med)],
               capsize=3, ecolor="#444", label="decomposition estimate of B")
    axs[0].plot(x, tru, "k_", ms=22, mew=2, label="true B")
    axs[0].set_xticks(x, [o.split("_")[0] for o in order]); axs[0].set_ylabel("heritable-propensity variance")
    axs[0].legend(fontsize=7); axs[0].set_title("Primary tumour: only the decomposition is unbiased", fontsize=9)
    fp = d.assign(sig=d.B_lo > 0).groupby("mech").sig.mean().reindex(order)
    axs[1].bar(x, fp.values, .55, color=[C[3] if o in ("B_switching", "E_coupled") else C[2] for o in order],
               edgecolor="white")
    for xi, v in zip(x, fp.values):
        axs[1].text(xi, v + .02, f"{v:.2f}", ha="center", fontsize=7)
    axs[1].set_xticks(x, [o.split("_")[0] for o in order]); axs[1].set_ylim(0, 1.05)
    axs[1].set_ylabel("P(CI excludes B = 0)")
    axs[1].set_title("A, C, D: true B = 0, so any bar is a false positive\nB, E: true B > 0, so the bar is power",
                     fontsize=8)
    for xi, o in zip(x, order):
        axs[1].text(xi, 0.92, "null" if o.split("_")[0] in ("A", "C", "D") else "true B > 0",
                    ha="center", fontsize=6.5, color="#555")
    save(fig, "m2_fig9_plateau_recovery")

# --- Fig 10: real data, composition vs within-clone decomposition
p = os.path.join(T, "m2_real_c3_composition.csv")
if have(p):
    d = pd.read_csv(p)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))
    lbl = [f"M{int(r.mouse)} {r.site}" for _, r in d.iterrows()]
    x = np.arange(len(d))
    ax.bar(x - .22, d.composition, .4, color=C[0], label="clone composition (uses primary expression only)")
    ax.bar(x + .22, d.within_clone, .4, color=C[1], hatch="//", edgecolor="white", label="within-clone shift")
    ax.plot(x, d.total, "k_", ms=20, mew=2, label="total primary → metastasis difference")
    ax.axhline(0, color="#444", lw=1)
    ax.set_xticks(x, lbl, rotation=20); ax.set_ylabel("EMT score difference")
    ax.set_title("Metastatic EMT enrichment is a change in which clones are present", fontsize=9)
    ax.legend(fontsize=7)
    save(fig, "m2_fig10_composition_vs_within_clone")

# --- Fig 11: technical control
p = os.path.join(T, "m2_real_c1_technical_control.csv")
if have(p):
    d = pd.read_csv(p)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, m in zip(axs, (1, 2)):
        g = d[d.mouse == m].pivot_table(index="feature", columns="site", values="within_clone_shift_sd")
        g = g.reindex(g.mean(axis=1).sort_values().index)
        for j, s in enumerate(g.columns):
            ax.scatter(g[s], np.arange(len(g)), color=C[j % 5], marker=MK[j % 5], label=s, zorder=3)
        ax.axvline(0, color="#444", lw=1)
        hl = list(g.index).index("e")
        ax.axhspan(hl - .5, hl + .5, color=C[3], alpha=.18, zorder=0)
        ax.set_yticks(np.arange(len(g)), [("EMT score" if i == "e" else i) for i in g.index])
        ax.set_xlabel("within-clone shift (primary SD)"); ax.set_title(f"mouse {m}", fontsize=9)
    axs[1].legend(fontsize=7)
    fig.suptitle("Many unrelated programmes shift as much as EMT between compartments", fontsize=9)
    save(fig, "m2_fig11_technical_control")

# --- Fig 12: real covariance curves and the identifiability limit
p = os.path.join(T, "m2_real_curves.csv"); q = os.path.join(T, "m2_real_decomposition.csv")
if have(p) and have(q):
    cu = pd.read_csv(p); de = pd.read_csv(q)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    for ax, m in zip(axs, (1, 2)):
        g = cu[(cu.mouse == m) & (cu.site.isin(["PT", "Met", "Liver"]))]
        for j, s in enumerate(["PT", "Met", "Liver"]):
            h = g[g.site == s].sort_values("depth")
            if len(h):
                ax.plot(h.depth, h.site_centred, color=C[j], marker=MK[j], label=f"{s} (site-centred)")
                ax.plot(h.depth, h.clone_centred, color=C[j], ls=":", lw=1.4, alpha=.8)
        ax.axhline(0, color="#888", lw=1); ax.set_xlabel("MRCA depth (shared indels)")
        ax.set_title(f"mouse {m}", fontsize=9)
        r = de[(de.mouse == m) & (de.site == "PT")]
        if len(r):
            rr = r.iloc[0]
            ax.annotate(f"PT plateau B = {rr.B:.3f}\n95% CI [{rr.B_lo:.3f}, {rr.B_hi:.3f}]\nstep-1 k = {rr.k:.2f}",
                        xy=(.03, .72), xycoords="axes fraction", fontsize=7,
                        bbox=dict(fc="white", ec="#ccc", lw=.6))
    axs[0].set_ylabel("state covariance")
    axs[1].legend(fontsize=7)
    fig.suptitle("Solid = site-centred (contains the plateau); dotted = clone-centred (plateau removed)", fontsize=8)
    save(fig, "m2_fig12_real_curves")

# --- Fig 13: model comparison and held-out prediction
p = os.path.join(T, "m2_real_model_comparison.csv"); q = os.path.join(T, "m2_real_heldout_site.csv")
if have(p):
    d = pd.read_csv(p)
    fig, axs = plt.subplots(1, 2, figsize=(8.8, 3.4))
    order = ["M0", "M1", "M2", "M4"]
    nm = {"M0": "M0 single state", "M1": "M1 clone propensity", "M2": "M2 organ induction",
          "M4": "M4 clone + organ"}
    w = .38
    for j, m in enumerate((1, 2)):
        g = d[d.mouse == m].set_index("model").reindex(order)
        axs[0].bar(np.arange(4) + (j - .5) * w, g.delta_aic, w, color=C[j], label=f"mouse {m}",
                   hatch=["", "//"][j], edgecolor="white")
    axs[0].set_xticks(np.arange(4), [nm[o] for o in order], rotation=18, ha="right")
    axs[0].set_ylabel("ΔAIC from best (lower = better)"); axs[0].legend(fontsize=7)
    axs[0].set_title("Clone-specific propensity dominates organ induction", fontsize=9)
    if have(q):
        h = pd.read_csv(q)
        axs[1].scatter(h.observed, h.pred_M0, color=C[1], marker="s", label="M0 single state", zorder=3)
        axs[1].scatter(h.observed, h.pred_M1, color=C[0], marker="o", label="M1 clone propensity", zorder=3)
        axs[1].plot([0, 1], [0, 1], color="#888", lw=1, ls="--")
        axs[1].set_xlabel("observed high-EMT fraction, held-out peritoneal site")
        axs[1].set_ylabel("predicted"); axs[1].legend(fontsize=7)
        e0, e1 = h.err_M0.mean(), h.err_M1.mean()
        axs[1].set_title(f"held-out error {e0:.2f} → {e1:.2f}", fontsize=9)
    save(fig, "m2_fig13_model_comparison")

# --- Fig 14: stress tests
p = os.path.join(ST, "e9_stress.csv")
if have(p):
    d = pd.read_csv(p); d = d[d.organ == 0]
    if len(d):
        fig, ax = plt.subplots(figsize=(6.2, 3.4))
        g = d.groupby(["eps", "clone_loss", "n_sample"]).B.median().reset_index()
        for j, (cl, ns) in enumerate([(0.0, 4000), (0.3, 4000), (0.0, 800), (0.3, 800)]):
            s = g[(g.clone_loss == cl) & (g.n_sample == ns)].sort_values("eps")
            if len(s):
                ax.plot(s.eps, s.B, color=C[j], marker=MK[j],
                        label=f"{ns} cells, {int(cl * 100)}% clone labels lost")
        t = d.truth_B.median()
        ax.axhline(t, color="#333", ls="--", lw=1, label=f"true B = {t:.3f}")
        ax.set_xlabel("state misclassification rate"); ax.set_ylabel("estimated plateau B")
        ax.legend(fontsize=7); ax.set_title("Misclassification attenuates B; it does not create it", fontsize=9)
        save(fig, "m2_fig14_stress")

print("figures written:", sorted(os.listdir(F)))
