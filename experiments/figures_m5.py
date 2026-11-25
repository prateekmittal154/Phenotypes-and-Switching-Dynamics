"""Milestone 5 figures. Usage: python figures_m5.py <tables_dir> <out_fig_dir>"""
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


LEVELS = [("k_L1", "T1 true tree\n+ true lengths"), ("k_L2", "T1 topology\nunit edges"),
          ("k_L3_greedy", "T3 greedy"), ("k_L3_nj", "T3 NJ"), ("k_L3_upgma", "T3 UPGMA"),
          ("k_L3_ens", "T3 ensemble")]

# Fig 27: where switching information is lost
if have("e17_hierarchy.csv"):
    d = rd("e17_hierarchy.csv")
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
    cols = [c for c, _ in LEVELS if c in d.columns]
    labs = [l for c, l in LEVELS if c in d.columns]
    data = [np.log(d[c].values / d.k_true.values) for c in cols]
    bp = axs[0].boxplot(data, tick_labels=labs, showfliers=False, patch_artist=True)
    for i, b in enumerate(bp["boxes"]):
        b.set_facecolor(C[0] if i < 2 else C[1]); b.set_alpha(.75)
    axs[0].axhline(0, color="#333", lw=1.2)
    axs[0].set_ylabel("log( $\\hat k$ / true k )"); axs[0].tick_params(axis="x", labelrotation=30, labelsize=7)
    axs[0].set_title("Bias in the switching rate at each observation level", fontsize=9)
    if have("e17_paired_ratios.csv"):
        p = rd("e17_paired_ratios.csv")
        x = np.arange(len(p))
        axs[1].bar(x, p.median_log_ratio, .55, yerr=p.iqr / 2, capsize=3, ecolor="#444",
                   color=[C[0] if v == "L2" else C[1] for v in p.level], edgecolor="white")
        axs[1].axhline(0, color="#333", lw=1.2)
        axs[1].set_xticks(x, p.level, rotation=20, fontsize=7.5)
        axs[1].set_ylabel("paired log( $\\hat k$ / $\\hat k_{T1}$ )")
        axs[1].set_title("Paired change against the true-tree fit (IQR)", fontsize=9)
    save(fig, "m5_fig27_information_loss")

# Fig 28: tree geometry, truth versus reconstruction
if have("e17_hierarchy.csv"):
    d = rd("e17_hierarchy.csv")
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))
    groups = ["T1", "greedy", "nj", "upgma"]
    for ax, stat, lab in zip(axs, ["mrca_norm", "colless", "close_frac"],
                             ["mean MRCA depth / tree height", "Colless imbalance",
                              "fraction of close-relative pairs"]):
        cols = [f"{stat}_{g}" for g in groups if f"{stat}_{g}" in d.columns]
        if not cols:
            continue
        vals = [d[c].dropna().values for c in cols]
        bp = ax.boxplot(vals, tick_labels=[c.split("_")[-1] for c in cols], showfliers=False,
                        patch_artist=True)
        for i, b in enumerate(bp["boxes"]):
            b.set_facecolor("#999" if i == 0 else C[1]); b.set_alpha(.8)
        ax.set_ylabel(lab, fontsize=8); ax.tick_params(axis="x", labelsize=7.5)
    axs[1].set_title("Grey = true sampled tree; orange = reconstructions", fontsize=9)
    save(fig, "m5_fig28_tree_geometry")

# Fig 29: topology recovery and recorder state
if have("e17_hierarchy.csv"):
    d = rd("e17_hierarchy.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.4))
    rc = [c for c in d.columns if c.startswith("rf_recall")]
    if rc:
        g = d.groupby("recorder")[rc].median()
        x = np.arange(len(rc)); w = .38
        for j, r in enumerate(g.index):
            axs[0].bar(x + (j - .5) * w, g.loc[r].values, w, color=C[j], hatch=["", "//"][j],
                       edgecolor="white", label=r)
        axs[0].set_xticks(x, [c.replace("rf_recall_", "") for c in rc])
        axs[0].set_ylabel("fraction of true clades recovered"); axs[0].set_ylim(0, 1)
        axs[0].legend(fontsize=7); axs[0].set_title("Topology recovery", fontsize=9)
    s = d.groupby("recorder")[["saturation", "missing"]].median()
    x2 = np.arange(len(s))
    axs[1].bar(x2 - .2, s.saturation, .38, color=C[2], label="saturation", edgecolor="white")
    axs[1].bar(x2 + .2, s["missing"], .38, color=C[3], hatch="//", label="missing", edgecolor="white")
    axs[1].set_xticks(x2, s.index, fontsize=8); axs[1].set_ylim(0, 1.05); axs[1].legend(fontsize=7)
    axs[1].set_title("Recorder state at readout", fontsize=9)
    save(fig, "m5_fig29_recorder_state")

# Fig 30: real data, estimates across published reconstructions
if have("m5_quinn_reconstruction_fits.csv"):
    d = rd("m5_quinn_reconstruction_fits.csv")
    fig, axs = plt.subplots(1, 2, figsize=(9.8, 3.6))
    order = [v for v in ["hybrid", "nj", "nj_priors", "ensemble"] if v in set(d.variant)]
    piv = d.pivot_table(index=["mouse", "clone"], columns="variant", values="pi")
    if "hybrid" in piv.columns:
        for j, v in enumerate([o for o in order if o != "hybrid"]):
            if v in piv.columns:
                axs[0].scatter(piv["hybrid"], piv[v], s=22, color=C[j], marker=MK[j], label=v, alpha=.8)
        axs[0].plot([0, 1], [0, 1], color="#888", ls="--", lw=1)
        axs[0].set_xlabel("occupancy $\\hat\\pi_c$, authors' primary tree")
        axs[0].set_ylabel("$\\hat\\pi_c$, alternative tree"); axs[0].legend(fontsize=7)
        axs[0].set_title("Occupancy is stable across reconstructions", fontsize=9)
    if have("m5_quinn_variance_by_reconstruction.csv"):
        v = rd("m5_quinn_variance_by_reconstruction.csv")
        mice = sorted(v.mouse.unique()); x = np.arange(len(mice)); w = .8 / max(len(order), 1)
        for j, var in enumerate(order):
            vals = [v[(v.mouse == m) & (v.variant == var)].V_pi_raw.mean() for m in mice]
            axs[1].bar(x + (j - (len(order) - 1) / 2) * w, vals, w, color=C[j % 5], edgecolor="white",
                       label=var)
        axs[1].set_xticks(x, mice, fontsize=8); axs[1].set_ylabel("$V_\\pi$ (raw)")
        axs[1].legend(fontsize=6.5); axs[1].set_title("Occupancy variance by reconstruction", fontsize=9)
    save(fig, "m5_fig30_real_reconstruction")

print("written:", sorted(x for x in os.listdir(F) if x.startswith("m5_") and x.endswith(".png")))
