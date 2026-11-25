"""Benchmark figures B1-B5, built only from executed experiment outputs.

Usage: python -I make_bench_figures.py <benchdir>
"""
import os, sys, json, glob
# If matplotlib's dependencies live in a separate site directory, point CLDYN_SITE
# at it; python -I drops PYTHONPATH, so it has to be added here.
_site = os.environ.get("CLDYN_SITE")
if _site and os.path.isdir(_site):
    sys.path.insert(0, _site)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

B = sys.argv[1]
MET = os.path.join(B, "results", "metrics")
TAB = os.path.join(B, "results", "tables")
FIG = os.path.join(B, "results", "figures")
os.makedirs(FIG, exist_ok=True)

C = dict(studyd="#0072B2", path="#E69F00", pastri="#009E73", scout="#CC79A7",
         grey="#8C8C8C", truth="#4D4D4D", red="#D55E00", sky="#56B4E9")
MCOL = {"Study D": C["studyd"], "PATH": C["path"], "PASTRI": C["pastri"],
        "SCOUT": C["scout"], "SCOUT (theta)": C["scout"], "SCOUT (alpha)": C["scout"]}
MM = 1 / 25.4
W2 = 180 * MM

plt.rcParams.update({
    "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.color": "#E8E8E8",
    "grid.linewidth": 0.5, "lines.linewidth": 1.3, "lines.markersize": 3.5,
    "axes.linewidth": 0.7, "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
})

SCEN_SHORT = {"S1": "stable\nstates", "S2": "rate\nvaries", "S3": "occupancy\nvaries",
              "S4": "both\nvary", "S5": "induction", "S6": "bottleneck",
              "S7": "unequal\nexpansion", "S8": "measurement\nerror"}


def panel(ax, letter, dx=-0.18, dy=1.07):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left")


def noaxis(ax):
    ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("wrote", name)


# ===================================================================== B1
def figB1():
    """Method and estimand comparison. A matrix, deliberately not a ranking."""
    fig = plt.figure(figsize=(W2, 108 * MM))
    ax = fig.add_axes([0, 0, 1, 1]); noaxis(ax)
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)

    meths = ["Study D", "PATH", "PASTRI", "SCOUT"]
    cols = [C["studyd"], C["path"], C["pastri"], C["scout"]]
    rows = [
        ("Native estimand", ["per-clone $\\pi_c$, $k_c$;\nvariance components",
                             "phylogenetic\nautocorrelation; $Q$",
                             "per-cycle transition\nmatrix $T$",
                             "OU $\\alpha$, $\\theta$, $\\sigma^2$"]),
        ("Phenotype", ["binary", "categorical /\ncontinuous", "categorical", "continuous"]),
        ("Uses branch\nlengths", ["yes", "via mean\npatristic distance", "no\n(topological depth)",
                                  "yes, or unit\nedges if absent"]),
        ("Rates vary\nwith depth", ["no", "no", "yes\n(by design)", "no"]),
        ("Uncertainty\nreported", ["profile intervals;\ncalibrated tests",
                                   "z-test for\ncorrelation only", "permutation test\n(no rate CI)",
                                   "likelihood / AICc"]),
        ("States when\nthe data cannot\nidentify a rate", ["reports the\nidentifiable window",
                                                           "returns a point\nestimate",
                                                           "returns a point\nestimate",
                                                           "returns a point\nestimate"]),
    ]
    x0, w = 22.0, 19.0
    ax.text(0, 96, "What each method estimates, and what it assumes",
            fontsize=8.5, fontweight="bold", va="top")
    ax.text(0, 91.5, "A compatibility matrix, not a ranking: the four methods target different quantities.",
            fontsize=6.3, color="#555", va="top")
    for j, (m, c) in enumerate(zip(meths, cols)):
        ax.add_patch(FancyBboxPatch((x0 + j * w, 80.5), w - 1.2, 6.0,
                                    boxstyle="round,pad=0.3", fc=c, alpha=0.18, ec=c, lw=0.8))
        ax.text(x0 + j * w + (w - 1.2) / 2, 83.5, m, fontsize=7.4, fontweight="bold",
                color=c, ha="center", va="center")
    y = 76.0
    for i, (lab, vals) in enumerate(rows):
        h = 11.5
        if i % 2 == 0:
            ax.add_patch(Rectangle((0, y - h + 1), 100, h - 0.6, fc="#F7F7F7", ec="none"))
        ax.text(0.5, y - h / 2 + 1, lab, fontsize=6.6, fontweight="bold", va="center",
                color="#333")
        for j, v in enumerate(vals):
            ax.text(x0 + j * w + (w - 1.2) / 2, y - h / 2 + 1, v, fontsize=6.0,
                    ha="center", va="center", color="#222")
        y -= h
    ax.plot([0, 100], [79.5, 79.5], color="#CCC", lw=0.8)
    save(fig, "Figure_B1_estimand_matrix")


# ===================================================================== B2
def figB2():
    """Ground-truth recovery across contrasting evolutionary mechanisms."""
    b3 = pd.read_csv(os.path.join(TAB, "TableB3_clone_level_true_tree.csv"))
    vs = pd.read_csv(os.path.join(TAB, "TableB3b_variance_components.csv"))
    vs = vs[vs.tree_variant == "true"]

    fig = plt.figure(figsize=(W2, 128 * MM))
    gs = fig.add_gridspec(3, 1, hspace=0.75, height_ratios=[1, 1, 1])

    # --- A: occupancy correlation with truth
    ax = fig.add_subplot(gs[0]); panel(ax, "A", dx=-0.07)
    occ = b3[b3.quantity == "occupancy"]
    scen = sorted(occ.scenario.unique())
    meths = ["Study D", "PATH", "PASTRI", "SCOUT (theta)"]
    x = np.arange(len(scen)); w = 0.2
    for j, m in enumerate(meths):
        v = [occ[(occ.scenario == s) & (occ.method == m)].corr_with_truth.mean() for s in scen]
        ax.bar(x + (j - 1.5) * w, v, w, color=MCOL[m], edgecolor="white", lw=0.4,
               label=m.replace(" (theta)", " ($\\theta$)"))
    ax.axhline(0, color="#444", lw=0.7)
    if "S2" in scen:
        ax.text(scen.index("S2"), 0.42, "no occupancy\nvariation to detect\n(true $V_\\pi$ = 0)",
                fontsize=5.6, ha="center", va="center", color="#777", style="italic")
    ax.set_xticks(x, [SCEN_SHORT[s] for s in scen], fontsize=5.8)
    ax.set_ylabel("correlation with\ntrue occupancy"); ax.set_ylim(-0.4, 1.15)
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Every method recovers clone occupancy", fontsize=7.4, pad=18)

    # --- B: switching-rate bias
    ax = fig.add_subplot(gs[1]); panel(ax, "B", dx=-0.07)
    rate = b3[b3.quantity == "switching rate"]
    meths = ["Study D", "PATH", "PASTRI"]
    for j, m in enumerate(meths):
        v, lo, hi = [], [], []
        for s in scen:
            r = rate[(rate.scenario == s) & (rate.method == m)]
            v.append(r.median_log_bias.mean() if len(r) else np.nan)
            lo.append(r.bias_lo.mean() if len(r) else np.nan)
            hi.append(r.bias_hi.mean() if len(r) else np.nan)
        v = np.array(v, float); lo = np.array(lo, float); hi = np.array(hi, float)
        err = np.vstack([np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)])
        ax.bar(x + (j - 1) * 0.26, v, 0.26, color=MCOL[m], edgecolor="white", lw=0.4,
               label=m, yerr=err, capsize=1.6, error_kw=dict(lw=0.6, ecolor="#444"))
    ax.axhline(0, color="#333", lw=0.9)
    ax.set_xticks(x, [SCEN_SHORT[s] for s in scen], fontsize=5.8)
    ax.set_ylabel("median log bias\nin switching rate")
    ax.legend(frameon=False, ncol=3, loc="upper left")
    ax.set_title("No method recovers the switching rate without bias; "
                 "zero is unbiased", fontsize=7.4)

    # --- C: variance components, the quantity actually reported in papers
    ax = fig.add_subplot(gs[2]); panel(ax, "C", dx=-0.07)
    key = vs[vs.scenario.isin(["S2", "S3", "S4"])]
    labs, xs, i = [], [], 0
    for s in ["S3", "S2", "S4"]:
        g = key[key.scenario == s]
        tv = g.true_V_pi.iloc[0] if len(g) else np.nan
        ax.bar(i, tv, 0.6, color=C["truth"], edgecolor="white", lw=0.4,
               label="truth" if i == 0 else None)
        labs.append("truth"); xs.append(i); i += 1
        for m in ["Study D", "PATH", "PASTRI"]:
            r = g[g.method == m]
            ax.bar(i, r.est_V_pi.iloc[0] if len(r) else np.nan, 0.6, color=MCOL[m],
                   edgecolor="white", lw=0.4, label=m if s == "S3" else None)
            labs.append(m.replace("Study D", "StudyD")); xs.append(i); i += 1
        i += 0.8
    ax.set_xticks(xs, labs, rotation=45, ha="right", fontsize=5.5)
    ax.set_ylabel("$V_\\pi$  (occupancy\nvariance across clones)")
    ymax = 0.075
    ax.set_ylim(0, ymax)
    for pos, lab2 in zip([1.5, 6.3, 11.1], ["occupancy varies (S3)", "rate varies only (S2)",
                                            "both vary (S4)"]):
        ax.text(pos, ymax * 0.93, lab2, fontsize=6.0, ha="center", va="top",
                color="#555", fontweight="bold")
    # the S2 truth value is exactly zero; mark it so the gap is not read as missing data
    ax.annotate("true value 0", xy=(5, 0), xytext=(5, ymax * 0.30), fontsize=5.4,
                ha="center", color="#555",
                arrowprops=dict(arrowstyle="-|>", lw=0.6, color="#888"))
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.06))
    ax.set_title("Occupancy variance is recovered by all three; the middle group shows a "
                 "rate-only mechanism", fontsize=7.4, pad=20)
    save(fig, "Figure_B2_ground_truth_recovery")


# ===================================================================== B3
def figB3():
    """Information and identifiability: PASTRI depth sweep + clone-size behaviour."""
    sw_files = glob.glob(os.path.join(MET, "gt_*_pastri_depth_sweep.csv"))
    fig = plt.figure(figsize=(W2, 62 * MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42)

    # --- A: PASTRI rate vs depth window
    ax = fig.add_subplot(gs[0]); panel(ax, "A")
    scored_f = os.path.join(MET, "pastri_depth_sweep_scored.csv")
    sw = None
    if sw_files:
        sw = pd.concat([pd.read_csv(f) for f in sw_files], ignore_index=True)
        sw["tree_variant"] = sw["tree_variant"].astype(str).str.strip().str.lower()
        sw["clone"] = sw["clone"].astype(str)
        sw["rep_dir"] = sw["rep_dir"].astype(str).str.rstrip("/")
        sw = sw[sw.tree_variant == "true"]
        tr = []
        for d in sw.rep_dir.unique():
            tj = os.path.join(d, "truth.json")
            if os.path.exists(tj):
                t = json.load(open(tj))
                for c in t["clones"]:
                    tr.append(dict(rep_dir=d, clone=f"C{c['clone']}", true_k=c["k"]))
        trdf = pd.DataFrame(tr)
        if len(trdf):
            trdf["rep_dir"] = trdf["rep_dir"].astype(str).str.rstrip("/")
            trdf["clone"] = trdf["clone"].astype(str)
            sw = sw.merge(trdf, on=["rep_dir", "clone"], how="left")
            sw["lr"] = np.log(sw.k.replace([np.inf, -np.inf], np.nan)) - np.log(sw.true_k)
            sw.to_csv(scored_f, index=False)
        else:
            # the replicate directories are not present (they are regenerable and are not
            # shipped); fall back to the scored table committed with the repository
            sw = None
    if sw is None and os.path.exists(scored_f):
        sw = pd.read_csv(scored_f)
    if sw is not None and "lr" in sw.columns:
        g = sw.groupby("fi_depth").lr.agg(["median", "count"]).reset_index()
        q = sw.groupby("fi_depth").lr.quantile([0.25, 0.75]).unstack()
        ax.plot(g.fi_depth, g["median"], "o-", color=C["pastri"])
        ax.fill_between(g.fi_depth, q[0.25], q[0.75], color=C["pastri"], alpha=0.18, lw=0)
        ax.axhline(0, color="#333", lw=0.9, ls="--")
        ax.set_xlabel("PASTRI depth window  $d \\leq$")
        ax.set_ylabel("log bias in rate")
        ax.set_title("PASTRI: bias grows as more\ndistal pairs enter", fontsize=7)

    # --- B: boundary / no-estimate rates
    ax = fig.add_subplot(gs[1]); panel(ax, "B")
    fr = pd.read_csv(os.path.join(TAB, "TableB4b_no_estimate_rates.csv"))
    fr = fr[(fr.tree_variant == "true") & (fr.method != "Study D (k at bound)")]
    agg = fr.groupby("method").agg(n=("n", "sum"), k=("n_no_estimate", "sum")).reset_index()
    agg["frac"] = agg.k / agg.n
    order = ["Study D", "PATH", "PASTRI", "SCOUT"]
    agg = agg.set_index("method").reindex(order).reset_index()
    ax.bar(np.arange(len(agg)), agg.frac, 0.6,
           color=[MCOL[m] for m in agg.method], edgecolor="white", lw=0.4)
    top = float(np.nanmax(agg.frac)) if np.isfinite(np.nanmax(agg.frac)) else 0.1
    ax.set_ylim(0, top * 1.38)
    for i, (f_, k_, n_) in enumerate(zip(agg.frac, agg.k, agg.n)):
        ax.text(i, f_ + top * 0.04, f"{int(k_)}/{int(n_)}", ha="center", fontsize=5.4)
    ax.set_xticks(np.arange(len(agg)), agg.method, fontsize=6, rotation=20, ha="right")
    ax.set_ylabel("clones with no\nfinite estimate")
    ax.set_title("Where a method declines\nto return a value", fontsize=7)

    # --- C: rate estimate vs true rate, pooled
    ax = fig.add_subplot(gs[2]); panel(ax, "C")
    cl = pd.read_csv(os.path.join(MET, "clone_level_all.csv"))
    cl = cl[cl.tree_variant == "true"]
    for m, col in [("Study D", "StudyD_k"), ("PATH", "PATH_k"), ("PASTRI", "PASTRI_k")]:
        v = pd.to_numeric(cl[col], errors="coerce")
        ok = np.isfinite(v) & (v > 0) & np.isfinite(cl.true_k)
        if ok.sum() < 5:
            continue
        bins = pd.qcut(cl.true_k[ok], 5, duplicates="drop")
        gg = pd.DataFrame(dict(t=cl.true_k[ok], e=v[ok], b=bins)).groupby("b", observed=True)
        ax.plot(gg.t.median(), gg.e.median(), "o-", color=MCOL[m], label=m)
    lims = [0.015, 6]
    ax.plot(lims, lims, color="#333", ls="--", lw=0.9, label="unbiased")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(*lims); ax.set_ylim(0.01, 60)
    ax.set_xlabel("true switching rate $k$")
    ax.set_ylabel("estimated rate")
    ax.legend(frameon=False, fontsize=5.4, loc="upper left")
    ax.set_title("Estimated against true rate\n(median per quintile)", fontsize=7)
    save(fig, "Figure_B3_information_identifiability")


# ===================================================================== B4
def figB4():
    """Clock and reconstruction sensitivity."""
    b4 = pd.read_csv(os.path.join(TAB, "TableB4_clock_sensitivity.csv"))
    cl = pd.read_csv(os.path.join(MET, "clone_level_all.csv"))
    fig = plt.figure(figsize=(W2, 62 * MM))
    gs = fig.add_gridspec(1, 3, wspace=0.42)

    # --- A: log ratio unit vs true edges, per method
    ax = fig.add_subplot(gs[0]); panel(ax, "A")
    r = b4[b4.quantity == "rate / decay"].set_index("method")
    order = [m for m in ["Study D", "PATH", "PASTRI", "SCOUT"] if m in r.index]
    v = r.loc[order, "median_log_ratio_unit_vs_true"].values
    lo = r.loc[order, "ci_lo"].values; hi = r.loc[order, "ci_hi"].values
    err = np.vstack([np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)])
    ax.bar(np.arange(len(order)), v, 0.6, color=[MCOL[m] for m in order],
           edgecolor="white", lw=0.4, yerr=err, capsize=2, error_kw=dict(lw=0.7, ecolor="#444"))
    ax.axhline(0, color="#333", lw=0.9)
    for i, vv in enumerate(v):
        if np.isfinite(vv):
            lab = "0 (unchanged)" if abs(vv) < 1e-9 else f"{np.exp(vv):.1f}x"
            ax.text(i, max(vv, 0) + 0.06, lab, ha="center", fontsize=5.4, color="#333")
    ax.set_xticks(np.arange(len(order)), order, rotation=20, ha="right", fontsize=6)
    ax.set_ylabel("log( unit edges / true times )")
    ax.set_ylim(0, max(1.1 * np.nanmax(hi), 2.0))
    ax.set_title("Effect of discarding\nbranch times", fontsize=7)

    # --- B: fraction of estimates literally unchanged
    ax = fig.add_subplot(gs[1]); panel(ax, "B")
    f = r.loc[order, "frac_identical"].values
    ax.bar(np.arange(len(order)), f, 0.6, color=[MCOL[m] for m in order],
           edgecolor="white", lw=0.4)
    for i, vv in enumerate(f):
        ax.text(i, vv + 0.02, f"{vv:.2f}", ha="center", fontsize=5.6)
    ax.set_xticks(np.arange(len(order)), order, rotation=20, ha="right", fontsize=6)
    ax.set_ylabel("fraction of estimates\nbit-identical")
    ax.set_ylim(0, 1.12)
    ax.set_title("PASTRI cannot see branch times\nby construction", fontsize=7)

    # --- C: occupancy is insensitive to the same change
    ax = fig.add_subplot(gs[2]); panel(ax, "C")
    ro = b4[b4.quantity == "occupancy"].set_index("method")
    order2 = [m for m in ["Study D", "PATH", "PASTRI"] if m in ro.index]
    v2 = ro.loc[order2, "median_abs_diff_unit_vs_true"].values
    ax.bar(np.arange(len(order2)), v2, 0.6, color=[MCOL[m] for m in order2],
           edgecolor="white", lw=0.4)
    for i, vv in enumerate(v2):
        if np.isfinite(vv):
            ax.text(i, vv + 0.00008, "0.000" if vv < 1e-9 else f"{vv:.4f}",
                    ha="center", fontsize=5.4, color="#333")
    ax.set_xticks(np.arange(len(order2)), order2, rotation=20, ha="right", fontsize=6)
    ax.set_ylabel("median |change| in\noccupancy estimate")
    ax.set_ylim(0, max(1.25 * np.nanmax(v2), 1e-4))
    ax.set_title("Occupancy barely moves", fontsize=7)
    save(fig, "Figure_B4_clock_sensitivity")


# ===================================================================== B5
def figB5():
    """Real-data comparison at animal level, plus cross-method agreement."""
    f = os.path.join(TAB, "TableB5_real_data.csv")
    if not os.path.exists(f):
        print("skipping B5: no real-data table yet")
        return
    d = pd.read_csv(f)
    ag_f = os.path.join(TAB, "TableB5b_cross_method_agreement.csv")
    ag = pd.read_csv(ag_f) if os.path.exists(ag_f) else None

    fig = plt.figure(figsize=(W2, 150 * MM))
    gs = fig.add_gridspec(3, 1, hspace=0.95, height_ratios=[1, 1, 1.1])

    order = ["Simeonov M1", "Simeonov M2", "Quinn 5k", "Quinn 10k", "Quinn 30k", "Quinn 100k"]
    animals = [a for a in order if a in set(d.animal)]
    x = np.arange(len(animals))

    # --- A: occupancy variance per animal
    ax = fig.add_subplot(gs[0]); panel(ax, "A", dx=-0.105, dy=1.12)
    meths = ["Study D", "PATH", "PASTRI", "SCOUT"]
    w = 0.2
    for j, m in enumerate(meths):
        v = [d[(d.animal == a) & (d.method == m)].V_pi_scaled.mean() for a in animals]
        ax.bar(x + (j - 1.5) * w, v, w, color=MCOL[m], edgecolor="white", lw=0.4, label=m)
    ax.set_xticks(x, animals, rotation=18, ha="right", fontsize=6)
    ax.set_ylabel("across-clone occupancy\nvariance (scaled per method)")
    ax.legend(frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.16))
    if "Quinn 30k" in animals:
        i30 = animals.index("Quinn 30k")
        ax.text(i30, 0.004, "1 clone:\nno across-clone\nvariance", fontsize=5.4,
                ha="center", va="bottom", color="#777", style="italic")
    ax.set_title("All four methods find clone-associated occupancy structure "
                 "in every eligible animal", fontsize=7.4, pad=22)

    # --- B: spread of rate estimates per animal
    ax = fig.add_subplot(gs[1]); panel(ax, "B", dx=-0.105, dy=1.12)
    for j, m in enumerate(["Study D", "PATH", "PASTRI"]):
        v = [d[(d.animal == a) & (d.method == m)].rate_iqr_log.mean() for a in animals]
        ax.bar(x + (j - 1) * 0.26, v, 0.26, color=MCOL[m], edgecolor="white", lw=0.4, label=m)
    ax.set_xticks(x, animals, rotation=18, ha="right", fontsize=6)
    ax.set_ylabel("IQR of log switching rate\nacross clones")
    ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.02))
    ax.set_title("The same cells give a three- to six-fold wider spread of rates under "
                 "Study D than under PATH", fontsize=7.4, pad=18)

    # --- C: cross-method rank agreement, the headline
    ax = fig.add_subplot(gs[2]); panel(ax, "C", dx=-0.105, dy=1.12)
    if ag is not None and len(ag):
        piv = ag.pivot_table(index="pair", columns="quantity", values="spearman",
                             aggfunc="mean")
        piv = piv.sort_values("occupancy", ascending=False)
        xp = np.arange(len(piv))
        ax.bar(xp - 0.2, piv.get("occupancy"), 0.38, color=C["studyd"],
               edgecolor="white", lw=0.4, label="occupancy")
        ax.bar(xp + 0.2, piv.get("rate"), 0.38, color=C["red"],
               edgecolor="white", lw=0.4, label="switching rate")
        mo = float(np.nanmedian(piv["occupancy"])); mr = float(np.nanmedian(piv["rate"]))
        ax.axhline(mo, color=C["studyd"], ls=":", lw=0.9)
        ax.axhline(mr, color=C["red"], ls=":", lw=0.9)
        ax.text(len(piv) - 0.4, mo + 0.02, f"median {mo:.2f}", fontsize=5.6,
                color=C["studyd"], ha="right")
        ax.text(len(piv) - 0.4, mr + 0.02, f"median {mr:.2f}", fontsize=5.6,
                color=C["red"], ha="right")
        ax.set_xticks(xp, piv.index, rotation=22, ha="right", fontsize=5.8)
        ax.set_ylabel("Spearman rank agreement\nbetween the two methods")
        ax.set_ylim(0, 1.08)
        ax.legend(frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, 1.02))
        ax.set_title("On the same cells, methods agree about occupancy and disagree "
                     "about kinetics", fontsize=7.4, pad=18)
    save(fig, "Figure_B5_real_data")


if __name__ == "__main__":
    figB1()
    try:
        figB2()
    except Exception as e:
        print("B2 failed:", e)
    try:
        figB3()
    except Exception as e:
        print("B3 failed:", e)
    try:
        figB4()
    except Exception as e:
        print("B4 failed:", e)
    try:
        figB5()
    except Exception as e:
        print("B5 failed:", e)
    print("figures done")
