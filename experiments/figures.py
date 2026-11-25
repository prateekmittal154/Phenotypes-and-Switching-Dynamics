"""Milestone 1 figures (only for executed analyses). Okabe-Ito colours plus a marker/linestyle per series."""
import json, os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.bounds.core import tau_min, tv, support_mismatch

T, F = os.path.join(ROOT, "results/tables"), os.path.join(ROOT, "results/figures")
os.makedirs(F, exist_ok=True)
C = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9"]
MK = ["o", "s", "^", "D", "v"]
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e6e6e6", "grid.linewidth": 0.6, "lines.linewidth": 2, "lines.markersize": 6})


def save(fig, name):
    fig.tight_layout(); fig.savefig(os.path.join(F, name + ".png"), dpi=200); fig.savefig(os.path.join(F, name + ".pdf"))
    plt.close(fig)


# Fig 1: identified set for tau as a function of Lambda (counterexample C1)
P = np.array([0.90, 0.09, 0.01]); Q = np.array([0.40, 0.30, 0.30])
L = np.logspace(0, 2.3, 300); t = [tau_min(P, Q, l) for l in L]
fig, ax = plt.subplots(figsize=(5.2, 3.4))
ax.fill_between(L, t, 1, color=C[0], alpha=.12, label="identified set [τ_min(Λ), 1]")
ax.plot(L, t, color=C[0], label="sharp lower bound τ_min(Λ)")
ax.axhline(tv(P, Q), color=C[1], ls="--", label=f"TV(P,Q) = {tv(P, Q):.2f} (Λ = 1)")
ax.axhline(support_mismatch(P, Q), color=C[2], ls=":", label="support mismatch = 0 (Λ = ∞)")
ax.scatter([67.5], [0], color=C[3], marker="D", zorder=5, label="pure selection, τ = 0 (needs Λ ≥ 67.5)")
ax.scatter([1], [0.61], color=C[1], marker="s", zorder=5, label="pure induction, τ = 0.61 (Λ = 1)")
ax.set_xscale("log"); ax.set_xlabel("assumed bound on relative selection Λ"); ax.set_ylabel("fraction of metastatic cells that changed state, τ")
ax.set_title("Same (P, Q), many mechanisms: τ is only bounded given Λ", fontsize=9); ax.legend(fontsize=7, loc="upper right")
save(fig, "m1_fig1_identified_set")

# Fig 2: false induction attribution under pure selection
e4 = pd.read_csv(os.path.join(T, "e4_scenarios.csv"))
fig, axs = plt.subplots(1, 2, figsize=(8, 3.2), sharey=True)
for ax, sc in zip(axs, ["S1", "S4_p0.0005"]):
    d = e4[e4.scenario == sc]
    for j, (e, lab) in enumerate([("tv", "TV(P̂,Q̂) > 0"), ("sm", "novel-state mass > 0"), ("plug", "plug-in τ_min(Λ) > 0"),
                                  ("robust", "robust lower bound τ_L > 0")]):
        ax.plot(d.n, d[f"{e}_pos"], color=C[j], marker=MK[j], label=lab)
    ax.axhline(0.05, color="#888", lw=1, ls=":")
    ax.set_xscale("log"); ax.set_xlabel("cells sampled per site (n_P = n_Q)")
    ax.set_title({"S1": "S1 pure selection (p = 0.005, Λ = 50)", "S4_p0.0005": "S4 rare state p = 5e-4 (Λ = 857)"}[sc], fontsize=9)
axs[0].set_ylabel("P(claim induction) when τ = 0"); axs[1].legend(fontsize=7)
save(fig, "m1_fig2_false_attribution")

# Fig 3: robust bound vs truth
fig, ax = plt.subplots(figsize=(5.2, 3.4))
for j, sc in enumerate(["S2", "S6", "S3"]):
    d = e4[e4.scenario == sc]
    ax.plot(d.n, d.robust_mean, color=C[j], marker=MK[j], label=f"{sc}: mean τ_L (Λ = {d.lam.iloc[0]:.0f})")
    ax.axhline(d.tau_true.iloc[0], color=C[j], ls="--", lw=1)
    ax.axhline(d.tau_min_pop.iloc[0], color=C[j], ls=":", lw=1)
ax.set_xscale("log"); ax.set_xlabel("cells per site"); ax.set_ylabel("τ")
ax.set_title("Robust lower bound (solid) → τ_min (dotted) ≤ τ_true (dashed)", fontsize=9); ax.legend(fontsize=7)
save(fig, "m1_fig3_bound_vs_truth")

# Fig 4: design calculus
e5 = pd.read_csv(os.path.join(T, "e5_design_required_n.csv")); rb = pd.read_csv(os.path.join(T, "e5_regional_bias.csv"))
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2))
for j, q in enumerate([0.05, 0.1, 0.3]):
    d = e5[(e5.q == q) & (e5.rho == 1)]
    axs[0].plot(d.lam, d.n_required_cells, color=C[j], marker=MK[j], label=f"metastatic frequency q = {q}")
axs[0].axhspan(1e3, 1e4, color="#999", alpha=.15, label="1k-10k cells (illustrative range)")
axs[0].set_xscale("log"); axs[0].set_yscale("log"); axs[0].set_xlabel("selection bound Λ")
axs[0].set_ylabel("primary cells needed"); axs[0].legend(fontsize=7)
axs[0].set_title("Required primary depth grows linearly in Λ / q", fontsize=9)
for j, n in enumerate(sorted(rb.n.unique())):
    d = rb[rb.n == n]; axs[1].plot(d.rho, d.false_induction, color=C[j], marker=MK[j], label=f"{n} cells, one region")
axs[1].set_xlabel("regional depletion ρ of the selected state"); axs[1].set_ylabel("P(claim induction), truth τ = 0")
axs[1].set_title("Regional sampling bias: more cells make it worse", fontsize=9); axs[1].legend(fontsize=7)
save(fig, "m1_fig4_design")

# Fig 5: founder bottleneck
sm = pd.read_csv(os.path.join(T, "e6_single_met_neutral.csv")); fl = pd.read_csv(os.path.join(T, "e6_fluctuation.csv"))
fig, axs = plt.subplots(1, 2, figsize=(8.4, 3.2))
for j, lam in enumerate(sorted(sm.lam.unique())):
    d = sm[sm.lam == lam]; axs[0].plot(d.nf, d.false_induction, color=C[j], marker=MK[j], label=f"assumed Λ = {lam:g}")
axs[0].set_xscale("log"); axs[0].set_xlabel("founders per metastasis N_f"); axs[0].set_ylabel("P(τ_L > 0), neutral seeding (τ = 0)")
axs[0].set_title("One metastasis: the mean-field bound fails at small N_f", fontsize=9); axs[0].legend(fontsize=7)
lab = {"heritable": "selection of pre-existing state", "founder_ind": "induction at founding (heritable)", "cell_ind": "per-cell induction"}
for j, m in enumerate(["heritable", "founder_ind", "cell_ind"]):
    d = fl[fl["mode"] == m].groupby("nf").var_s.median()
    axs[1].plot(d.index, d.values + 1e-6, color=C[j], marker=MK[j], ls=["-", "--", ":"][j], label=lab[m])
axs[1].set_xscale("log"); axs[1].set_yscale("log"); axs[1].set_xlabel("founders per metastasis N_f")
axs[1].set_ylabel("variance of state-s fraction\nacross metastases"); axs[1].legend(fontsize=7)
axs[1].set_title("Multi-metastasis 'fluctuation test' (8 metastases)", fontsize=9)
save(fig, "m1_fig5_founder_bottleneck")

# Fig 6: clone-resolved vs population bound
c5 = pd.read_csv(os.path.join(T, "e4_s5_clones.csv"))
g = c5.groupby("eps")[["clone_robust", "pop_robust_lam5", "pop_robust_lam5000"]].mean()
fig, ax = plt.subplots(figsize=(5.2, 3.4)); x = np.arange(len(g)); wdt = .26
for j, (col, labx) in enumerate([("clone_robust", "clone-resolved τ_L (Λ_in = 5)"),
                                 ("pop_robust_lam5", "population τ_L, Λ = 5 (misspecified)"),
                                 ("pop_robust_lam5000", "population τ_L, Λ = 5000 (valid)")]):
    ax.bar(x + (j - 1) * wdt, g[col], wdt - .02, color=C[j], label=labx, hatch=["", "//", ".."][j], edgecolor="white")
    for xi, v in zip(x + (j - 1) * wdt, g[col]):
        ax.text(xi, v + .003, f"{v:.3f}", ha="center", fontsize=6, color="#333")
ax.axhline(c5.tau_true.iloc[0], color="#333", ls="--", lw=1, label=f"τ_true = {c5.tau_true.iloc[0]:.3f}")
ax.set_xticks(x, [f"clone-label error {e:g}" for e in g.index]); ax.set_ylabel("mean lower bound")
ax.set_title("S5: clone labels turn a vacuous bound into an informative one", fontsize=9); ax.legend(fontsize=7)
save(fig, "m1_fig6_clones")

# Fig 7: dynamics identifiability
e7 = pd.read_csv(os.path.join(T, "e7_dynamics_fits.csv")); e7 = e7[e7.ncell > 0]
fig, axs = plt.subplots(1, 3, figsize=(9, 3.0))
truth = dict(a=0.08, c=0.02, db=0.6)
for ax, par in zip(axs, ["a", "c", "db"]):
    data, labs = [], []
    for des in ["2 snapshots", "4 snapshots", "8 snapshots"]:
        for nc in (2000, 20000):
            data.append(e7[(e7.design == des) & (e7.ncell == nc)][par].values); labs.append(f"{des[0]}t/{nc // 1000}k")
    ax.boxplot(data, tick_labels=labs, showfliers=False); ax.axhline(truth[par], color=C[1], ls="--", lw=1)
    ax.set_title({"a": "switch rate a (0→1)", "c": "switch rate c (1→0)", "db": "selection Δb"}[par], fontsize=9)
    ax.tick_params(axis="x", labelrotation=45)
save(fig, "m1_fig7_dynamics")
print("figures written")
