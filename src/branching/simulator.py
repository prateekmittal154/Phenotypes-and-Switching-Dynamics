"""Lineage-aware stochastic birth-death-switching-dissemination simulator (2 states, 1 primary + 2 metastatic organs).

Cells are individuals. Each step of length dt, every live cell independently does one of:
- divide (prob beta[o,s] dt),
- die (delta[o,s] dt),
- switch state (rate a[c,o] for 0->1, b[c,o] for 1->0),
- if in the primary, disseminate (mu * kappa[s] dt).

Events are drawn from one uniform, so dt * total rate must stay below 1. This is the standard
discrete-time approximation of the continuous-time multitype branching process. Its error is O(dt) per event, and
dt = 0.02 with rates <= 2/day keeps every per-step probability below 0.1.

A disseminated cell picks a metastatic organ o with probability q[o] and founds with probability phi[o], otherwise it
dies. This gives stochastic founding on top of the branching extinction of small founder lineages.

Ancestry: every division creates two daughter nodes (parent pointer, birth time, organ). Evolving-barcode edits occur
only after t_dox: each daughter is edited with probability p_edit. That gives the 'visible' internal nodes from which
observed trees are built, mimicking macsGESTALT, whose trees resolve only branches carrying an edit.

Ground truth returned:
- the full node table (parent, time, edited);
- the sampled cells with clone, organ and true state;
- the founder records (time, clone, state at dissemination, organ);
- per-clone parameters;
- switch-event counts.
"""
from dataclasses import dataclass, field
import numpy as np

S, O = 2, 3  # states (0 low-EMT, 1 high-EMT); organs (0 primary, 1 met site A, 2 met site B)


@dataclass
class Params:
    C: int = 30                      # clones (single-cell founders labelled at t = 0)
    rho: float = 0.3                 # founder P(state 1) at labelling
    beta: np.ndarray = field(default_factory=lambda: np.full((O, S), 0.9))
    delta: np.ndarray = field(default_factory=lambda: np.full((O, S), 0.2))
    a: float = 0.3                   # base switching rate 0 -> 1
    b: float = 0.7                   # base switching rate 1 -> 0
    eta_sd: float = 0.0              # clone-specific log-propensity SD (M1): a_c = a e^eta, b_c = b e^-eta
    organ_gain: np.ndarray = field(default_factory=lambda: np.ones(O))  # multiplies a in organ o (M2 induction)
    organ_loss: np.ndarray = field(default_factory=lambda: np.ones(O))  # multiplies b in organ o
    mu: float = 0.004                # dissemination rate per primary cell per day
    kappa: np.ndarray = field(default_factory=lambda: np.array([1.0, 1.0]))  # state competence for dissemination
    q: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.5, 0.5]))
    phi: np.ndarray = field(default_factory=lambda: np.array([0.0, 0.3, 0.3]))  # founding success
    T: float = 11.0
    dt: float = 0.02
    t_dox: float = 4.0
    p_edit: float = 0.35
    max_cells: int = 250_000
    eta: np.ndarray = None           # realised clone propensities (filled in)


def simulate(p: Params, rng):
    eta = rng.normal(0, p.eta_sd, p.C) if p.eta_sd > 0 else np.zeros(p.C)
    p.eta = eta
    # live-cell arrays
    node = np.arange(p.C); clone = np.arange(p.C); organ = np.zeros(p.C, int)
    state = (rng.random(p.C) < p.rho).astype(int)
    N_par, N_t, N_edit, N_org = [-1] * p.C, [0.0] * p.C, [False] * p.C, [0] * p.C
    founders = []; n_switch = np.zeros(O, int)
    t = 0.0
    a_co = p.a * np.exp(eta)[:, None] * p.organ_gain[None, :]
    b_co = p.b * np.exp(-eta)[:, None] * p.organ_loss[None, :]
    while t < p.T - 1e-9 and len(node) < p.max_cells:
        br = p.beta[organ, state]; dr = p.delta[organ, state]
        sr = np.where(state == 0, a_co[clone, organ], b_co[clone, organ])
        mr = np.where(organ == 0, p.mu * p.kappa[state], 0.0)
        u = rng.random(len(node)) / p.dt
        div = u < br
        die = (~div) & (u < br + dr)
        sw = (~div) & (~die) & (u < br + dr + sr)
        mig = (~div) & (~die) & (~sw) & (u < br + dr + sr + mr)
        # switching
        np.add.at(n_switch, organ[sw], 1)
        state = np.where(sw, 1 - state, state)
        # dissemination and founding
        if mig.any():
            idx = np.where(mig)[0]
            dest = rng.choice(O, len(idx), p=p.q)
            ok = rng.random(len(idx)) < p.phi[dest]
            for i, o, k in zip(idx, dest, ok):
                if k:
                    founders.append((t, int(clone[i]), int(state[i]), int(o), int(node[i])))
            organ = organ.copy(); organ[idx[ok]] = dest[ok]
            die = die.copy(); die[idx[~ok]] = True
        # divisions: two daughter nodes
        keep = ~die & ~div
        d_idx = np.where(div)[0]
        nd = len(d_idx)
        if nd:
            start = len(N_par)
            par = np.repeat(node[d_idx], 2)
            edited = (rng.random(2 * nd) < p.p_edit) & (t >= p.t_dox)
            N_par.extend(par.tolist()); N_t.extend([t] * (2 * nd)); N_edit.extend(edited.tolist())
            N_org.extend(np.repeat(organ[d_idx], 2).tolist())
            new_nodes = np.arange(start, start + 2 * nd)
            node = np.concatenate([node[keep], new_nodes])
            clone = np.concatenate([clone[keep], np.repeat(clone[d_idx], 2)])
            organ = np.concatenate([organ[keep], np.repeat(organ[d_idx], 2)])
            state = np.concatenate([state[keep], np.repeat(state[d_idx], 2)])
        else:
            node, clone, organ, state = node[keep], clone[keep], organ[keep], state[keep]
        t += p.dt
    nodes = dict(parent=np.array(N_par), t=np.array(N_t), edited=np.array(N_edit), organ=np.array(N_org))
    return dict(node=node, clone=clone, organ=organ, state=state, nodes=nodes, founders=founders,
                n_switch=n_switch, eta=eta, t_end=t)


def sample_cells(sim, n_per_organ, rng, eps=0.0, clone_loss=0.0):
    """Partial sampling per organ, expression misclassification eps and lost clone labels (clone_loss)."""
    keep = []
    for o in range(O):
        idx = np.where(sim["organ"] == o)[0]
        keep.append(rng.choice(idx, min(n_per_organ[o], len(idx)), replace=False) if len(idx) else idx)
    k = np.concatenate(keep)
    k = k[rng.random(len(k)) >= clone_loss]
    obs = sim["state"][k].copy()
    flip = rng.random(len(k)) < eps
    obs[flip] = 1 - obs[flip]
    return dict(node=sim["node"][k], clone=sim["clone"][k], organ=sim["organ"][k], state_true=sim["state"][k],
                state=obs)


def edit_paths(sample, nodes):
    """Observed lineage path for each sampled cell: the tuple of edited ancestor nodes, root to leaf.

    This mirrors the evolving-barcode resolution: only edited branches are visible.
    """
    par, ed = nodes["parent"], nodes["edited"]
    paths = []
    for n in sample["node"]:
        pth = []
        while n >= 0:
            if ed[n]:
                pth.append(int(n))
            n = par[n]
        paths.append(tuple(reversed(pth)))
    return paths


def true_paths(sample, nodes):
    """Full (unobservable) ancestor path, root to leaf, for exact-coalescence checks."""
    par = nodes["parent"]; out = []
    for n in sample["node"]:
        pth = []
        while n >= 0:
            pth.append(int(n)); n = par[n]
        out.append(tuple(reversed(pth)))
    return out
