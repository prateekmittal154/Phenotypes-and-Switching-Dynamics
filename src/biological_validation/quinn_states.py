"""Derive cell states for the Quinn (GSE161363) A549 xenograft system, independently of EMT and of lineage.

The Milestone 4 brief requires state definitions derived within this system rather than imported from the
Simeonov EMT analysis, to avoid circular validation. The procedure is:

  pass 1  stream the matrix market file, accumulating per-gene detection counts and summed counts;
  pass 2  re-stream, extracting the top `n_hvg` genes by detection breadth into a dense cells x genes array;
  then    log-normalise per cell, standardise per gene, regress out log library size, and take the first
          principal component of the PRIMARY-SITE cells only as the state axis.

The axis is therefore unsupervised and fitted on primary-site cells, so it cannot be tuned to any
metastasis-versus-primary contrast. Binarisation is at the primary-site median. What PC1 correlates with is
reported rather than assumed.

Usage: python -I quinn_states.py <quinn_dir> <out_dir>
"""
import gzip, os, sys
import numpy as np, pandas as pd

MICE = {"5k": ("GSM4905335", "GSM4905334_trees_5k"),
        "30k": ("GSM4905337", "GSM4905336_trees_30k"),
        "10k": ("GSM4905339", "GSM4905338_trees_10k"),
        "100k": ("GSM4905341", "GSM4905340_trees_100k")}
# Primary implantation site is the lung; mediastinal nodes, liver and blood are disseminated sites.
PRIMARY_PREFIX = ("LL", "RL", "RE", "RW", "NE", "NW", "SE", "SW", "30kLL", "30kRL")
N_HVG = 2000
# Mitochondrial and ribosomal genes dominate technical variance in these libraries and are excluded from the
# state axis. A cell-cycle score is computed separately as a second, biologically interpretable state
# definition, so that conclusions do not rest on one unsupervised axis.
TECH_PREFIX = ("MT-", "MTRNR", "RPL", "RPS")
S_GENES = ["MCM2", "PCNA", "TYMS", "RRM2", "MCM6", "CDC45", "UNG", "GINS2", "CDC6", "RFC2"]
G2M_GENES = ["TOP2A", "MKI67", "CCNB2", "CDK1", "BUB1", "AURKA", "PLK1", "TPX2", "CENPF", "UBE2C"]


def _open(p):
    return gzip.open(p, "rt")


def pass1(mtx, n_genes):
    det = np.zeros(n_genes, np.int32)
    with _open(mtx) as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=8_000_000):
            np.add.at(det, chunk.g.values - 1, 1)
    return det


def pass2(mtx, keep_genes, n_cells):
    col = {g: i for i, g in enumerate(keep_genes)}
    X = np.zeros((n_cells, len(keep_genes)), np.float32)
    with _open(mtx) as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=8_000_000):
            m = chunk.g.isin(col)
            sub = chunk[m]
            if len(sub):
                X[sub.c.values - 1, [col[g] for g in sub.g.values]] += sub.v.values
    return X


def lib_sizes(mtx, n_cells):
    tot = np.zeros(n_cells, np.float64)
    with _open(mtx) as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=8_000_000):
            np.add.at(tot, chunk.c.values - 1, chunk.v.values)
    return tot


def run_mouse(qdir, mouse, out):
    gsm, _ = MICE[mouse]
    bcs = pd.read_csv(os.path.join(qdir, f"{gsm}_barcodes.{mouse}.tsv.gz"), header=None)[0].values
    genes = pd.read_csv(os.path.join(qdir, f"{gsm}_genes.{mouse}.tsv.gz"), sep="\t", header=None)
    meta = pd.read_csv(os.path.join(qdir, f"{gsm}_meta.{mouse}.tsv.gz"), sep="\t")
    mtx = os.path.join(qdir, f"{gsm}_matrix.{mouse}.mtx.gz")
    with _open(mtx) as fh:
        for line in fh:
            if not line.startswith("%"):
                ng, nc, _ = map(int, line.split()); break
    assert nc == len(bcs) and ng == len(genes), (ng, nc, len(genes), len(bcs))
    det = pass1(mtx, ng)
    sym_all = genes.iloc[:, 1].astype(str).values
    tech = np.array([s0.startswith(TECH_PREFIX) for s0 in sym_all])
    det_f = det.astype(float).copy(); det_f[tech] = -1
    keep = np.argsort(-det_f)[:N_HVG] + 1
    cc_idx = np.where(np.isin(sym_all, S_GENES + G2M_GENES))[0] + 1
    mt_idx = np.where(np.array([str(x).startswith("MT-") for x in sym_all]))[0] + 1
    extra = sorted(set(cc_idx.tolist()) | set(mt_idx.tolist()))
    X = pass2(mtx, set(keep.tolist()) | set(extra), nc)
    cols = {g: i for i, g in enumerate(sorted(set(keep.tolist()) | set(extra)))}
    tot = lib_sizes(mtx, nc)
    keep_pos = [cols[g] for g in keep]
    sym = genes.iloc[keep - 1, 1].values
    # barcode naming differs per mouse: 5k/30k/100k match the matrix barcodes directly, 10k needs the
    # "<site>." prefix stripped. Pick whichever convention matches the most cells.
    import re as _re
    plain = np.array([_re.sub(r"^[^.]*\.", "", b) for b in bcs])
    mset = set(meta.cellBC)
    n_direct = sum(b in mset for b in bcs); n_strip = sum(b in mset for b in plain)
    key = bcs if n_direct >= n_strip else plain
    d = pd.DataFrame(dict(cellBC=key, _row=np.arange(nc), lib=tot))
    d = d.merge(meta, on="cellBC", how="inner").drop_duplicates("cellBC")
    idx = d._row.values
    Xi_all = X[idx]
    Xi = Xi_all[:, keep_pos]
    d = d.drop(columns=["_row"])
    L = np.maximum(d.lib.values, 1)[:, None]
    Y = np.log1p(1e4 * Xi / L)
    Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-8)
    loglib = np.log(d.lib.values)
    loglib = (loglib - loglib.mean()) / (loglib.std() + 1e-9)
    Y = Y - np.outer(loglib, (loglib @ Y) / (loglib @ loglib))      # regress out library size
    prim = d.sampleID.astype(str).str.replace("^30k", "", regex=True).isin(PRIMARY_PREFIX).values
    if prim.sum() < 50:
        prim = np.ones(len(d), bool)
    U, S, Vt = np.linalg.svd(Y[prim] - Y[prim].mean(0), full_matrices=False)
    pc1 = Y @ Vt[0]
    d["pc1"] = pc1
    d["is_primary"] = prim
    d["state2"] = (pc1 > np.median(pc1[prim])).astype(int)
    # technical diagnostics and a second, cell-cycle state definition
    L1 = np.maximum(d.lib.values, 1)
    mt_pos = [cols[g] for g in mt_idx if g in cols]
    cc_pos = [cols[g] for g in cc_idx if g in cols]
    d["mt_frac"] = Xi_all[:, mt_pos].sum(1) / L1 if mt_pos else 0.0
    d["n_genes"] = (Xi_all > 0).sum(1)
    if cc_pos:
        cc = np.log1p(1e4 * Xi_all[:, cc_pos] / L1[:, None])
        cc = (cc - cc.mean(0)) / (cc.std(0) + 1e-8)
        d["cc_score"] = cc.mean(1)
        d["state2_cc"] = (d.cc_score > np.median(d.cc_score[prim])).astype(int)
    load = pd.Series(Vt[0], index=sym).sort_values()
    d.to_csv(os.path.join(out, f"quinn_cells_{mouse}.csv.gz"), index=False)
    pd.DataFrame(dict(gene=load.index, loading=load.values)).to_csv(
        os.path.join(out, f"quinn_pc1_loadings_{mouse}.csv"), index=False)
    var_expl = float(S[0] ** 2 / (S ** 2).sum())
    print(f"{mouse}: {len(d)} cells with lineage, {d.LineageGroup.nunique()} lineage groups, "
          f"{int(prim.sum())} primary-site cells, PC1 var {var_expl:.3f}", flush=True)
    print("   PC1 top negative:", ", ".join(load.index[:8]), flush=True)
    print("   PC1 top positive:", ", ".join(load.index[-8:]), flush=True)
    return dict(mouse=mouse, n_cells=len(d), n_lineage_groups=int(d.LineageGroup.nunique()),
                n_primary=int(prim.sum()), pc1_var=var_expl,
                sites=int(d.sampleID.nunique()),
                corr_pc1_loglib=float(np.corrcoef(pc1, np.log(d.lib))[0, 1]),
                corr_pc1_mt=float(np.corrcoef(pc1, d.mt_frac)[0, 1]),
                corr_pc1_ngenes=float(np.corrcoef(pc1, d.n_genes)[0, 1]),
                corr_pc1_cc=float(np.corrcoef(pc1, d.cc_score)[0, 1]) if "cc_score" in d else np.nan)


def main(qdir, out):
    os.makedirs(out, exist_ok=True)
    rows = [run_mouse(qdir, m, out) for m in MICE]
    pd.DataFrame(rows).to_csv(os.path.join(out, "quinn_state_summary.csv"), index=False)
    print("done")


if __name__ == "__main__":
    main(*sys.argv[1:3])
