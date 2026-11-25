"""Stream GSE173958 10x matrices and keep only a fixed gene panel for the authors' cancer cells.

Usage: python -I extract_panel.py <raw_dir> <mendeley_proc_dir> <out_dir>
Writes out_dir/panel_M{1,2}.csv.gz: rows = colData cells; columns = metadata + log1p(1e4*count/nCount_RNA) per gene.
Streaming keeps memory small (the full matrices are ~1.6 GB gzipped).
"""
import gzip, os, sys
import numpy as np, pandas as pd

STATE_E = ["Cdh1", "Epcam", "Krt8", "Krt18", "Krt19", "Krt7", "Cldn4", "Cldn7", "Esrp1", "Ocln"]
STATE_M = ["Vim", "Zeb1", "Zeb2", "Snai2", "Twist1", "Cdh2", "Fn1", "Sparc", "Tgfbi", "Col1a1"]
HELDOUT = {
    "prolif": ["Mki67", "Top2a", "Mcm2", "Pcna", "Ccnb1", "Cdk1"],
    "glyco": ["Eno1", "Pgk1", "Ldha", "Pkm", "Tpi1", "Aldoa", "Gapdh", "Slc2a1", "Hk2"],
    "hypoxia": ["Vegfa", "Bnip3", "Ndrg1", "Car9", "Egln3", "P4ha1", "Adm"],
    "s100": ["S100a4", "S100a6", "S100a10", "S100a11"],
    "stress": ["Fos", "Jun", "Atf3", "Egr1", "Hspa1a", "Hspa1b"],
    "stem": ["Cd44", "Sox9", "Prom1", "Aldh1a1", "Lgr5"],
    "ifn": ["Isg15", "Ifit1", "Ifit3", "Irf7", "Stat1"],
    "adhesion": ["Lamb3", "Lamc2", "Itgb4", "Itga6", "Col17a1"],
    "chromatin": ["Hmga2", "Hmga1", "Ezh2", "Kdm5b"],
}
PANEL = STATE_E + STATE_M + [g for v in HELDOUT.values() for g in v]
SAMPLES = {1: {"PT": "GSM5283482_M1-PT", "PTab": "GSM5283483_M1-SS", "Met": "GSM5283484_M1-Met",
               "Liver": "GSM5283485_M1-Liver", "Lung": "GSM5283486_M1-Lung", "Blood": "GSM5283487_M1-CTC"},
           2: {"PT": "GSM5283488_M2-PT", "Met": "GSM5283489_M2-Met", "Liver": "GSM5283490_M2-Liver",
               "Lung": "GSM5283491_M2-Lung"}}


def extract(raw, prefix, keep_cells):
    feats = pd.read_csv(os.path.join(raw, prefix + "-features.tsv.gz"), sep="\t", header=None)
    bcs = pd.read_csv(os.path.join(raw, prefix + "-barcodes.tsv.gz"), header=None)[0].str.replace("-1", "", regex=False)
    gidx = {i + 1: g for i, g in enumerate(feats[1]) if g in set(PANEL)}
    cidx = {j + 1: b for j, b in enumerate(bcs) if b in keep_cells}
    out = np.zeros((len(cidx), len(PANEL)), np.float32)
    crow = {j: k for k, j in enumerate(cidx)}; gcol = {i: PANEL.index(g) for i, g in gidx.items()}
    with gzip.open(os.path.join(raw, prefix + "-matrix.mtx.gz"), "rt") as fh:
        for line in fh:
            if not line.startswith("%"):
                break
        for chunk in pd.read_csv(fh, sep=" ", header=None, names=["g", "c", "v"], chunksize=5_000_000):
            m = chunk.g.isin(gidx.keys()) & chunk.c.isin(cidx.keys())
            for g, c, v in chunk[m].itertuples(index=False):
                out[crow[c], gcol[g]] += v
    return pd.DataFrame(out, columns=PANEL, index=[cidx[j] for j in cidx])


def main(raw, proc, outd):
    os.makedirs(outd, exist_ok=True)
    for m, samples in SAMPLES.items():
        cd = pd.read_csv(os.path.join(proc, f"M{m}/input-files/mouse{m}_cds_colData.txt"), sep="\t")
        cd["prefix"] = cd.rowname.str.split("_").str[1]; cd["bc"] = cd.rowname.str.split("_").str[2]
        parts = []
        for pre, gsm in samples.items():
            sub = cd[cd.prefix == pre]
            X = extract(raw, gsm, set(sub.bc))
            X["rowname"] = f"M{m}_{pre}_" + X.index
            parts.append(X)
            print(m, pre, len(sub), "matched", len(X), flush=True)
        X = pd.concat(parts).set_index("rowname")
        d = cd.set_index("rowname").join(X, how="inner")
        lib = d["nCount_RNA"].values[:, None]
        d[PANEL] = np.log1p(1e4 * d[PANEL].values / lib)
        d.to_csv(os.path.join(outd, f"panel_M{m}.csv.gz"))
        missing = [g for g in PANEL if (d[g] > 0).sum() == 0]
        print(f"M{m}: {len(d)} cells; genes never detected: {missing}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:4])
