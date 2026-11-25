"""E12b: recompute real-recorder windows on an INFORMATION basis, not the ad-hoc M3 criterion.

Milestone 3 defined the usable range of k as 1/(2 L_max) .. 1/(2 L_min), which is the range over which rho is
neither saturated nor flat. That is a hard identifiability statement but it is arbitrary as a practical
criterion. E12 shows the Fisher information retains >=50% of its peak over k*L in roughly [0.10, 1.0], a ~10-fold
band. The efficiency window for a recorder whose within-clone pairs span lineage distances [L_min, L_max] is
therefore k in [0.10 / L_max, 1.0 / L_min], which is about one decade wider than the M3 figure.

Both are reported. Uses the depth spectra already computed in Milestone 3.
"""
import os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAB = os.path.join(ROOT, "results/tables")
LO_KL, HI_KL = 0.10, 1.00          # from e12_optimum_width.csv, the >=50%-of-peak band

d = pd.read_csv(os.path.join(TAB, "m3_depth_spectra_summary.csv"))
d["L_max"] = (d.depth_max - d.depth_p10).clip(lower=0.5)
d["L_min"] = (d.depth_max - d.depth_p90).clip(lower=0.5)
d["k_lo_M3"] = 1 / (2 * d.L_max)
d["k_hi_M3"] = 1 / (2 * d.L_min)
d["k_lo_info"] = LO_KL / d.L_max
d["k_hi_info"] = HI_KL / d.L_min
d["decades_M3"] = np.log10(d.k_hi_M3 / d.k_lo_M3)
d["decades_info"] = np.log10(d.k_hi_info / d.k_lo_info)
d.to_csv(os.path.join(TAB, "e12b_windows_recomputed.csv"), index=False)
g = d.groupby(["system", "dataset"])[["depth_max", "k_lo_M3", "k_hi_M3", "decades_M3",
                                      "k_lo_info", "k_hi_info", "decades_info"]].median()
print(g.round(3).to_string())
print("\nmedian width, M3 criterion:   %.2f decades" % d.decades_M3.median())
print("median width, info criterion: %.2f decades" % d.decades_info.median())
