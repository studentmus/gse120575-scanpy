"""Step 2: QC, following Sade-Feldman et al. 2018 STAR Methods.

Paper (values are log2(TPM+1)):
  cells: drop if < 1000 detected genes, or mean housekeeping expression < 2.5,
         or zero expression of both PTPRC (CD45) and CD3E.
  genes: keep genes with log2(TPM+1) > 4.5 in >= 10 cells (paper: per sample,
         protein-coding only - simplified here, see docs/02_qc.md).
The GEO matrix already contains the 16,291 post-QC cells, so the cell filters
are a check: we expect (almost) nothing to be removed.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "data" / "gse120575_raw.h5ad"
OUT = ROOT / "data" / "gse120575_qc.h5ad"
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# Approximation: the paper does not list its housekeeping genes in the Methods.
HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "PPIA", "RPLP0", "TPT1", "EEF1A1", "TMSB4X"]

adata = sc.read_h5ad(IN)
X = adata.X

adata.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
adata.obs["housekeeping"] = np.asarray(adata[:, HOUSEKEEPING].X.mean(axis=1)).ravel()
cd45_cd3e = adata[:, ["PTPRC", "CD3E"]].X.toarray()
adata.obs["cd45_or_cd3e"] = (cd45_cd3e > 0).any(axis=1)

print(adata.obs[["n_genes", "housekeeping"]].describe().round(2))

fails = {
    "n_genes < 1000": adata.obs["n_genes"] < 1000,
    "housekeeping < 2.5": adata.obs["housekeeping"] < 2.5,
    "PTPRC = CD3E = 0": ~adata.obs["cd45_or_cd3e"],
}
for name, mask in fails.items():
    print(f"{name:20s} {mask.sum():6d} cells")
keep = ~np.logical_or.reduce(list(fails.values()))
print(f"cells kept: {keep.sum()} / {adata.n_obs}")

fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
axes[0].hist(adata.obs["n_genes"], bins=80, color="grey")
axes[0].axvline(1000, color="red", ls="--")
axes[0].set(xlabel="detected genes per cell", ylabel="cells")
axes[1].hist(adata.obs["housekeeping"], bins=80, color="grey")
axes[1].axvline(2.5, color="red", ls="--")
axes[1].set(xlabel="mean housekeeping log2(TPM+1)")
fig.tight_layout()
fig.savefig(FIG / "qc_histograms.png", dpi=150)

adata = adata[keep].copy()

# Gene filter: expressed at a meaningful level (> 4.5, i.e. TPM > ~22) in >= 10 cells.
expressed = np.asarray((adata.X > 4.5).sum(axis=0)).ravel()
adata = adata[:, expressed >= 10].copy()
print(f"genes kept: {adata.n_vars}")

adata.write_h5ad(OUT, compression="gzip")
print(f"saved {OUT}")
