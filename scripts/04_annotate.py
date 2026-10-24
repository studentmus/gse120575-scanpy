"""Step 4: annotate Leiden clusters with canonical immune markers.

Labels are assigned by hand from the marker table / dotplot below (the
CLUSTER_LABELS dict), i.e. a transparent, reviewable decision.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "data" / "gse120575_clustered.h5ad"
OUT = ROOT / "data" / "gse120575_annotated.h5ad"
sc.settings.figdir = ROOT / "figures"

MARKERS = {
    "T": ["CD3E", "CD2"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4", "IL7R"],
    "Treg": ["FOXP3", "IL2RA"],
    "memory/naive-like": ["TCF7", "SELL", "CCR7"],
    "exhausted": ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "ENTPD1"],
    "cycling": ["MKI67", "TOP2A"],
    "NK": ["NCAM1", "KLRD1", "GNLY"],
    "B": ["CD19", "MS4A1"],
    "plasma": ["MZB1", "IGJ"],  # IGJ = JCHAIN (old symbol)
    "mono/macro": ["CD14", "CD68", "LYZ"],
    "DC": ["FCER1A", "CLEC10A", "LILRA4"],
}

# Assigned by hand from the printed marker table (mean log2(TPM+1) per cluster).
# T-cell rules: HAVCR2/PDCD1 >= ~3.5 -> exhausted; TCF7 >= 3.5 with IL7R -> memory-like;
# T clusters with no clear profile (CD8 ~5-6, CD4 ~2, moderate PD1) stay "T mixed".
CLUSTER_LABELS = {
    "0": "NK",                       # GNLY 9.5, KLRD1 6.4, CD3E low
    "1": "T mixed",
    "2": "T memory-like (TCF7+)",    # CD8A 6.2, TCF7 3.5, IL7R 5.0, HAVCR2 0.6
    "3": "T memory-like (TCF7+)",    # TCF7 3.7, IL7R 5.9, CCR7 3.4
    "4": "CD8 T cytotoxic",          # CD8A 7.1, KLRD1 4.7, HAVCR2 2.2
    "5": "T mixed",
    "6": "T mixed",
    "7": "Treg",                     # CD4 6.2, FOXP3 3.0, IL2RA 3.7
    "8": "B",                        # CD19 5.8, MS4A1 7.2
    "9": "T mixed",
    "10": "CD8 T exhausted",         # PDCD1 5.5, HAVCR2 6.1, ENTPD1 3.5, TCF7 0.7
    "11": "T mixed",
    "12": "pDC",                     # LILRA4 10.0, IGJ 11.5
    "13": "CD8 T exhausted",         # HAVCR2 3.7, PDCD1 3.7, LAG3 3.2
    "14": "Unclear (T+myeloid)",     # CD3E and CD68 both ~3-4: likely doublets, 83 cells
    "15": "T mixed",
    "16": "T cycling",               # MKI67 5.9, TOP2A 6.8
    "17": "Mono/macrophage",         # LYZ 10.8, CD68 7.7, CD14 6.6
    "18": "Mono/macrophage",         # CD14 11.0, CD68 9.6
    "19": "T mixed",
    "20": "T mixed",
    "21": "B",                       # MS4A1 9.1
    "22": "Plasma",                  # MZB1 9.6, IGJ 9.8
    "23": "B",                       # CD19 6.0, CCR7 7.2 (naive-like B)
}

adata = sc.read_h5ad(IN)
genes = [g for gs in MARKERS.values() for g in gs]

# Mean log2(TPM+1) per cluster, from the unscaled values kept in .raw
expr = pd.DataFrame(adata.raw[:, genes].X.toarray(), columns=genes, index=adata.obs_names)
table = expr.groupby(adata.obs["leiden"], observed=True).mean().round(1)
pd.set_option("display.width", 250)
print(table.to_string())

sc.pl.dotplot(adata, MARKERS, groupby="leiden", use_raw=True,
              standard_scale="var", save="_markers.png", show=False)

if CLUSTER_LABELS:
    adata.obs["cell_type"] = adata.obs["leiden"].map(CLUSTER_LABELS).astype("category")
    print(adata.obs["cell_type"].value_counts().to_string())
    sc.pl.umap(adata, color="cell_type", save="_cell_types.png", show=False)
    adata.write_h5ad(OUT, compression="gzip")
    print(f"saved {OUT}")
