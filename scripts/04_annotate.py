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
# T-cell rules: CD8A high with HAVCR2 and PDCD1 >= ~4 -> exhausted; TCF7 >= 3.5 with
# IL7R and low HAVCR2 -> memory-like; MKI67/TOP2A high -> cycling. T clusters with no
# clear profile (CD8 ~5-6, CD4 ~2, moderate PD1, or CD4-CD8- like 11) stay "T mixed".
CLUSTER_LABELS = {
    "0": "NK",                       # GNLY 10.3, KLRD1 6.5, CD3E lower than T clusters
    "1": "CD8 T exhausted",          # CD8A 8.8, PDCD1 5.3, HAVCR2 4.9, TCF7 0.9
    "2": "T mixed",
    "3": "T memory-like (TCF7+)",    # TCF7 3.6, IL7R 5.9, CCR7 3.1, HAVCR2 0.8
    "4": "Treg",                     # CD4 6.6, FOXP3 3.3, IL2RA 4.2
    "5": "T cycling",                # MKI67 5.3, TOP2A 6.2 (mostly CD8)
    "6": "T mixed",
    "7": "CD8 T cytotoxic",          # CD8A 8.2, KLRD1 3.2, HAVCR2 2.8, TCF7 1.6
    "8": "T mixed",
    "9": "T mixed",
    "10": "pDC",                     # LILRA4 10.1, IGJ 11.9
    "11": "T mixed",                 # CD3E+, CD8A/CD4 ~0, GNLY 4.6: CD4-CD8- T
    "12": "CD8 T exhausted",         # CD8A 7.6, HAVCR2 4.2, PDCD1 4.0, LAG3 4.0
    "13": "Mono/macrophage",         # CD14 9.7, CD68 8.9
    "14": "Mono/macrophage",         # LYZ 10.8, CD68 7.7, CD14 6.7
    "15": "B",                       # CD19 6.7, MS4A1 9.5
    "16": "Plasma",                  # MZB1 9.4, IGJ 10.8
    "17": "B",                       # CD19 6.2, CCR7 6.5 (naive-like B)
    "18": "T mixed",
    "19": "T memory-like (TCF7+)",   # CD8A 6.7, TCF7 3.5, IL7R 4.8, HAVCR2 0.7
    "20": "T mixed",
    "21": "T cycling",               # MKI67 5.7, TOP2A 6.8
    "22": "B",                       # CD19 6.1, MS4A1 7.7
    "23": "T memory-like (TCF7+)",   # TCF7 4.7, IL7R 3.9, SELL 3.3
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
