"""Step 3: HVG -> PCA -> neighbours -> UMAP -> Leiden.

X is already log2(TPM+1) (see docs/01_load.md), so no normalize_total / log1p.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import scanpy as sc

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "data" / "gse120575_qc.h5ad"
OUT = ROOT / "data" / "gse120575_clustered.h5ad"
sc.settings.figdir = ROOT / "figures"

SEED = 0
adata = sc.read_h5ad(IN)

# 1. Highly variable genes. flavor="seurat" works on log data (seurat_v3 needs counts).
# It un-logs with expm1 internally; declaring base 2 makes it invert log2 correctly.
adata.uns["log1p"] = {"base": 2}
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
print(f"HVGs: {adata.var['highly_variable'].sum()}")

# 2. Keep the unscaled log values for plotting markers; scale a copy for PCA.
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, random_state=SEED)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="_variance.png", show=False)

# 3. kNN graph on the first 30 PCs, then UMAP (for viewing) and Leiden (for clusters).
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, random_state=SEED)
print(adata.obs["leiden"].value_counts().sort_index().to_string())

sc.pl.umap(adata, color=["leiden", "patient", "response", "timepoint"],
           ncols=2, wspace=0.4, save="_overview.png", show=False)

adata.write_h5ad(OUT, compression="gzip")
print(f"saved {OUT}")
