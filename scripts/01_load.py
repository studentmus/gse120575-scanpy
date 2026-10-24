"""Step 1: load the GSE120575 expression matrix + per-cell metadata into AnnData.

Input  (data/, from GEO):
  GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz
  GSE120575_patient_ID_single_cells.txt.gz
Output:
  data/gse120575_raw.h5ad   (cells x genes, X = values as published)
"""
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

DATA = Path(__file__).resolve().parents[1] / "data"
EXPR = DATA / "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
META = DATA / "GSE120575_patient_ID_single_cells.txt.gz"
OUT = DATA / "gse120575_raw.h5ad"


def read_expression(path):
    # Two header lines: cell IDs, then sample labels (e.g. "Pre_P1").
    # Data rows: gene name + one value per cell + a trailing empty field.
    with gzip.open(path, "rt", encoding="latin-1") as f:  # GEO files are not UTF-8
        cells = f.readline().rstrip("\n").split("\t")[1:]
        f.readline()
        df = pd.read_csv(f, sep="\t", header=None, index_col=0)
    df = df.iloc[:, : len(cells)]  # drop the empty trailing column
    df.columns = cells
    df.index.name = None
    return df


def read_metadata(path):
    # GEO template: ~20 lines of free text, then a table starting at "Sample name".
    with gzip.open(path, "rt", encoding="latin-1") as f:  # GEO files are not UTF-8
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("Sample name"))
    rows = [l.rstrip("\n").split("\t") for l in lines[start:]]
    rows = [r for r in rows if r[0].startswith("Sample")]  # stop at trailing notes
    header = rows[0]
    meta = pd.DataFrame([r[: len(header)] for r in rows[1:]], columns=header)
    meta = meta.rename(columns={
        "title": "cell_id",
        "characteristics: patinet ID (Pre=baseline; Post= on treatment)": "sample",
        "characteristics: response": "response",
        "characteristics: therapy": "therapy",
    })[["cell_id", "sample", "response", "therapy"]]
    # "Pre_P1" -> timepoint "Pre", patient "P1"; "Post_P1_2" is a 2nd biopsy of P1
    parts = meta["sample"].str.split("_", expand=True)
    meta["timepoint"], meta["patient"] = parts[0], parts[1]
    return meta.set_index("cell_id")


expr = read_expression(EXPR)
meta = read_metadata(META)
print(f"expression: {expr.shape[0]} genes x {expr.shape[1]} cells")
print(f"metadata:   {meta.shape[0]} cells")

missing = set(expr.columns) - set(meta.index)
assert not missing, f"{len(missing)} cells without metadata"

adata = sc.AnnData(
    X=sparse.csr_matrix(expr.T.values.astype(np.float32)),
    obs=meta.loc[expr.columns].astype("category"),
    var=pd.DataFrame(index=expr.index.astype(str)),
)
adata.var_names_make_unique()  # a few gene symbols occur twice
print(adata)
print(adata.obs.groupby(["response", "timepoint"], observed=True)
      .agg(cells=("patient", "size"), samples=("sample", "nunique")))

adata.write_h5ad(OUT, compression="gzip")
print(f"saved {OUT}")
