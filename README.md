# Immune cell composition and checkpoint response in melanoma (GSE120575)

A minimal, reproducible Scanpy re-analysis of
**Sade-Feldman et al., *Cell* 2018** — "Defining T cell states associated with response to
checkpoint immunotherapy in melanoma" ([PMC6641984](https://pmc.ncbi.nlm.nih.gov/articles/PMC6641984/)).

## Question
Do tumours from patients who respond to checkpoint blockade (anti-PD1 / anti-CTLA4) contain a
different mix of immune cell types than tumours from non-responders, before and during therapy?

## Data
[GEO GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575): Smart-seq2
scRNA-seq of CD45+ cells sorted from 48 melanoma biopsies (32 patients), taken at baseline
(Pre) or on treatment (Post), labelled Responder / Non-responder per biopsy.

| | Pre | Post |
|---|---|---|
| Responder | 9 samples, 2,725 cells | 8 samples, 2,839 cells |
| Non-responder | 10 samples, 3,203 cells | 21 samples, 7,524 cells |

**The GEO matrix is log2(TPM+1), not raw counts.** So:
- there is no `normalize_total` / `log1p`;
- the log base (2) is declared to Scanpy before HVG selection;
- count-based methods (HVG `seurat_v3`, scVI, pseudobulk DESeq2) don't apply.

## Pipeline
| Step | Script | What it does |
|---|---|---|
| 1 | [`01_load.py`](scripts/01_load.py) | GEO matrix + metadata → AnnData |
| 2 | [`02_qc.py`](scripts/02_qc.py) | Paper's QC: ≥1000 genes, housekeeping ≥2.5, PTPRC or CD3E > 0; genes > 4.5 in ≥10 cells → 16,190 cells, 27,686 genes |
| 3 | [`03_cluster.py`](scripts/03_cluster.py) | 2,000 HVGs (`seurat` flavor, log base 2) → scale → 50 PCs → kNN (15 neighbours, 30 PCs) → UMAP + Leiden (res 1.0, 24 clusters) |
| 4 | [`04_annotate.py`](scripts/04_annotate.py) | 24 clusters → 11 cell types by canonical markers (labels and evidence in code) |
| 5 | [`05_proportions.py`](scripts/05_proportions.py) | Cell-type fraction per sample, R vs NR by timepoint, Mann–Whitney per type |

## Result
![Cell-type proportions per sample](figures/celltype_proportions.png)

Each dot is one biopsy. Directions agree with the paper:

- **Exhausted CD8 T cells:** higher in non-responders (median Pre 22% vs 11%, p = 0.025; Post 18% vs 6%, p = 0.032).
- **Cycling T cells:** higher in non-responders (Pre 5.4% vs 0.9%, p = 0.002). This is the only difference that survives Bonferroni correction across 22 tests.
- **TCF7+ memory-like T cells:** higher in responders before therapy (26% vs 12%, p = 0.006, not significant after correction).
- **B cells:** responders tend to have more (p = 0.06 Pre, 0.07 Post; not significant at this n).

p-values are unadjusted two-sided Mann–Whitney tests unless stated. All values are in [`results/`](results/).

![UMAP of annotated cell types](figures/umap_cell_types.png)

## Limitations
- **Small n:** 8–21 samples per group, 22 tests, so the results are descriptive.
- **Non-independence:** some patients contribute several biopsies, and response is labelled per biopsy, so one patient can carry both labels. A mixed model with a patient random effect would be stricter.
- **Mixed therapies:** anti-PD1, anti-CTLA4 and the combination are pooled.
- **Annotation is cluster-level and manual.** 4,589 cells stay "T mixed" (CD4/CD8 not resolved). T-cell subclustering, as in the paper, would split CD8 states further.
- **QC is approximate:**
  - the paper's housekeeping gene list isn't published;
  - its gene filter was per-sample and protein-coding only;
  - 101 cells (0.6%) of the published matrix fail the literal rules and were removed.
- **Compositional data:** fractions within a sample are dependent. scCODA or propeller would model this properly.
- **No batch correction:** T cells mix well across patients; some myeloid clusters are patient-specific.

## Reproduce
Data files are not in the repo. Download them into `data/`:
```bash
mkdir data && cd data
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz
```
Create the pinned environment and run all steps:
```bash
conda env create -f environment.yml
conda activate scanpy
for s in 01_load 02_qc 03_cluster 04_annotate 05_proportions; do python scripts/$s.py; done
```
- On Windows, use `.\run_all.ps1` (it runs each step via `conda run`).
- Steps 3–4 use `random_state=0`.

## License
[MIT](LICENSE)
