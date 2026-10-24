"""Step 5: cell-type proportions per sample, responders vs non-responders, by timepoint.

The unit is the sample (biopsy), not the cell: proportions are computed within
each sample, so samples with many cells don't dominate.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "data" / "gse120575_annotated.h5ad"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

COLORS = {"Responder": "#2a78d6", "Non-responder": "#eb6834"}
INK, MUTED = "#1f1f1e", "#6b6a64"

obs = sc.read_h5ad(IN).obs

# Cells per (sample, cell type) -> fraction within each sample
counts = pd.crosstab(obs["sample"], obs["cell_type"])
props = counts.div(counts.sum(axis=1), axis=0)
info = obs.groupby("sample", observed=True)[["response", "timepoint", "patient"]].first()
long = (props.join(info).reset_index()
        .melt(id_vars=["sample", "response", "timepoint", "patient"],
              var_name="cell_type", value_name="fraction"))
long.to_csv(RES / "celltype_proportions_per_sample.csv", index=False)

# Descriptive Mann-Whitney R vs NR per cell type and timepoint (unadjusted)
tests = []
for (ct, tp), d in long.groupby(["cell_type", "timepoint"]):
    r = d.loc[d["response"] == "Responder", "fraction"]
    nr = d.loc[d["response"] == "Non-responder", "fraction"]
    tests.append({"cell_type": ct, "timepoint": tp, "n_R": len(r), "n_NR": len(nr),
                  "median_R": r.median(), "median_NR": nr.median(),
                  "p_mannwhitney": mannwhitneyu(r, nr).pvalue})
tests = pd.DataFrame(tests)
tests.to_csv(RES / "celltype_proportions_tests.csv", index=False)
print(tests.round(3).to_string(index=False))

# Figure: one small panel per cell type; x = timepoint, colour = response
cell_types = props.mean().sort_values(ascending=False).index
ncol = 4
nrow = int(np.ceil(len(cell_types) / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.6 * nrow), squeeze=False)
rng = np.random.default_rng(0)
offsets = {"Responder": -0.18, "Non-responder": 0.18}

for ax, ct in zip(axes.flat, cell_types):
    d = long[long["cell_type"] == ct]
    for xi, tp in enumerate(["Pre", "Post"]):
        for resp, dx in offsets.items():
            y = d.loc[(d["timepoint"] == tp) & (d["response"] == resp), "fraction"]
            x = xi + dx
            ax.boxplot(y, positions=[x], widths=0.28, showfliers=False,
                       medianprops=dict(color=COLORS[resp], lw=2),
                       boxprops=dict(color=MUTED, lw=1), whiskerprops=dict(color=MUTED, lw=1),
                       capprops=dict(color=MUTED, lw=1))
            ax.scatter(x + rng.uniform(-0.07, 0.07, len(y)), y, s=14, color=COLORS[resp],
                       edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_title(ct, fontsize=10, color=INK)
    ax.set_xticks([0, 1], ["Pre", "Post"])
    ax.set_xlim(-0.6, 1.6)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 5, 10]))
    ax.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0, decimals=0))
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(axis="y", color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#c9c8c1")
for ax in list(axes.flat)[len(cell_types):]:
    ax.set_visible(False)

n = info.groupby(["response", "timepoint"], observed=True).size()
handles = [plt.Line2D([], [], marker="o", ls="", color=c,
                      label=f"{r} (samples: Pre {n[r, 'Pre']}, Post {n[r, 'Post']})")
           for r, c in COLORS.items()]
fig.legend(handles=handles, loc="upper center", ncol=2, frameon=False, fontsize=9,
           bbox_to_anchor=(0.5, 1.0))
fig.suptitle("Share of each immune cell type per tumour sample (GSE120575)",
             y=1.04, fontsize=12, color=INK)
fig.supylabel("fraction of CD45+ cells in sample", fontsize=9, color=MUTED)
fig.tight_layout()
fig.savefig(ROOT / "figures" / "celltype_proportions.png", dpi=150, bbox_inches="tight")
print("saved figures/celltype_proportions.png")
