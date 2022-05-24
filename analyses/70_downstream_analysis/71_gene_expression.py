# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [conda env:conda-2021-hairy-cell-leukemia-wolf-scanpy]
#     language: python
#     name: conda-env-conda-2021-hairy-cell-leukemia-wolf-scanpy-py
# ---

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import scanpy as sc
import progeny
import pandas as pd
from statsmodels.stats.multitest import multipletests
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy_helpers as sh
import scipy
import altair as alt

sc.set_figure_params(figsize=(5, 5))

# %% [markdown]
# ## Load data

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")
artifact_dir = "../../data/70_downstream_analyses/gene_expression"

# %%
# !mkdir -p {artifact_dir}

# %%
adata.obs.columns

# %%
# plot fos/jun markers and covariates
sc.pl.embedding(
    adata,
    basis="umap_uncorrected",
    color=[
        "pct_counts_mito",
        "FOS",
        "JUN",
        "RPL37",
        "HSPA1A",
        "ZFP36",
        "n_genes_by_counts",
        "total_counts",
    ],
)

# %% [markdown]
# # Top Genes healhty vs malignant B cell

# %%
de_genes = pd.read_csv(
    "../../data/50_de_analysis/deseq2_results/deseq2_res_bulk_healthy_malignant/healthy_b_cell_malignant_b_cell_IHWallGenes.tsv",
    sep="\t",
)

# %%
top_genes = (
    de_genes.assign(
        direction=lambda x: ["up" if _ > 0 else "down" for _ in x["log2FoldChange"]]
    )
    .groupby("direction")
    .apply(lambda x: x.head(15))
)

# %%
pb_b_cells = sh.pseudobulk.pseudobulk(
    adata[adata.obs["cell_type"].isin(["malignant B cell", "healthy B cell"]), :],
    layer="raw_counts",
    groupby=["patient", "cell_type"],
)

# %%
sc.pp.normalize_total(pb_b_cells, target_sum=1e6)

# %%
sc.pp.log1p(pb_b_cells, base=2)

# %% [markdown]
# ## Matrixplot with zscores

# %%
pb_b_cells.layers["z_scores"] = scipy.stats.zscore(pb_b_cells.X, axis=0)

# %%
fig = sc.pl.matrixplot(
    pb_b_cells,
    var_names=top_genes["gene_id"],
    groupby=["cell_type", "patient"],
    cmap="bwr",
    vmin=-2,
    vmax=2,
    layer="z_scores",
    return_fig=True
)
fig.savefig(f"{artifact_dir}/healthy_malignant_b_cells_pseudobulk_zscore_heatmap.pdf", bbox_inches="tight")

# %% [markdown]
# ### fold change plot

# %%
top_genes

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_b_cells[:, top_genes["gene_id"]],
    groupby="cell_type",
    paired_by="patient",
    metric="diff",
    de_res_df=top_genes,
    metric_name="log2(fold change)"
)

# %%
ch.display()

# %%
ch.save(f"{artifact_dir}/healthy_malignant_b_cells_pseudobulk_fold_changes.svg")

# %%
