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

sc.settings.set_figure_params(figsize=(4, 4))
import pandas as pd
import numpy as np
import itertools
import scanpy_helpers as sh

# %%
adata_all = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")
adata_malignant_b = sc.read_h5ad(
    "../../data/40_cluster_analysis/adata_malignant_b_cells.h5ad"
)
artifact_dir = "../../data/50_de_analysis/pseudobulk"

# %%
adata = adata_all[adata_all.obs["cell_type"] == "malignant B cell", :].copy()

# %%
adata.obs

# %%
adata.obs.groupby(
    ["patient", "timepoint", "sample", "response"], observed=True
).size().reset_index(name="size").sort_values(["timepoint", "patient"])

# %% [markdown]
# ## Comparison between short- and long-term responders
# * Only T0
# * All timepoints

# %%
sh.pseudobulk.write_pseudobulk(
    sh.pseudobulk.pseudobulk(
        adata, groupby=["patient", "response", "sex", "age"], layer="raw_counts"
    ),
    f"{artifact_dir}/bulk_response_all_timepoints",
)

# %%
sh.pseudobulk.write_pseudobulk(
    sh.pseudobulk.pseudobulk(
        adata[adata.obs["timepoint"] == "T0", :],
        groupby=["patient", "response", "sex", "age"],
        layer="raw_counts",
    ),
    f"{artifact_dir}/bulk_response_t0",
)

# %% [markdown]
# ## DE analysis of timepoints (pseudobulk, only P2-3)

# %%
tmp_adata = adata[adata.obs["patient"].isin(["P2", "P3"]), :]
tmp_adata.obs["timepoint"] = [
    "pre-treatment" if t == "T0" else "post-treatment"
    for t in tmp_adata.obs["timepoint"]
]

# %%
sh.pseudobulk.write_pseudobulk(
    sh.pseudobulk.pseudobulk(
        tmp_adata, groupby=["patient", "timepoint", "sex", "age"], layer="raw_counts"
    ),
    f"{artifact_dir}/bulk_timepoints",
)

# %% [markdown]
# ## Comparison of JUN+ cluster vs rest of malignant B cells

# %%
adata_malignant_b.obs

# %%
sh.pseudobulk.write_pseudobulk(
    sh.pseudobulk.pseudobulk(
        adata_malignant_b, groupby=["patient", "cell_phenotype"], layer="raw_counts"
    ),
    f"{artifact_dir}/bulk_fos_jun_vs_rest",
)

# %% [markdown]
# ## Comparison of healthy and malignant B cells

# %%
adata_healthy_malignant = adata_all[
    adata_all.obs["cell_type"].isin(["malignant B cell", "healthy B cell"]), :
].copy()
adata_healthy_malignant.obs["cell_type"] = (
    adata_healthy_malignant.obs["cell_type"].str.lower().str.replace(" ", "_")
)
sh.pseudobulk.write_pseudobulk(
    sh.pseudobulk.pseudobulk(
        adata_healthy_malignant,
        groupby=["cell_type", "patient"],
        layer="raw_counts",
    ),
    f"{artifact_dir}/bulk_healthy_malignant",
)

# %%
