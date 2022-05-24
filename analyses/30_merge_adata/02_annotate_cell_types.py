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
import scanpy as sc
sc.set_figure_params(figsize=(4,4))
import scanpy_helpers as sh

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")

# %% [markdown] tags=[]
# ### Cell-type annotation

# %%
sc.pl.umap(adata, color=["EPCAM", "CD3E", "CDK1", "CD79A", "SDC1", "MZB1", "FLT3", "FAM129C"], ncols=2)

# %%
sc.pl.umap(adata, color="leiden_scvi", legend_loc="on data")

# %%
sh.annotation.annotate_cell_types(adata, {
    "T cell": [9],
    "Plasma cell": [10],
    "healthy B cell": [6],
    "malignant B cell (dividing)": [8],
    "malignant B cell": [0, 1, 3, 5, 4, 2, 7]
}, column="leiden_scvi")

# %%
adata.write_h5ad("../../data/30_merge_adata/adata_scvi_annotated.h5ad")
