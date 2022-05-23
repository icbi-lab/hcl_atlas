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
import scanpy_helpers as sh
import matplotlib.pyplot as plt

# %%
adata_malignant_b = sc.read_h5ad("../../data/40_cluster_analysis/adata_malignant_b_cells.h5ad")

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")

# %%
artifact_dir = "../../data/70_downstream_analyses/overview_plots"

# %%
# !mkdir -p {artifact_dir}

# %%
sh.colors.set_scale_anndata(adata, "patient")
sh.colors.set_scale_anndata(adata, "timepoint")
sh.colors.set_scale_anndata(adata, "cell_type")
sh.colors.set_scale_anndata(adata, "response")

# %%
with plt.rc_context({"figure.figsize": (7, 7), "figure.dpi": 300}):
    for basis in ["umap", "umap_uncorrected"]:
        for color in ["patient", "timepoint", "cell_type", "response"]:
            fig = sc.pl.embedding(adata, color=color, frameon=False, basis=basis, size=7, return_fig=True)
            fig.savefig(f"{artifact_dir}/overview_{basis}_{color}.pdf", bbox_inches="tight")
            fig.show()

# %%
