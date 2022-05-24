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
import itertools

sc.settings.set_figure_params(vector_friendly=True, dpi=1200)

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi_annotated.h5ad")

# %%
artifact_dir = "../../data/70_downstream_analyses/overview_plots"

# %%
# !mkdir -p {artifact_dir}

# %%
sh.colors.set_scale_anndata(adata, "patient")
sh.colors.set_scale_anndata(adata, "timepoint")
sh.colors.set_scale_anndata(adata, "cell_type")
sh.colors.set_scale_anndata(adata, "response")

# %% tags=[]
with plt.rc_context({"figure.figsize": (7, 7), "figure.dpi": 300}):
    for basis in ["umap", "umap_uncorrected"]:
        for color in ["patient", "timepoint", "cell_type", "response"]:
            fig = sc.pl.embedding(
                adata,
                color=color,
                frameon=False,
                basis=basis,
                size=7,
                return_fig=True,
                show=False,
            )
            fig.savefig(
                f"{artifact_dir}/overview_{basis}_{color}.pdf",
                bbox_inches="tight",
                dpi=1200,
            )
            fig.show()

# %% [markdown]
# ## Marker genes

# %%
marker_genes = {
    "T cell": ["CD3E"],
    "dividing cells": ["CDK1"],
    "B lineage": ["CD79A"],
    "B cell healthy": ["FAM129C"],
    "B cell malignant": ["FLT3"],
    "Plasma cell": ["MZB1", "SDC1"],
}

# %%
fig = sc.pl.dotplot(adata, var_names=marker_genes, groupby="cell_type", return_fig=True)
fig.savefig(f"{artifact_dir}/cell_type_markers_dotplot.pdf", bbox_inches="tight")

# %%
fig = sc.pl.umap(
    adata,
    color=list(itertools.chain.from_iterable(marker_genes.values())),
    return_fig=True,
)
fig.savefig(f"{artifact_dir}/cell_type_markers_umap.pdf", dpi=600, bbox_inches="tight")

# %% [markdown]
# ## cell stats

# %%
sh.util.cell_type_fractions(
    adata, ["patient", "timepoint", "cell_type"], ["patient", "timepoint"]
).to_csv(f"{artifact_dir}/cell_type_fractions_patient_timepoint.csv")

# %%
sh.util.cell_type_fractions(
    adata, ["patient", "cell_type"], ["patient"]
).to_csv(f"{artifact_dir}/cell_type_fractions_patient.csv")

# %%
sh.util.cell_type_fractions(
    adata, ["patient", "timepoint", "response"], ["patient", "response"]
).sort_values(["patient", "response", "timepoint

# %%
