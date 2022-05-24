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
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import progeny
import pandas as pd
import scipy.stats
import scanpy_helpers as sh

sc.set_figure_params(figsize=(4, 4))

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")
artifact_dir = "../../data/40_cluster_analysis"

# %%
sh.colors.set_scale_anndata(adata, "patient")
sh.colors.set_scale_anndata(adata, "timepoint")
sh.colors.set_scale_anndata(adata, "cell_type")
sh.colors.set_scale_anndata(adata, "response")

# %%
sc.pl.umap(adata, color="cell_type")

# %%
adata = adata[adata.obs["cell_type"] == "malignant B cell", :]

# %%
sc.pl.embedding(adata, "umap_uncorrected", color=["patient", "timepoint"])

# %%
patient_adatas = dict()
# for patient in ["P1", "P2", "P3"]:
for patient in adata.obs["patient"].sort_values().unique():
    print(patient)
    tmp_adata = adata[adata.obs["patient"] == patient, :].copy()
    #     progeny.run(tmp_adata, scale=True)
    #     for col in tmp_adata.obsm["progeny"].columns:
    #         tmp_adata.obs[f"PW:{col}"] = tmp_adata.obsm["progeny"][col]
    sc.tl.pca(tmp_adata)
    sc.pp.neighbors(tmp_adata)
    sc.tl.umap(tmp_adata)
    sc.tl.leiden(tmp_adata, resolution=1)
    fig = sc.pl.umap(tmp_adata, color=["patient", "timepoint", "leiden"], ncols=4, return_fig=True, frameon=False)
    fig.savefig(f"{artifact_dir}/umap_{patient}.pdf", bbox_inches="tight", dpi=1200)
    patient_adatas[patient] = tmp_adata

# %%
for patient, tmp_adata in patient_adatas.items():
    tmp_adata.write_h5ad(f"{artifact_dir}/adata_{patient}.h5ad")

# %%
for tmp_adata in patient_adatas.values():
    sc.tl.rank_genes_groups(tmp_adata, method="wilcoxon", groupby="leiden")

# %% [markdown]
# It appears that in all patients, there is a cluster that is FOSB+/JUN+/DUSP1+.

# %%
for tmp_patient, tmp_adata in patient_adatas.items():
    print(tmp_patient)
    fig = sc.pl.rank_genes_groups_dotplot(tmp_adata, n_genes=5, title=tmp_patient, return_fig=True)
    fig.savefig(f"{artifact_dir}/dotplot_{tmp_patient}.pdf", bbox_inches="tight")

# %%
fig = sc.pl.embedding(adata, "umap_uncorrected", color=["DUSP1", "FOSB", "JUN"], frameon=False, return_fig=True)
fig.savefig(f"{artifact_dir}/umap_fos_dusp_jun_expression.pdf", bbox_inches="tight", dpi=1200)

# %%
for tmp_patient, tmp_adata in patient_adatas.items():
    print(tmp_patient)
    sc.pl.umap(tmp_adata, color=["leiden", "FOSB", "DUSP1", "KLF4"])

# %%
# for each patient, the clusters that are FOS+ JUN+
fos_annotation = {"P1": [6], "P2": [6, 2], "P3": [4], "P4": [0], "P5": [6], "P6": [4]}

# %%
adata.obs["cell_phenotype"] = "malignant_b"
for tmp_patient, tmp_adata in patient_adatas.items():
    adata.obs.loc[
        tmp_adata[
            tmp_adata.obs["leiden"].isin([str(x) for x in fos_annotation[tmp_patient]])
        ].obs_names,
        "cell_phenotype",
    ] = "fos_malignant_b"

# %%
sc.pl.embedding(adata, "umap_uncorrected", color=["cell_phenotype", "FOSB"])

# %%
adata.write_h5ad(f"{artifact_dir}/adata_malignant_b_cells.h5ad")
