# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: SSH apollo-15 scVI GPU (ssh)
#     language: ''
#     name: rik_ssh_apollo_15_scvigpussh
# ---

# %%
import sys
print(sys.path)

# %%
import scanpy as sc
from glob import glob
import anndata
import scvi
import progeny

sc.set_figure_params(figsize=(4, 4))
import pandas as pd

# %%
patient_data = pd.read_csv("../../tables/patient_data.csv")

# %%
adatas = [
    sc.read_h5ad(f)
    for f in glob("../../data/20_scrnaseq_qc/01_qc_and_filtering/**/*.h5ad")
]

# %%
adata = anndata.concat(adatas, index_unique="_")
adata.layers["raw_counts"] = adata.X.copy()

# %%
patient_data

# %%
adata.obs = (
    adata.obs.reset_index()
    .merge(patient_data, how="left", left_on="patient", right_on="patient_id")
    .set_index("cell_id")
)

# %%
adata_scvi = adata.copy()

# %%
adata

# %%
adata.obs

# %%
sc.pp.highly_variable_genes(
    adata, flavor="seurat_v3", n_top_genes=6000, batch_key="patient"
)

# %%
sc.pp.normalize_total(adata)

# %%
sc.pp.log1p(adata)

# %%
sc.tl.pca(adata)

# %%
sc.pp.neighbors(adata)

# %%
sc.tl.umap(adata)

# %%
sc.tl.leiden(adata)

# %%
sc.tl.paga(adata)

# %%
sc.pl.paga(adata)

# %%
sc.tl.umap(adata, init_pos="paga")

# %%
sc.pl.umap(
    adata, color=["sample", "patient", "timepoint", "leiden"], ncols=2, wspace=0.6
)

# %%
adata.obs["sample"].value_counts()

# %% [markdown]
# # scVI

# %%
adata_scvi = adata_scvi[adata_scvi.obs["sample"] != "HCL_P1_T0", :].copy()
adata_scvi._sanitize()

# %%
sc.pp.highly_variable_genes(
    adata_scvi, flavor="seurat_v3", n_top_genes=6000, batch_key="patient", subset=True
)

# %%
scvi.data.setup_anndata(adata_scvi, batch_key="sample")

# %%
scvi.data.view_anndata_setup(adata_scvi)

# %%
# vae = scvi.model.SCVI(adata_scvi, use_cuda=True)

# %%
# vae.train()

# %%
# vae.save("scvi_model")

# %%
vae = scvi.model.SCVI.load("scvi_model", adata=adata_scvi)

# %%
adata_scvi.obsm["X_scVI"] = vae.get_latent_representation()

# %%
sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
sc.tl.leiden(adata_scvi)
sc.tl.umap(adata_scvi)

# %%
sc.pl.umap(
    adata_scvi, color=["sample", "patient", "timepoint", "leiden"], ncols=2, wspace=0.6
)

# %%
adata = adata[adata_scvi.obs_names, :].copy()

# %%
adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"]
adata.obsm["X_umap"] = adata_scvi.obsm["X_umap"]
adata.obsm["X_scVI"] = adata_scvi.obsm["X_scVI"]

# %%
sc.tl.leiden(adata, resolution=0.5)
adata.obs["leiden_scvi"] = adata.obs["leiden"]

# %%
# !mkdir -p "../../data/30_merge_adata/"

# %% tags=[]
adata.write_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")
