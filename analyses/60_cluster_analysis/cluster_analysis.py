# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
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

sc.set_figure_params(figsize=(4, 4))

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")

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
    sc.pl.umap(tmp_adata, color=["patient", "timepoint", "leiden"], ncols=4)
    patient_adatas[patient] = tmp_adata

# %%
for patient, tmp_adata in patient_adatas.items():
    tmp_adata.write_h5ad(f"../../data/60_cluster_analysis/60_adata_{patient}.h5ad")

# %%
for tmp_adata in patient_adatas.values():
    sc.tl.rank_genes_groups(tmp_adata, method="wilcoxon", groupby="leiden")

# %% [markdown]
# It appears that in all patients, there is a cluster that is FOSB+/JUN+/DUSP1+.

# %%
for tmp_patient, tmp_adata in patient_adatas.items():
    print(tmp_patient)
    sc.pl.rank_genes_groups_dotplot(tmp_adata, n_genes=5, title=tmp_patient)

# %%
for tmp_patient in ["P1", "P2", "P3"]:
    tmp_adata = patient_adatas[tmp_patient]
    sc.tl.rank_genes_groups(tmp_adata, method="wilcoxon", groupby="timepoint")
    sc.pl.rank_genes_groups_dotplot(tmp_adata, n_genes=10, title=tmp_patient)

# %%
sc.pl.embedding(adata, "umap_uncorrected", color=["DUSP1", "FOSB", "JUN"])

# %%
for tmp_patient, tmp_adata in patient_adatas.items():
    print(tmp_patient)
    sc.pl.umap(tmp_adata, color=["leiden", "FOSB", "DUSP1", "KLF4"])

# %%
# for each patient, the clusters that are FOS+ JUN+
fos_annotation = {"P1": [5], "P2": [9, 2], "P3": [8], "P4": [0], "P5": [9], "P6": [3]}

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
adata.write_h5ad("../../data/60_cluster_analysis/adata_malignant_b_cells.h5ad")

# %% [markdown]
# ## Plot cluster abundance by timepoint

# %%
for patient in ["P1", "P2", "P3"]:
    cells_by_timepoint = (
        patient_adatas[patient]
        .obs.groupby(["timepoint", "leiden"])
        .size()
        .reset_index(name="n_cells")
        .groupby(["timepoint"])
        .apply(lambda df: df.assign(n_cells_per_timepoint=np.sum(df["n_cells"])))
        .pipe(
            lambda df: df.assign(frac_cells=df["n_cells"] / df["n_cells_per_timepoint"])
        )
    )
    ax = sns.lineplot(
        data=cells_by_timepoint, x="timepoint", y="frac_cells", hue="leiden"
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0, title="cluster")
    ax.set_ylabel("fraction of cells")
    ax.set_title(patient)
    plt.show()

# %%
import matplotlib.colors

# %%
for patient, tmp_adata in patient_adatas.items():
    sc.pl.dotplot(
        tmp_adata,
        var_names=["STAG3"],
        groupby="leiden",
        dot_max=1,
        swap_axes=True,
        title=patient,
    )

# %%
adata.obs

# %% [markdown]
# ### Progeny analysis of cell phenotypes

# %%
model = progeny.load_model(
    organism="Human",  # If working with mouse, set to Mouse
    top=1000,  # For sc we recommend ~1k target genes since there are dropouts
)

# %%
progeny.run(
    adata,
    model,
    center=True,  # Center gene expression by mean per cell
    num_perm=100,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
)

# %%
adata_pw = progeny.extract(adata)

# %%
sc.pl.matrixplot(
    adata_pw,
    var_names=adata.obsm["progeny"].columns,
    groupby="cell_phenotype",
    swap_axes=True,
    cmap="bwr",
    vmin=-1,
    vmax=1,
)

# %%
adata_pw.obs["cell_phenotype_by_patient"] = [
    f"{cpt}_{patient}"
    for patient, cpt in zip(adata_pw.obs["patient"], adata_pw.obs["cell_phenotype"])
]

# %%
sc.pl.matrixplot(
    adata_pw,
    var_names=adata.obsm["progeny"].columns,
    groupby="cell_phenotype_by_patient",
    swap_axes=True,
    cmap="bwr",
    vmin=-1,
    vmax=1,
    dendrogram=False,
)

# %%
progeny_df = adata.obsm["progeny"].join(adata.obs.loc[:, ["cell_phenotype", "patient"]])
progeny_by_sample = (
    progeny_df.groupby(["cell_phenotype", "patient"], observed=True)
    .agg("mean")
    .reset_index()
    .sort_values("patient")
)

# %%
progeny_by_sample

# %%
progeny_by_sample_melt = progeny_by_sample.melt(
    id_vars=["cell_phenotype", "patient"],
    var_name="pathway",
    value_name="progeny_score",
)
progeny_by_sample_melt["response"] = pd.Categorical(
    progeny_by_sample_melt["cell_phenotype"],
    categories=["fos_malignant_b", "malignant_b"],
)

# %%
progeny_by_sample_melt

# %%
df_for_test = progeny_by_sample.drop("patient", axis="columns").set_index(
    "cell_phenotype"
)

# %%
_, pvals = scipy.stats.ttest_rel(
    df_for_test.loc["fos_malignant_b", :], df_for_test.loc["malignant_b", :], axis=0
)
pval_dict = {k: v for k, v in zip(df_for_test.columns, pvals)}

# %%
fig, axes = plt.subplots(4, 4, figsize=(4 * 3, 4 * 4), tight_layout=True)
axes = axes.flatten()
for pw, ax in zip(adata_pw.var.index, axes):
    sns.stripplot(
        x="cell_phenotype",
        data=progeny_by_sample_melt.query(f"pathway == '{pw}'"),
        y="progeny_score",
        ax=ax,
        hue="patient",
        size=10,
        linewidth=2,
    )
    sns.lineplot(
        x="cell_phenotype",
        data=progeny_by_sample_melt.query(f"pathway == '{pw}'"),
        hue="patient",
        y="progeny_score",
        ax=ax,
    )
    sns.boxplot(
        x="cell_phenotype",
        data=progeny_by_sample_melt.query(f"pathway == '{pw}'"),
        y="progeny_score",
        ax=ax,
    )

    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=10)
    ax.legend().set_visible(False)
    ax.set_ylabel(pw)
    ax.set_title("unadj. p={:.2f}, t-test".format(pval_dict[pw]))
axes[-1].set_visible(False)
axes[-2].set_visible(False)
axes[-3].legend().set_visible(True)
axes[-3].legend(bbox_to_anchor=(1.1, 1.05))
fig.tight_layout()
plt.show()

# %%
