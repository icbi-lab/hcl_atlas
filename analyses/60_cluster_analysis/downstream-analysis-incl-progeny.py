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
import progeny
import pandas as pd
from statsmodels.stats.multitest import multipletests
import scipy.stats
import matplotlib.pyplot as plt
import seaborn as sns

sc.set_figure_params(figsize=(5, 5))

# %% [markdown]
# ## Load data

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi.h5ad")

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

# %%
adata = adata[adata.obs["timepoint"] == "T0", :]

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
def rank_pw_groups_all(adata, groupby):
    res = (
        pd.concat(
            [
                progeny.rank_pws_groups(adata, groupby=groupby, group=ct)
                for ct in adata.obs[groupby].unique()
            ]
        )
        .reset_index()
        .sort_values("name")
    )
    _, padj, _, _ = multipletests(res["pval"].values, alpha=0.05, method="fdr_bh")
    res["pval_adj"] = padj
    return res


# %% [markdown]
# ## IG loci
#
# Only look at H cells! 

# %%
sc.pl.embedding(adata, "umap_uncorrected", color="cell_type")

# %%
adata_h = adata[adata.obs["cell_type"] == "malignant B cell", :]

# %%
sc.pl.embedding(adata_h, basis="umap_uncorrected", color=["leiden", "patient"])

# %%
sc.pl.DotPlot(
    adata_h,
    groupby="patient",
    var_names=adata.var_names[adata.var_names.str.startswith("IGH")],
).show()

# %%
sc.pl.DotPlot(
    adata_h,
    groupby="patient",
    var_names=adata.var_names[adata.var_names.str.startswith("IGLV")],
).show()

# %% [markdown]
# ### P6 diff between clones

# %%
adata_p6 = adata_h[adata_h.obs["patient"] == "P6"]
adata_p6.obs["clone"] = [
    {"6": "clone 1", "7": "clone 2"}.get(x, None) for x in adata_p6.obs["leiden"]
]

# %%
sc.pl.DotPlot(
    adata_p6,
    groupby="clone",
    var_names=adata.var_names[adata.var_names.str.startswith("IGH")],
).show()

# %%
sc.pl.DotPlot(
    adata_p6,
    groupby="clone",
    var_names=adata.var_names[adata.var_names.str.startswith("IGLV")],
).show()

# %%
adata.obs["p6_clone"] = adata_p6.obs["clone"]

# %%
sc.pl.embedding(adata, basis="umap_uncorrected", color="p6_clone")

# %% [markdown]
# ### CD27 expression
# CD27 sould only be expressed on H cells

# %%
sc.pl.dotplot(adata, var_names=["CD27"], groupby="cell_type")

# %% [markdown]
# CD27 expression on healthy B cells, by patient

# %%
sc.pl.dotplot(
    adata[adata.obs["cell_type"] == "healthy B cell", :],
    var_names=["CD27"],
    groupby="patient",
)

# %% [markdown]
# ### Marker dotplot from B-cell atlas

# %%
b_cell_atlas_genes = [
    "HOPX",
    "PDE4D",
    "IGHE",
    "SELL",
    "EMP3",
    "CIB1",
    "PSAP",
    "CD72",
    "DAPP1",
    "LTB",
    "HCK",
    "ZEB2",
    "RHOB",
    "TNFRSF1B",
    "FCRL3",
    "FCRL5",
    "FGR",
    "MPP6",
    "TAGLN2",
    "IGHA2",
    "AHNAK",
    "S100A4",
    "CRIP2",
    "ITGB1",
    "JCHAIN",
    "VIM",
    "PLPP5",
    "FCER2",
    "IL4R",
    "CRIP1",
    "LGALS1",
    "CTSH",
    "S100A10",
    "IGHG2",
    "VPREB3",
    "PPP1R14A",
    "PCDH9",
    "PLD4",
    "IGHM",
    "MT-ATP8",
    "IGHD",
    "SOX4",
    "AL139020.1",
    "IGLL5",
    "TCL1A",
]

# %%
sc.pl.dotplot(
    adata,
    var_names=b_cell_atlas_genes,
    swap_axes=False,
    groupby="cell_type",
)

# %% [markdown]
# ## Gene heatmap

# %%
tmp_adata  = adata.copy()
tmp_adata.obs["patient_cell_type"] = [
    f"malignant B {p}" if ct == "malignant B cell" else ct for p, ct in zip(adata.obs["patient"], adata.obs["cell_type"])
]

# %%
sc.pl.heatmap(
    tmp_adata,
    var_names="""CD27
ITGAX
BCL2
CCND1
CDKN1B
FLT3
IL3RA
PRKCA
IL1R2
GAS7
JUND
FGFR1
CXCR4
CXCR5
CCR7
ITGB1
TIMP1
TIMP4
TNFRSF1A
TNFRSF1B
LSP1
CDC42
RAC1
RHOA
SYT1
FGF2
""".split(),
    groupby="patient_cell_type",
    swap_axes=True, figsize=(10, 8)
)

# %% [markdown]
# ### Heatmap DE malignant vs. healthy

# %%
de_res = pd.read_csv("../../data/70_de_analysis/72_run_de/deseq2_res_sc_healthy_vs_malignant_b_cells/malignant_healthy_IHWallGenes.tsv", sep="\t")

# %%
sc.pl.heatmap(
    tmp_adata,
    var_names=de_res.loc[de_res["log2FoldChange"] > 0, :].sort_values("padj")["gene_id"][:30].values,
    groupby="patient_cell_type",
    swap_axes=True, figsize=(10, 8)
)

# %% [markdown]
# ## Compare pathway scores by cell-type
#
#  * statistical unit = "cell"
#  * all cell-types
#  * "T0" only

# %%
sc.pl.matrixplot(
    adata_pw,
    var_names=adata.obsm["progeny"].columns,
    groupby="cell_type",
    swap_axes=True,
    cmap="bwr",
    title="cell_type",
    vmin=-1,
    vmax=1,
)

# %% [markdown]
# ### healthy and malignant B cells by patient
#  * excluding P3 since it only has 3 healthy B cells

# %%
adata_pw.obs.loc[
    adata_pw.obs["cell_type"].isin(["healthy B cell", "malignant B cell"]), :
].groupby(["patient", "cell_type"], observed=True).size().reset_index()

# %%
adata_pw_ct_patient = adata_pw[
    adata_pw.obs["cell_type"].isin(["healthy B cell", "malignant B cell"])
    & (adata.obs["patient"] != "P3"),
    :,
]
adata_pw_ct_patient.obs["ct_patient"] = [
    f"{p} {ct}"
    for p, ct in zip(
        adata_pw_ct_patient.obs["patient"], adata_pw_ct_patient.obs["cell_type"]
    )
]

# %%
sc.pl.matrixplot(
    adata_pw_ct_patient,
    var_names=adata.obsm["progeny"].columns,
    groupby="ct_patient",
    swap_axes=True,
    cmap="bwr",
    title="ct_patient",
    vmin=-1,
    vmax=1,
)

# %%
sc.pl.matrixplot(
    adata_pw_ct_patient,
    var_names=adata.obsm["progeny"].columns,
    groupby="ct_patient",
    swap_axes=True,
    cmap="bwr",
    title="ct_patient",
    vmin=-1,
    vmax=1,
    dendrogram=True,
)

# %%
sc.pl.umap(
    adata_pw,
    color=["cell_type", "PI3K", "WNT", "Androgen", "Estrogen", "MAPK", "EGFR"],
    cmap="coolwarm",
    size=15,
    ncols=3,
    wspace=0.5,
    vmin=-5,
    vmax=5,
)

# %% [markdown]
# ### statistical tests for pathway differences between cell-types

# %%
pw_test_cell_type = rank_pw_groups_all(adata_pw, "cell_type").sort_values("pval")
pw_test_cell_type

# %% [markdown]
# ## Pathway differences between short-term and long-term responders
#  * statistical unit = "patient"
#  * use "malignant B cells" only 
#  * use "T0" only

# %%
adata = adata[adata.obs["cell_type"] == "malignant B cell", :]

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
progeny_df = adata.obsm["progeny"].join(adata.obs.loc[:, ["patient", "response"]])

# %%
progeny_by_patient = (
    progeny_df.groupby(["patient", "response"], observed=True)
    .agg("mean")
    .reset_index()
    .sort_values("patient")
)

# %%
progeny_by_patient_melt = progeny_by_patient.melt(
    id_vars=["patient", "response"], var_name="pathway", value_name="progeny_score"
)
progeny_by_patient_melt["response"] = pd.Categorical(
    progeny_by_patient_melt["response"], categories=["short_term", "long_term"]
)

# %%
df_for_test = progeny_by_patient.drop("patient", axis="columns").set_index("response")

# %%
_, pvals = scipy.stats.ttest_ind(
    df_for_test.loc["short_term", :], df_for_test.loc["long_term", :], axis=0
)
pval_dict = {k: v for k, v in zip(df_for_test.columns, pvals)}

# %%
fig, axes = plt.subplots(5, 3, figsize=(3 * 3, 5 * 3), tight_layout=True)
axes = axes.flatten()
for pw, ax in zip(adata_pw.var.index, axes):
    sns.stripplot(
        x="response",
        data=progeny_by_patient_melt.query(f"pathway == '{pw}'"),
        y="progeny_score",
        ax=ax,
        hue="patient",
        size=10,
        linewidth=2,
    )
    sns.boxplot(
        x="response",
        data=progeny_by_patient_melt.query(f"pathway == '{pw}'"),
        y="progeny_score",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=90, labelsize=10)
    ax.legend().set_visible(False)
    ax.set_ylabel(pw)
    ax.set_title("unadj. p={:.2f}, t-test".format(pval_dict[pw]))
axes[-1].set_visible(False)
axes[-2].legend().set_visible(True)
axes[-2].legend(bbox_to_anchor=(1.1, 1.05))
fig.tight_layout()
plt.show()

# %% [markdown]
# ### matrixplot by response

# %%
sc.pl.matrixplot(
    adata_pw,
    var_names=adata.obsm["progeny"].columns,
    groupby="response",
    swap_axes=True,
    cmap="bwr",
    title="cell_type",
    vmin=-0.5,
    vmax=0.5,
)

# %% [markdown]
# ### matrixplot by patient

# %%
sc.pl.matrixplot(
    adata_pw,
    var_names=adata.obsm["progeny"].columns,
    groupby="patient",
    swap_axes=True,
    cmap="bwr",
    title="patient",
    vmin=-0.5,
    vmax=0.5,
)

# %%
