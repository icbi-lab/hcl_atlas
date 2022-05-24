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
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import progeny
import pandas as pd
import scipy.stats
import scanpy_helpers as sh

sc.set_figure_params(figsize=(4, 4))

# %%
adata = sc.read_h5ad("../../data/30_merge_adata/adata_scvi_annotated.h5ad")

# %%
adata_malignant_b = sc.read_h5ad(
    "../../data/40_cluster_analysis/adata_malignant_b_cells.h5ad"
)

# %%
artifact_dir = "../../data/70_downstream_analyses/pathway_analysis"

# %%
# !mkdir -p {artifact_dir}

# %% [markdown]
# # Progeny

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
    num_perm=0,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
)

# %% [markdown]
# # Progeny of FOS+ cluster vs rest

# %%
adata_pw = progeny.extract(adata)

# %%
adata_pw.obs["cell_type"].unique()

# %%
pb_pw = sh.pseudobulk.pseudobulk(
    adata_pw[adata_pw.obs["cell_type"].isin(["malignant B cell", "healthy B cell"]), :],
    aggr_fun=np.mean,
    groupby=["patient", "cell_type"],
)
pb_pw.obs["cell_type"] = pb_pw.obs["cell_type"].str.replace(" cell", "")

# %%
pb_pw._sanitize()

# %%
sc.tl.pca(pb_pw)
sc.tl.dendrogram(
    pb_pw, groupby=["patient", "cell_type"], use_rep="X_pca", linkage_method="average"
)

# %%
fig = sc.pl.matrixplot(
    pb_pw,
    var_names=pb_pw.var_names,
    groupby=["patient", "cell_type"],
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
    vmin=-1,
    vmax=1,
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/progeny_healthy_vs_malignant_b_clustered_heatmap.pdf",
    bbox_inches="tight",
)

# %%
pb_pw.obs

# %%
res = (
    sh.compare_groups.lm.test_lm(
        pb_pw,
        "~ C(cell_type, Treatment('healthy B')) + patient",
        groupby="cell_type",
        contrasts="Treatment('healthy B')",
    )
    .pipe(sh.util.fdr_correction)
    .sort_values("fdr")
)

# %%
res

# %%
ch = sh.compare_groups.pl.plot_lm_result_altair(res, p_cutoff=1)
ch.save(f"{artifact_dir}/progeny_healthy_vs_malignant_b_lm_heatmap.svg")
ch.display()

# %%
sh.colors.plot_palette("patient")

# %%
# reorder for proper ordering of legend
pb_pw = pb_pw[
    pb_pw.obs.sort_values(["patient", "cell_type"]).index,
]

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_pw,
    groupby="cell_type",
    paired_by="patient",
    metric="diff",
    var_names=res["variable"].tolist(),
    de_res_df=res,
    pvalue_col="fdr",
    var_col="variable",
)
ch.save(f"{artifact_dir}/progeny_healthy_vs_malignant_b_fold_change_bar_chart.svg")
ch.display()

# %%
res

# %%
sns.set_palette(sns.color_palette(sh.colors.COLORS.patient.values()))
fig = sh.pairwise.plot_paired(
    pb_pw,
    groupby="cell_type",
    paired_by="patient",
    var_names=res["variable"].tolist(),
    ylabel="PROGENy score",
    pvalues=res["fdr"].tolist(),
    pvalue_template="FDR={:.3f}",
    boxplot_kwargs={"color": "white"},
    panel_size=(3, 4.5),
    return_fig=True,
    rotate_x=90,
)
fig.savefig(
    f"{artifact_dir}/progeny_healthy_vs_malignant_b_paired_boxplot.pdf",
    bbox_inches="tight",
)

# %% [markdown]
# ## SR vs LR (timepoint T0)

# %%
pb_pw_sr_lr = sh.pseudobulk.pseudobulk(
    adata_pw[
        (adata_pw.obs["timepoint"] == "T0")
        & (adata_pw.obs["cell_type"] == "malignant B cell"),
        :,
    ],
    groupby=["patient", "response"],
    aggr_fun=np.mean,
)

# %%
pb_pw_sr_lr.obs

# %%
# reorder for proper ordering of legend
pb_pw_sr_lr = pb_pw_sr_lr[
    pb_pw_sr_lr.obs.sort_values(["patient", "response"]).index,
]

# %%
res = sh.compare_groups.lm.test_lm(
    pb_pw_sr_lr,
    "~ C(response, Treatment('long_term'))",
    groupby="response",
    contrasts="Treatment('long_term')",
).sort_values("pvalue")

# %%
res

# %%
sns.set_palette(
    sns.color_palette(
        [sh.colors.COLORS.patient[p] for p in pb_pw_sr_lr.obs["patient"].unique()]
    )
)
fig = sh.pairwise.plot_paired(
    pb_pw_sr_lr,
    groupby="response",
    var_names=res["variable"].tolist(),
    ylabel="PROGENy score",
    pvalues=res["pvalue"].tolist(),
    boxplot_kwargs={"color": "white"},
    panel_size=(3, 4),
    rotate_x=90,
    hue="patient",
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/progeny_long_term_vs_short_term_t0_paired_boxplot.pdf",
    bbox_inches="tight",
)

# %% [markdown]
# ## malignant B vs Fos malignant B

# %%
progeny.run(
    adata_malignant_b,
    model,
    center=True,  # Center gene expression by mean per cell
    num_perm=0,  # Simulate m random activities
    norm=True,  # Normalize by number of edges to correct for large regulons
    scale=True,  # Scale values per feature so that values can be compared across cells
)

# %%
adata_pw_malignant_b = progeny.extract(adata_malignant_b)

# %%
pb_pw_malignant_b = sh.pseudobulk.pseudobulk(
    adata_pw_malignant_b,
    aggr_fun=np.mean,
    groupby=["patient", "cell_phenotype"],
)

# %%
pb_pw_malignant_b.obs["cell_phenotype"] = pd.Categorical(
    pb_pw_malignant_b.obs["cell_phenotype"],
    categories=["malignant_b", "fos_malignant_b"],
)

# %%
pb_pw_malignant_b._sanitize()
sc.tl.pca(pb_pw_malignant_b)
sc.tl.dendrogram(
    pb_pw_malignant_b,
    groupby=["patient", "cell_phenotype"],
    use_rep="X_pca",
    linkage_method="average",
)

# %%
fig = sc.pl.matrixplot(
    pb_pw_malignant_b,
    var_names=pb_pw_malignant_b.var_names,
    groupby=["patient", "cell_phenotype"],
    cmap="bwr",
    swap_axes=True,
    dendrogram=True,
    vmin=-2,
    vmax=2,
    return_fig=True,
)
fig.savefig(
    f"{artifact_dir}/progeny_malignant_b_vs_fos_malignant_b_clustered_heatmap.pdf",
    bbox_inches="tight",
)

# %%
res = (
    sh.compare_groups.lm.test_lm(
        pb_pw_malignant_b,
        "~ C(cell_phenotype, Treatment('malignant_b')) + patient",
        groupby="cell_phenotype",
        contrasts="Treatment('malignant_b')",
    )
    .pipe(sh.util.fdr_correction)
    .sort_values("fdr")
)

# %%
res

# %%
ch = sh.compare_groups.pl.plot_lm_result_altair(res, p_cutoff=1)
ch.save(f"{artifact_dir}/progeny_malignant_b_vs_fos_malignant_b_lm_heatmap.svg")
ch.display()

# %%
# reorder for proper ordering of legend
pb_pw_malignant_b = pb_pw_malignant_b[
    pb_pw_malignant_b.obs.sort_values(["patient", "cell_phenotype"]).index,
]
pb_pw

# %%
ch = sh.pairwise.plot_paired_fc(
    pb_pw_malignant_b,
    groupby="cell_phenotype",
    paired_by="patient",
    metric="diff",
    var_names=res["variable"].tolist(),
    de_res_df=res,
    pvalue_col="fdr",
    var_col="variable",
)
ch.save(f"{artifact_dir}/progeny_malignant_b_vs_fos_malignant_b_fold_change_barchart.svg")
ch.display()

# %%
sns.set_palette(sns.color_palette(sh.colors.COLORS.patient.values()))
fig = sh.pairwise.plot_paired(
    pb_pw_malignant_b,
    groupby="cell_phenotype",
    paired_by="patient",
    var_names=res["variable"].tolist(),
    ylabel="PROGENy score",
    pvalues=res["fdr"].tolist(),
    pvalue_template="FDR={:.3f}",
    boxplot_kwargs={"color": "white"},
    panel_size=(3, 4.5),
    rotate_x=90,
    return_fig=True
)
fig.savefig(f"{artifact_dir}/progeny_malignant_b_vs_fos_malignant_b_paired_boxplot.pdf", bbox_inches="tight")

# %%
