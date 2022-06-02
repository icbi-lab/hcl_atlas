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
adata = sc.read_h5ad("../../data/40_cluster_analysis/adata_malignant_b_cells.h5ad")
artifact_dir = "../../data/70_downstream_analyses/cluster_analysis"

# %%
# !mkdir -p {artifact_dir}

# %%
patient_adatas = {
    p: sc.read_h5ad(f"../../data/40_cluster_analysis/adata_{p}.h5ad")
    for p in adata.obs["patient"].unique()
}

# %% [markdown]
# # Cluster abundance by timepoint

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
    fig = ax.get_figure()
    fig.savefig(f"{artifact_dir}/cluster_abundance_{patient}.pdf", bbox_inches="tight")
    plt.show()

# %% [markdown]
# ## single panel

# %%
fractions_timepoint = sh.util.cell_type_fractions(
    adata, ["patient", "timepoint", "cell_phenotype"], ["patient", "timepoint"]
)

# %%
fractions_timepoint.to_csv(f"{artifact_dir}/fos_jun_fractions_per_timepoint.csv")

# %%
fractions_timepoint

# %%
sns.set_palette(sns.color_palette(sh.colors.COLORS.patient.values()))
tmp_df = fractions_timepoint.loc[lambda x: x["cell_phenotype"] == "fos_malignant_b"]
ax = sns.lineplot(
    x="timepoint",
    y="frac_cells",
    data=tmp_df,
    hue="patient",
    legend=False,
)
sns.stripplot(
    x="timepoint",
    y="frac_cells",
    data=tmp_df,
    hue="patient",
    linewidth=1,
    ax=ax,
    size=10,
    jitter=False,
    alpha=0.8,
)
ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
ax.set_ylabel("fraction FOS/JUN+ cells")
ax.get_figure().savefig(f"{artifact_dir}/fos_jun_abundance.pdf", bbox_inches="tight")

# %% [markdown]
# ## fractions per patient

# %%
fractions_per_patient = sh.util.cell_type_fractions(adata, ["patient", "cell_phenotype"], ["patient"])

# %%
fractions_per_patient

# %%
fractions_timepoint.to_csv(f"{artifact_dir}/fos_jun_fractions_per_patient.csv")
