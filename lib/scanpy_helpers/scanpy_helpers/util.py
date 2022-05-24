from tqdm.auto import tqdm
import contextlib
import os
import statsmodels.stats.multitest
import numpy as np


def cell_type_fractions(adata, groupby, groupby_frac):
    return (
        adata.obs.groupby(groupby, observed=True)
        .size()
        .reset_index(name="n_cells")
        .groupby(groupby_frac)
        .apply(lambda x: x.assign(frac_cells=x["n_cells"] / np.sum(x["n_cells"])))
    ).sort_values(groupby)


def log2_fc(
    df, mean_col="intercept", diff_col="coef", key_added="log2_fc", inplace=False
):
    """Add a fold change based on the base-mean and the change in means.
    In a linear model with intercept, the intercept represents the base mean and
    the coefficient the change in means.
    """
    if not inplace:
        df = df.copy()

    def _logfc(mean, diff):
        if mean <= 0 and diff > 0:
            # the intercept may be negative when there is a categorical covariate.
            # We treat it as 0.
            # an increase from 0 -> infinite fold change
            return np.inf
        elif mean + diff <= 0:
            # a decrease to 0 -> -infinite fold change
            return -np.inf
        else:
            return np.log2(diff + mean) - np.log2(mean)

    # The intercept is the mean, the coef the deviation from the mean.
    # Thereby, fold-change = (intercept + coef) / intercept
    df[key_added] = [
        _logfc(mean, diff) for mean, diff in zip(df[mean_col], df[diff_col])
    ]

    if not inplace:
        return df


def fdr_correction(df, pvalue_col="pvalue", *, key_added="fdr", inplace=False):
    """Adjust p-values in a data frame with test results using FDR correction."""
    if not inplace:
        df = df.copy()

    df[key_added] = statsmodels.stats.multitest.fdrcorrection(df[pvalue_col].values)[1]

    if not inplace:
        return df


def split_anndata(adata, groupby):
    """Split an anndata object into a dict of anndata objects based on a column in obs"""
    categories = adata.obs[groupby].unique()
    return {cat: adata[adata.obs[groupby] == cat, :].copy() for cat in tqdm(categories)}


def chunk_adatas(ad, chunksize=200):
    """Generate chunks of adata objects (by variable)"""
    for i in range(0, ad.shape[1], chunksize):
        yield ad[:, i : i + chunksize].copy()


def suppress_stdout(func):
    """Decorator to suppress stdout"""

    def wrapper(*a, **ka):
        with open(os.devnull, "w") as devnull:
            with contextlib.redirect_stdout(devnull):
                return func(*a, **ka)

    return wrapper


def _choose_mtx_rep(adata, use_raw=False, layer=None):
    is_layer = layer is not None
    if use_raw and is_layer:
        raise ValueError(
            "Cannot use expression from both layer and raw. You provided:"
            f"'use_raw={use_raw}' and 'layer={layer}'"
        )
    if is_layer:
        return adata.layers[layer]
    elif use_raw:
        return adata.raw.X
    else:
        return adata.X
