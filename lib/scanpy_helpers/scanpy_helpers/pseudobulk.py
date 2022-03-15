from typing import Sequence, Union
from anndata import AnnData, ImplicitModificationWarning
import numpy as np
import pandas as pd
from operator import and_
from functools import reduce
import warnings
from .util import _choose_mtx_rep


def pseudobulk(
    adata,
    *,
    groupby: Union[str, Sequence[str]],
    aggr_fun=np.sum,
    min_obs=10,
    layer=None,
    use_raw=None,
) -> AnnData:
    """
    Calculate Pseudobulk of groups

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        One or multiple columns to group by
    aggr_fun
        Callback function to calculate pseudobulk. Must be a numpy ufunc supporting
        the `axis` attribute.
    min_obs
        Exclude groups with less than `min_obs` observations

    Returns
    -------
    New anndata object with same vars as input, but reduced number of obs.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    combinations = adata.obs.loc[:, groupby].drop_duplicates()

    X = _choose_mtx_rep(adata, use_raw, layer)

    # precompute masks
    masks = {}
    for col in groupby:
        masks[col] = {}
        for val in combinations[col].unique():
            masks[col][val] = adata.obs[col] == val

    expr_agg = []
    obs = []

    for comb in combinations.itertuples(index=False):
        mask = reduce(and_, (masks[col][val] for col, val in zip(groupby, comb)))
        if np.sum(mask) < min_obs:
            continue
        expr_row = aggr_fun(X[mask, :], axis=0)
        obs_row = comb._asdict()
        obs_row["n_obs"] = np.sum(mask)
        # convert matrix to array if required (happens when aggregating spares matrix)
        try:
            expr_row = expr_row.A1
        except AttributeError:
            pass
        obs.append(obs_row)
        expr_agg.append(expr_row)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ImplicitModificationWarning)
        return AnnData(
            X=np.vstack(expr_agg),
            var=adata.var,
            obs=pd.DataFrame.from_records(obs),
        )


def write_pseudobulk(pb, prefix):
    """Write pseudobulk samplesheet and expression in the form required for the deseq2icbi script"""
    pb.obs.index.name = "sample_id"
    pb.obs.to_csv(f"{prefix}_samplesheet.csv", index=True)
    expr_df = pd.DataFrame(
        pb.X.T,
        columns=pb.obs_names,
        index=pb.var_names,
    )
    expr_df.index.name = "gene_id"
    expr_df.insert(0, "gene_name", expr_df.index)
    expr_df.to_csv(f"{prefix}_bulk_df.tsv", sep="\t")
