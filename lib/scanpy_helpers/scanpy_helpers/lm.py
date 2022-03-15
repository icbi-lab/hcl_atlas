"""High-level wrappers around linear models to find differences between groups. """


from typing import List, Mapping, Optional, Sequence, Tuple, Literal
import pandas as pd
from tqdm.auto import tqdm
import statsmodels.formula.api as smf
from tqdm.contrib.concurrent import process_map
import re
from anndata import AnnData
from .pseudobulk import pseudobulk
import numpy as np
import pandas as pd
from .util import chunk_adatas
import statsmodels.stats
import itertools
from multiprocessing import cpu_count
from threadpoolctl import threadpool_limits
import altair as alt


def plot_lm_result_altair(
    df,
    p_cutoff=0.1,
    p_col="fdr",
    x="variable",
    y="group",
    color="coef",
    title="heatmap",
):
    df_filtered = df.loc[lambda _: _[p_col] < p_cutoff, :]
    df_subset = df.loc[
        lambda _: _[x].isin(df_filtered[x].unique()) & _[y].isin(df[y].unique())
    ]
    if not df_subset.shape[0]:
        print("No values to plot")
        return

    def _get_significance(fdr):
        if fdr < 0.001:
            return "< 0.001"
        elif fdr < 0.01:
            return "< 0.01"
        elif fdr < 0.1:
            return "< 0.1"
        else:
            return np.nan

    df_subset["FDR"] = pd.Categorical([_get_significance(x) for x in df_subset["fdr"]])

    value_max = max(abs(np.nanmin(df_subset[color])), abs(np.nanmax(df_subset[color])))
    return (
        alt.Chart(df_subset, title=title)
        .mark_rect()
        .encode(
            x=x,
            y=y,
            color=alt.Color(
                color,
                scale=alt.Scale(
                    scheme="redblue",
                    reverse=True,
                    domain=[-value_max, value_max],
                ),
            ),
        )
        + alt.Chart(df_subset.loc[lambda x: ~x["FDR"].isnull()])
        .mark_point(color="white", filled=True, stroke="black", strokeWidth=0)
        .encode(
            x=x,
            y=y,
            size=alt.Size(
                "FDR:N",
                scale=alt.Scale(
                    domain=["< 0.001", "< 0.01", "< 0.1"],
                    range=4 ** np.array([3, 2, 1]),
                ),
            ),
        )
    ).configure_mark(opacity=1)


def lm_test_all(
    adatas: Mapping[str, AnnData],
    *,
    groupby: List[str],
    column_to_test: str,
    contrasts: str = "Sum",
    lm_covariate_str: str = "",
    min_categories: int = 2,
) -> Optional[pd.DataFrame]:
    """High-level function to compare groups multiple single cell objects.

    Performs pseudobulk aggregation and a statistical test with linear model

    Parameters
    ----------
    adatas
        Dictionary id -> adata
    groupby
        Grouping variables, will be used to generate pseudobulk (e.g. ["patient", "dataset"])
    column_to_test
        column in obs that contains the variables to test. Must be the same within each group.
        (e.g. a property that applies to all cells of a patient)
    contrasts
        Which contrast-coding to use. T = treatment-coding, S =sum-to-zero
        See https://www.statsmodels.org/devel/contrasts.html
    lm_covariate_str
        Covariate string that is literally included in the linear model formula.
        E.g. `+ dataset`.
    min_categories
        Only perform test if there are at least `min_categories` unique entries in `column_to_test`.

    Returns
    -------
    Data frame with coefficients and pvalues.
    """
    res_list = []
    for ct, tmp_adata in tqdm(adatas.items()):
        tmp_bdata = pseudobulk(
            tmp_adata,
            groupby=groupby + [column_to_test],
            aggr_fun=np.mean,
        )
        if tmp_bdata.obs[column_to_test].nunique() >= min_categories:
            tmp_res = test_lm(
                tmp_bdata,
                f"~ C({column_to_test}, {contrasts}) {lm_covariate_str}",
                column_to_test,
                contrasts=contrasts,
                progress=False,
            )
            res_list.append(tmp_res.assign(cell_type=ct))

    return (
        pd.concat(res_list)
        # .dropna(how="any")
        .assign(
            fdr=lambda x: statsmodels.stats.multitest.fdrcorrection(x["pvalue"].values)[
                1
            ]
        ).sort_values("pvalue")
    )


def test_lm(
    pseudobulk: AnnData,
    formula: str,
    groupby: str,
    *,
    contrasts: str = "Sum",
    progress: bool = True,
    n_jobs: int = None,
    chunksize=200,
):
    """
    Use a linear model to find differences between groups

    In this case we use sum-to-zero or deviation coding to find
    deviations from the mean of means

    Parameters
    ----------
    pseudobulk
        AnnData object with pseudobulk
    formula
        Formula for the linear model. Currently MUST have the following structure
        `~ 0 + C({column}, Sum) + ...`.
    groupby
        Column in adata that contains the groups
    contrasts
        Which contrast-coding to use. T = treatment-coding, S =sum-to-zero
        See https://www.statsmodels.org/devel/contrasts.html
    progress
        Show the tqdm progress bar
    n_jobs
        Run test in parallel. Set to 1 to disable parallelism.
    chunksize
        Splits up anndata in chunks of var and processes chunks in parallel.

    Returns
    -------
    lms
        List of linear models
    df
        Pandas data frame with coefficients and pvalues

    """
    if n_jobs is None:
        n_jobs = cpu_count()

    # likely overhead is larger until several times chunksize, but not tested in detail.
    if pseudobulk.shape[1] < chunksize * 2 or n_jobs == 1:
        return _test_lm(
            pseudobulk, formula, groupby, contrasts=contrasts, progress=progress
        )
    else:
        return pd.concat(
            process_map(
                _test_lm,
                list(chunk_adatas(pseudobulk, chunksize=chunksize)),
                itertools.repeat(formula),
                itertools.repeat(groupby),
                itertools.repeat(contrasts),
                itertools.repeat(False),
                max_workers=n_jobs,
            )
        )


def _test_lm(
    pseudobulk: AnnData,
    formula: str,
    groupby: str,
    contrasts: str,
    progress: bool = False,
) -> pd.DataFrame:

    var_names = pseudobulk.var_names
    formula = "Q('{col}') " + formula

    df = pseudobulk.obs.join(
        pd.DataFrame(pseudobulk.X, index=pseudobulk.obs_names, columns=var_names)
    )

    pseudobulk.obs[groupby] = pd.Categorical(pseudobulk.obs[groupby])
    all_groups = pseudobulk.obs[groupby].unique()

    def test_all_params(res, all_groups):
        """If we use sum-to-zero coding, perform the test for the omitted factor"""
        # only using the categories gets rid of NAN
        if contrasts == "S":
            groups_to_test = all_groups.categories.values.tolist()
        else:
            groups_to_test = [
                g for g in all_groups.categories.values.tolist() if g not in contrasts
            ]
        assert "nan" not in groups_to_test
        # print(all_groups)
        # print(res.summary())
        keys = [
            f"C({groupby}, {contrasts})[{contrasts[0]}.{g}]" for g in groups_to_test
        ]

        intercept = res.params["Intercept"]

        if contrasts.startswith("Sum"):
            # print(keys)
            # print(res.params)
            coefs = res.params[keys[:-1]].to_dict()
            pvals = res.pvalues[keys[:-1]].to_dict()
            # test the level that was omitted for redundancy
            coefs[keys[-1]] = -sum(coefs.values())
            pvals[keys[-1]] = float(
                res.f_test(" + ".join([f"{k}" for k in keys[:-1]]) + " = 0").pvalue
            )
        else:
            coefs = res.params[keys].to_dict()
            pvals = res.pvalues[keys].to_dict()

        return coefs, pvals, intercept

    results: List[pd.DataFrame] = []
    lms = []
    var_iter = tqdm(var_names) if progress else var_names
    for col in var_iter:
        try:
            with threadpool_limits(1):
                mod = smf.ols(formula=formula.format(col=col), data=df)
                res = mod.fit()
                coefs, pvals, intercept = test_all_params(res, all_groups)
                res_df = (
                    pd.DataFrame.from_dict(coefs, orient="index", columns=["coef"])
                    .assign(intercept=intercept)
                    .join(
                        pd.DataFrame.from_dict(
                            pvals, orient="index", columns=["pvalue"]
                        )
                    )
                    .assign(
                        variable=col,
                        group=lambda x: [
                            re.search(f"\[{contrasts[0]}\.(.*)\]", k).groups()[0]
                            for k in x.index
                        ],
                    )
                )
            results.append(res_df)
            lms.append(res)
        except ValueError:
            pass

    return pd.concat(results)
