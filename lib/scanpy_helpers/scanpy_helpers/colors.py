import altair as alt
import pandas as pd


def set_scale_anndata(adata, column, palette=None):
    if palette is None:
        palette = column

    adata._sanitize()

    tmp_cols = getattr(COLORS, palette)
    adata.uns[f"{column}_colors"] = [
        tmp_cols[cat] for cat in adata.obs[column].cat.categories
    ]


def altair_scale(variable):
    return alt.Scale(
        domain=list(getattr(COLORS, variable).keys()),
        range=list(getattr(COLORS, variable).values()),
    )


def plot_palette(variable):
    """Display a palette"""
    tmp_cols = getattr(COLORS, variable)
    return (
        alt.Chart(
            pd.DataFrame.from_dict(tmp_cols, orient="index", columns=["color"])
            .reset_index()
            .rename(columns={"index": variable})
        )
        .mark_rect(height=40, width=30)
        .encode(
            x=alt.X(variable),
            color=alt.Color(variable, scale=altair_scale(variable), legend=None),
        )
    )


def plot_all_palettes():
    return alt.vconcat(
        *[plot_palette(v) for v in COLORS.__dict__.keys() if not v.startswith("_")]
    ).resolve_scale(color="independent")


class COLORS:
    timepoint = {"T0": "#edf8b1", "T1": "#7fcdbb", "T2": "#2c7fb8"}
    # okabe palette
    patient = {
        "P1": "#E69F00",
        "P2": "#56B4E9",
        "P3": "#009E73",
        "P4": "#0072B2",
        "P5": "#D55E00",
        "P6": "#CC79A7",
    }
    # Dark2
    cell_type = {
        "healthy B cell": "#1b9e77",
        "malignant B cell": "#d95f02",
        "malignant B cell (dividing)": "#e7298a",
        "T cell": "#7570b3",
        "Plasma cell": "#66a61e",
    }
    response = {"short_term": "#377eb8", "long_term": "#4daf4a"}

    # sex = {
    #     "male": "#80b1d3",
    #     "female": "#fccde5",
    #     "unknown": "#dddddd",
    # }
    pass
