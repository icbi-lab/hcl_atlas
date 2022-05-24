import numpy as np
import scanpy as sc


def annotate_cell_types(
    adata,
    cell_type_map,
    *,
    key_added="cell_type",
    default_cell_type="other",
    column="leiden",
):
    """Generate a cell-type column from a Mapping cell_type -> [list of clusters]"""
    res = np.full(adata.shape[0], default_cell_type, dtype=object)
    for ct, clusters in cell_type_map.items():
        clusters = [str(x) for x in clusters]
        res[adata.obs[column].isin(clusters)] = ct

    adata.obs[key_added] = res
    sc.pl.umap(adata, color=key_added)
