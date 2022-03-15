#!/usr/bin/env python
"""Compute the correlation (coexpression) of a certain gene with all others.

Usage:
    gene_correlation.py ADATA GENE CORES OUTPUT_FILE

e.g.
    gene_correlation.py adata.h5ad KLRB1 32 ../correlation.xlsx
"""
import scanpy as sc
import sys
import warnings
from multiprocessing import Pool
import scipy.stats
from tqdm import tqdm
import pandas as pd
import numpy as np

warnings.filterwarnings("ignore", category=scipy.stats.PearsonRConstantInputWarning)

adata_path = sys.argv[1]
query_gene = sys.argv[2].strip()
cpus = int(sys.argv[3])
output_file = sys.argv[4]

adata = sc.read_h5ad(adata_path)


def get_cor(gene):
    return scipy.stats.pearsonr(
        adata[:, query_gene].X.toarray()[:, 0], adata[:, gene].X.toarray()[:, 0]
    )


with Pool(cpus) as p:
    res = list(
        tqdm(p.imap(get_cor, adata.var_names, chunksize=50), total=adata.shape[1])
    )


res_df = pd.DataFrame().from_records(res)
res_df.index = adata.var_names
res_df.columns = ["pearson_r", "pvalue"]
res_df["log_p"] = -np.log10(res_df["pvalue"])
res_df.sort_values("pearson_r", ascending=False, inplace=True)

res_df.to_csv(output_file)
