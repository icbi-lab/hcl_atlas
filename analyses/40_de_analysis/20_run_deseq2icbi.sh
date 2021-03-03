#!/bin/bash

mkdir -p ../../data/40_de_analysis/deseq2_res

../../lib/deseq2icbi/runDESeq2_ICBI.R ../../data/40_de_analysis/deseq_samplesheet.csv ../../data/40_de_analysis/bulk_df.tsv \
    --result_dir=../../data/40_de_analysis/deseq2_res \
    --c1=short_term \
    --c2=long_term \
    --sample_col=sample \
    --condition_col=response \
    --gene_id_type=SYMBOL