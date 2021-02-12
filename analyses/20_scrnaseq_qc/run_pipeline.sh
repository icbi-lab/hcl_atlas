#!/usr/bin/bash

cd ../.. && \
nextflow run /home/sturm/projects/2020/single-cell-analysis-nf/main.nf \
    --input tables/samplesheet_scrnaseq_preprocessing.csv \
    --outdir data/20_scrnaseq_qc \
    -resume \
    -profile icbi_long \
    --skip_solo \
    -w /data/scratch/sturm/projects/2021/hairy_cell_leukemia_wolf/work

