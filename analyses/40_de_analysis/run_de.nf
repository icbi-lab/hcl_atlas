#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process DE_BULK {
    publishDir "../../data/40_de_analysis/20_run_de/", mode: "copy"

    cpus 4
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-de2"

    input:
        tuple val(id), path(files)

    output:
        path("deseq2_res*")

    script:
    """
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    mkdir deseq2_res_${id}
    runDESeq2_ICBI.R ${files.get(1)} ${files.get(0)} \\
        --result_dir=deseq2_res_${id} \\
        --c1=short_term \\
        --c2=long_term \\
        --sample_col=patient \\
        --condition_col=response \\
        --gene_id_type=SYMBOL \\
        --n_cpus=${task.cpus}
    """
}

process DE_SINGLE_CELL {
    publishDir "../../data/40_de_analysis/20_run_de/", mode: "copy"

    errorStrategy 'ignore'
    cpus 11
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-de2"

    input:
        tuple val(id), path(files)

    output:
        path("deseq2_res*")

    script:
    """
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus}  \\
        MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus}  \\
        MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus}
    mkdir deseq2_res_${id}
    runDESeq2_ICBI.R ${files.get(1)} ${files.get(0)} \\
        --result_dir=deseq2_res_${id} \\
        --c1=short_term \\
        --c2=long_term \\
        --sample_col=cell_id \\
        --condition_col=response \\
        --gene_id_type=SYMBOL \\
        --fc_cutoff=0.58 \\
        --fdr_cutoff=0.000001 \\
        --skip_pca
    """
}

workflow {
    DE_BULK(
        Channel.fromFilePairs("../../data/40_de_analysis/10_make_pseudobulk/bulk*{bulk_df.tsv,samplesheet.csv}", checkIfExists:true),
    )
    DE_SINGLE_CELL(
        Channel.fromFilePairs("../../data/40_de_analysis/10_make_pseudobulk/sc*{.tsv,samplesheet.csv}", checkIfExists:true)
    )
}
