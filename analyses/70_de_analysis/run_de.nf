#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process DESEQ2 {
    tag { meta.id }
    publishDir "../../data/70_de_analysis/72_run_de/", mode: "copy"

    cpus { meta['singlecell'] ? 11 : 4 }
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-de2"
    errorStrategy 'finish'

    input:
        tuple val(meta), path(expr), path(samplesheet)

    output:
        path("deseq2_res*")

    script:
    def sc_options = meta['singlecell'] ? "--fc_cutoff=0.58 --fdr_cutoff=0.000001 --skip_pca" : ""
    def covariates = meta['covariates'] ? "--covariate_formula=${meta['covariate_formula']}" : ""
    def sc_rm = meta['singlecell'] ? "rm deseq2_res_sc_**/*detectedGenes*_min_10_reads*.tsv" : ""
    def paired = meta['paired_grp'] ? "--paired_grp=${meta['paired_grp']}" : ""
    """
    mkdir deseq2_res_${meta.id}
    runDESeq2_ICBI.R $samplesheet $expr \\
        --result_dir=deseq2_res_${meta.id} \\
        --c1=${meta.c1} \\
        --c2=${meta.c2} \\
        --sample_col=${meta.sample_col} \\
        --condition_col=${meta.condition_col} \\
        --gene_id_type=SYMBOL \\
        --n_cpus=${task.cpus} \\
        ${paired} \\
        ${sc_options} \\
        ${covariates}

    ${sc_rm}
    """
}

workflow {
    def input_path = "../../data/70_de_analysis/71_make_pseudobulk/"
    DESEQ2(
        Channel.from(
            [
                [singlecell: false, id: "bulk_response_all_timepoints", c1: "short_term", c2: "long_term", sample_col: "patient", condition_col: "response"],
                [singlecell: false, id: "bulk_response_t0", c1: "short_term", c2: "long_term", sample_col: "patient", condition_col: "response"],
                [singlecell: false, id: "bulk_timepoints", c1: "post-treatment", c2: "pre-treatment", sample_col: "sample", condition_col: "timepoint", paired_grp: "patient"],
                [singlecell: true, id: "sc_response_all_timepoints_P1_vs_long_term", c1: "short_term", c2: "long_term", sample_col: "cell_id", condition_col: "response"],
                [singlecell: true, id: "sc_response_all_timepoints_P2_vs_long_term", c1: "short_term", c2: "long_term", sample_col: "cell_id", condition_col: "response"],
                [singlecell: true, id: "sc_response_all_timepoints_P3_vs_long_term", c1: "short_term", c2: "long_term", sample_col: "cell_id", condition_col: "response"],
                [singlecell: true, id: "sc_response_T0_P2_vs_long_term", c1: "short_term", c2: "long_term", sample_col: "cell_id", condition_col: "response"],
                [singlecell: true, id: "sc_response_T0_P3_vs_long_term", c1: "short_term", c2: "long_term", sample_col: "cell_id", condition_col: "response"],
                [singlecell: true, id: "sc_healthy_vs_malignant_b_cells", c1: "malignant", c2: "healthy", sample_col: "cell_id", condition_col: "cell_type", covariate_formula: "+patient"],
                [singlecell: true, id: "sc_fosb_b_cells", c1: "fos_malignant_b", c2: "malignant_b", sample_col: "cell_id", condition_col: "cell_phenotype", covariate_formula: "+patient+timepoint"],
            ]
        ).map {
            it -> [
                it,
                file(input_path + it['id'] + (it['singlecell'] ? "" : "_bulk_df") + ".tsv", checkIfExists: true),
                file(input_path + it['id'] + "_samplesheet.csv", checkIfExists: true)
            ]
        }
    )
}
