#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process DESEQ2 {
    tag { meta.id }
    publishDir "../../data/50_de_analysis/deseq2_results", mode: "copy"

    cpus { meta['singlecell'] ? 11 : 4 }
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-de2"
    // container "${projectDir}/../../data/containers/2021-hairy-cell-leukemia-wolf-de2.sif"
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
        --c1="${meta.c1}" \\
        --c2="${meta.c2}" \\
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
    def input_path = "../../data/50_de_analysis/pseudobulk/"
    DESEQ2(
        Channel.from(
            [
                [singlecell: false, id: "bulk_response_all_timepoints", c1: "short_term", c2: "long_term", sample_col: "sample_id", condition_col: "response"],
                [singlecell: false, id: "bulk_response_t0", c1: "short_term", c2: "long_term", sample_col: "sample_id", condition_col: "response"],
                [singlecell: false, id: "bulk_timepoints", c1: "post-treatment", c2: "pre-treatment", sample_col: "sample_id", condition_col: "timepoint", paired_grp: "patient"],
                [singlecell: false, id: "bulk_fos_jun_vs_rest", c1: "DUSP1ʰᶦ╱FOSBʰᶦ╱JUNDʰᶦ", c2: "DUSP1ˡᵒ╱FOSBˡᵒ╱JUNDˡᵒ", sample_col: "sample_id", condition_col: "cell_phenotype", paired_grp: "patient"],
                [singlecell: false, id: "bulk_healthy_malignant", c1: "healthy B cell", c2: "HCL cell", sample_col: "sample_id", condition_col: "cell_type", paired_grp: "patient"],
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
