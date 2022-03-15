#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SPLIT_ADATA {
	cpus 1
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-scanpy"

    input:
        tuple val(id), path(adata)

    output:
        path("*.h5ad"), emit: adata

    script:
    """
    #!/usr/bin/env python

    import scanpy as sc

    adata = sc.read_h5ad("$adata")

    for p in adata.obs["patient"].unique():
        tmp_adata = adata[(adata.obs["patient"] == p) & (adata.obs["cell_type"] == "malignant B cell"), :].copy()
        tmp_adata.write_h5ad(f"adata_malignant_b_{p}.h5ad")
    """
}

process RUN_COEXPRESSION {
    cpus 22
    conda "/data/scratch/sturm/conda/envs/2021-hairy-cell-leukemia-wolf-scanpy"
    publishDir "../../data/70_de_analysis/78_run_coexpression", mode: "copy"

    input:
        path(adata)

    output:
        path("*.csv"), emit: coexpression

    script:
    """
    gene_correlation.py $adata LGR5 ${task.cpus} ${adata.baseName}.csv
    """
}

workflow{
    SPLIT_ADATA(["adata_scvi", file("../../data/30_merge_adata/adata_scvi.h5ad")])
    RUN_COEXPRESSION(SPLIT_ADATA.out.adata.flatten())
}
