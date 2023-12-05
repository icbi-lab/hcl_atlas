# HCL atlas

This repository contains the source code to analyze the single-cell data from Hairy Cell Leukemia published in 

> Jan-Paul Bohn et al., Single-cell RNA-Seq-based deconvolution of hairy cell leukemia reveals novel disease drivers and identifies DUSP1 as potential therapeutic target, *Submitted*

## Description of analysis steps

* **10_prepare_adata**: Load BD Rhapsody WTA analysis pipeline outputs into AnnData objects and add metadata.
* **20_scrnaseq_qc**: Use a nextflow pipeline (stored in `lib/single-cell-analysis-nf`) to perform threshold-based filtering of single-cell data and apply `SOLO` for doublet detection.
* **30_merge_adata**: Merge samplese into a single AnnData object, train a `scVI` model for batch effect removal, and annotate cell-types based on unsupervised clustering
* **40_cluster_analysis**: Identify and investigate subclusters representing cell-states that go beyond the major cell-types
* **50_de_analysis**: Generate pseudobulk and perform differential gene expression analysis using DESeq2 (based on a wrapper script stored in `lib/deseq2_workflow`)
* **70_downstream_analysis**: Perform pathway analyses and generate figures for publication based on the data generated in the previous steps

## Data availability

Preprocessed data (UMI counts) as well as intermediate and final results are available from Zenodo: https://zenodo.org/doi/10.5281/zenodo.10262912
For each analysis above below, there is one output folder provided as zip file on zenodo. The output of previous steps is used 
as input for subsequent steps. 

Additionally, we provide the conda environments used for the analysis as singularity containers (also on Zenodo). 

## Contact
For any requests regarding the single-cell data analysis, please use the [issue tracker](https://github.com/icbi-lab/hairy_cell_leukemia/issues).
For anything else, you can reach out to the corresponding authors listed in the manuscript. 

## Notes on reproducibility
We aim at making the analysis reproducible by providing the input data and all software dependencies as a singularity container. 
However, some analysis steps inherently produce different results when ran on different hardware. See also https://github.com/scverse/scanpy/issues/2014
for a discussion. 

This particularly affects the generation of the neighborhood graph, which impacts unsupervised clustering and therefore cell-type annotation. 
Instead of running the entire analysis from scratch, you can also start from some of the intermediate results including the cell-type
annotation we provide on Zenodo. 
