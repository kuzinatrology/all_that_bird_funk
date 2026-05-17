## Pipeline Requirements

To run the pipeline, the following input data are required:

1. CDS sequences, complete genomes, genome annotation files (`.gtf`), and proteomes for bird species.
2. Life-history traits of the species, including body mass and maximum lifespan (MLS).
3. Amino acid sequences of the genes of interest for which homologs will be searched. In this pipeline, the OpenGenes gene set is used.
4. A species tree.

## Analysis

1. Homologous sequence search — `cicl_mmseq_easy_search.sh`
2. Correlation analysis between CpG density and longevity quotient (LQ) — `cpg_density_count.ipynb` and `pgls.sh`
3. dN/dS analysis — `dn_ds_pipeline.ipynb` and the corresponding `.sh` files referenced in the pipeline
