# Single-Cell and Spatial Transcriptomics Analysis of Alzheimers Disease

Dataset: GSE138852 (Grubman et al. 2019, Nature Neuroscience)
Tissue: Entorhinal Cortex | 6 AD patients + 6 Controls | 12958 nuclei
Spatial: 10x Visium Mouse Brain Anterior | 2611 spots | 32285 genes
Tools: R, Seurat v5, ggplot2, dplyr
Environment: Linux HPC SLURM

## Key Findings

- Astrocytes depleted 25.7% to 5.7% in AD
- OPCs depleted 14.7% to 2.8% in AD
- 507 DE genes in astrocytes: C3 CD44 LINGO1 up; HES4 HES5 WIF1 down
- 12 genes validated by both DE and pseudotime analyses
- 10 brain regions mapped from 10x Visium spatial data

## Pipeline

- Day 1: QC normalization clustering cell type annotation
- Day 2: Marker genes differential expression AD vs Control
- Day 3: Spatial transcriptomics spatially variable genes region mapping
- Day 4: Pseudotime trajectory cross-dataset validation

## Reference

Grubman et al. 2019. Nature Neuroscience. https://doi.org/10.1038/s41593-019-0539-4
