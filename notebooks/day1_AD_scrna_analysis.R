# =============================================================================
# Day 1 — Single Cell RNA-seq Analysis of Alzheimer's Disease
# Dataset: GSE138852 (Grubman et al. 2019, Nature Neuroscience)
# Tissue: Entorhinal Cortex | 6 AD patients vs 6 Controls
# Author: [Your Name]
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1 — Load Libraries
# -----------------------------------------------------------------------------
# Think of libraries like tools in a toolbox
# Seurat   = main scRNA-seq analysis toolkit
# ggplot2  = plotting
# dplyr    = data manipulation (like awk/grep but for dataframes)
# patchwork = combine multiple plots into one figure

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

# Set working directory — all file paths will be relative to this
setwd("/scratch/mallya/bysanir/scrna_project")

# Also set personal R library path (important on HPC)
.libPaths("/scratch/mallya/bysanir/R/library")

# Create output directory for plots
dir.create("results/day1", recursive = TRUE, showWarnings = FALSE)
cat("Output directory created: results/day1/\n")


# -----------------------------------------------------------------------------
# SECTION 2 — Load the Data
# -----------------------------------------------------------------------------
# counts.csv  = gene expression matrix (genes x cells)
#               rows = genes (~10,850)
#               cols = cell barcodes (~13,214)
#               values = UMI counts (how many times each gene was detected)
#
# covariates  = metadata about each cell
#               diagnosis (AD vs Control)
#               cell type (from paper's annotation)
#               donor ID

cat("Loading counts matrix...\n")
counts <- read.csv(
  "data/alzheimers/GSE138852_counts.csv.gz",
  row.names = 1     # first column = gene names, use as row names
)
cat("Counts dimensions (genes x cells):", dim(counts), "\n")
# Expected output: 10850 13214

cat("\nLoading metadata...\n")
meta <- read.csv(
  "data/alzheimers/GSE138852_covariates.csv.gz",
  row.names = 1     # first column = cell barcodes, use as row names
)
cat("Metadata dimensions (cells x columns):", dim(meta), "\n")
cat("Metadata columns:", colnames(meta), "\n")
# Expected output: 13214 5
# Columns: oupSample.batchCond oupSample.cellType oupSample.cellType_batchCond
#          oupSample.subclustID oupSample.subclustCond

# Quick look at metadata
cat("\nFirst 3 rows of metadata:\n")
print(head(meta, 3))

# Verify cell names match between counts and metadata
# Like checking if two files have the same IDs before joining them
cat("\nCell names match between counts and metadata:",
    all(colnames(counts) == rownames(meta)), "\n")
# Expected: TRUE — if FALSE, analysis will fail


# -----------------------------------------------------------------------------
# SECTION 3 — Create Seurat Object
# -----------------------------------------------------------------------------
# A Seurat object is a container that holds:
#   - the count matrix
#   - metadata
#   - all downstream analysis results (PCA, UMAP, clusters etc)
# Think of it like a project directory that holds all related files

cat("\nCreating Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts   = counts,          # the gene x cell count matrix
  meta.data = meta,           # cell metadata
  project  = "AD_EntorhinalCortex",
  min.cells = 3,              # keep genes detected in >= 3 cells
                              # removes extremely rare/noisy genes
  min.features = 200          # keep cells with >= 200 genes detected
                              # removes empty droplets
)

cat("Seurat object created:\n")
print(seurat_obj)
# Expected:
# An object of class Seurat
# 10850 features across 13214 samples within 1 assay
# Active assay: RNA (10850 features, 0 variable features)


# -----------------------------------------------------------------------------
# SECTION 4 — Quality Control (QC)
# -----------------------------------------------------------------------------
# Goal: remove low quality cells before analysis
# Three types of bad cells to remove:
#
# 1. Empty droplets   → very few genes (nFeature < 200)
#    What happened: 10x droplet captured no nucleus, just ambient RNA
#
# 2. Doublets         → too many genes (nFeature > 1500)
#    What happened: two nuclei were captured in one droplet
#    They look like one "super cell" with double the genes
#
# 3. Dying cells      → high mitochondrial % (percent.mt > 5%)
#    What happened: nucleus is damaged, cytoplasmic RNA leaked out
#    Only mitochondrial RNA (more robust) remains
#    NOTE: for snRNA-seq, %mito is naturally very low (<1%)
#    because nuclei don't contain mitochondria

# Calculate % mitochondrial genes per cell
# Human mitochondrial genes all start with "MT-"
# e.g. MT-CO1, MT-CO2, MT-ND1, MT-ATP6 etc
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"    # regex: genes starting with MT-
)

# Check QC metric distributions
cat("\n--- QC Metrics Summary ---\n")
cat("\nnFeature_RNA (genes per cell):\n")
print(summary(seurat_obj$nFeature_RNA))
# Expected: Min~274, Median~646, Max~1632

cat("\nnCount_RNA (UMIs per cell):\n")
print(summary(seurat_obj$nCount_RNA))
# Expected: Min~328, Median~917, Max~2808

cat("\npercent.mt (% mitochondrial):\n")
print(summary(seurat_obj$percent.mt))
# Expected: Min~0, Median~0.24, Max~9.95
# Very low mito % is expected for snRNA-seq (nuclei have no mitochondria)

# Visualize QC metrics as violin plots
# Split by condition (AD vs Control) to check both groups have similar quality
p_qc <- VlnPlot(
  seurat_obj,
  features  = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by  = "oupSample.batchCond",   # AD vs Control
  ncol      = 3,
  pt.size   = 0                        # don't show individual points (too many)
)
ggsave("results/day1/QC_violinplot.png",
       plot = p_qc, width = 12, height = 5, dpi = 150)
cat("\nQC violin plot saved: results/day1/QC_violinplot.png\n")

# Scatter plot: nCount vs nFeature (helps spot doublets and empty droplets)
p_scatter <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "oupSample.batchCond"
)
ggsave("results/day1/QC_scatter.png",
       plot = p_scatter, width = 8, height = 6, dpi = 150)
cat("QC scatter plot saved: results/day1/QC_scatter.png\n")

# Filter cells based on QC thresholds
# Keep cells where ALL three conditions are TRUE
cells_before <- ncol(seurat_obj)
seurat_obj <- subset(seurat_obj,
  subset = nFeature_RNA > 200 &    # remove empty droplets
           nFeature_RNA < 1500 &   # remove doublets
           percent.mt < 5          # remove dying cells
)
cells_after <- ncol(seurat_obj)

cat("\n--- Filtering Results ---\n")
cat("Cells before filtering:", cells_before, "\n")
cat("Cells after filtering: ", cells_after, "\n")
cat("Cells removed:         ", cells_before - cells_after, "\n")
cat("% removed:             ", round((cells_before - cells_after)/cells_before*100, 2), "%\n")
# Expected: 256 cells removed (2%) — low because authors pre-filtered


# -----------------------------------------------------------------------------
# SECTION 5 — Normalization with SCTransform
# -----------------------------------------------------------------------------
# Problem: different cells have different sequencing depths
#   Cell A: 500 total UMIs → GeneX has 10 counts → 10/500 = 2% of reads
#   Cell B: 2000 total UMIs → GeneX has 10 counts → 10/2000 = 0.5% of reads
#   Without correction, GeneX appears 4x higher in Cell A — but it's just
#   a technical artifact of sequencing depth, not real biology
#
# SCTransform solution:
#   1. Models the relationship between UMI count and gene expression
#   2. Corrects for sequencing depth differences
#   3. Selects top 3000 most variable genes (most informative for cell type ID)
#   4. Optionally regresses out confounding variables (like %mito)
#
# SCTransform is better than the old NormalizeData + ScaleData approach
# because it uses a statistical model (negative binomial regression)
# instead of simple scaling

cat("\nRunning SCTransform normalization...\n")
cat("This will take 3-5 minutes...\n")
seurat_obj <- SCTransform(
  seurat_obj,
  vars.to.regress = "percent.mt",  # remove variation due to mito content
  verbose = TRUE
)

cat("SCTransform complete!\n")
cat("Variable features selected:", length(VariableFeatures(seurat_obj)), "\n")
# Expected: 3000 highly variable genes

# Look at top variable genes
top_variable <- head(VariableFeatures(seurat_obj), 20)
cat("\nTop 20 most variable genes:\n")
print(top_variable)


# -----------------------------------------------------------------------------
# SECTION 6 — PCA (Principal Component Analysis)
# -----------------------------------------------------------------------------
# Problem: you have 3000 variable genes = 3000 dimensions
#   Impossible to visualize or cluster in 3000 dimensions
#
# PCA solution:
#   Finds combinations of genes (principal components) that capture
#   the most variation in the data
#   PC1 = direction of most variation
#   PC2 = direction of second most variation, perpendicular to PC1
#   etc.
#
#   Compresses 3000 genes → 50 PCs
#   Most biological signal is in first 20-30 PCs
#   Later PCs are mostly noise
#
# Each PC is a weighted combination of genes
# Genes with high weights = drivers of variation in that PC

cat("\nRunning PCA...\n")
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# Look at which genes drive each PC
cat("\nTop genes driving first 5 PCs:\n")
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)
# PC1: separates oligos (PLP1, MBP) from neurons/astrocytes (DPP10, SLC1A2)
# PC3: picks up microglia (CD74, MEF2C)
# PC5: picks up astrocytes (GFAP)

# Elbow plot — helps decide how many PCs to use downstream
# Use PCs before the "elbow" where variance explained levels off
p_elbow <- ElbowPlot(seurat_obj, ndims = 50)
ggsave("results/day1/PCA_elbow.png",
       plot = p_elbow, width = 8, height = 5, dpi = 150)
cat("Elbow plot saved: results/day1/PCA_elbow.png\n")
cat("Interpretation: use PCs up to where the curve flattens (~PC20-30)\n")

# PCA plot colored by condition
p_pca <- DimPlot(seurat_obj, reduction = "pca",
                 group.by = "oupSample.batchCond") +
  ggtitle("PCA — AD vs Control")
ggsave("results/day1/PCA_plot.png",
       plot = p_pca, width = 8, height = 6, dpi = 150)
cat("PCA plot saved: results/day1/PCA_plot.png\n")


# -----------------------------------------------------------------------------
# SECTION 7 — UMAP
# -----------------------------------------------------------------------------
# Problem: even 50 PCs are too many to visualize
#
# UMAP solution:
#   Non-linear dimensionality reduction
#   Compresses 50 PCs → 2 dimensions (X, Y) for plotting
#   Tries to keep similar cells close together
#   and different cells far apart
#
# Important: UMAP is for visualization only
#   Distances between clusters are not always meaningful
#   Use PCA space for actual computation (clustering, trajectory etc)
#
# dims = 1:30 means use first 30 PCs as input to UMAP
# (matches elbow plot — most signal in first 30 PCs)

cat("\nRunning UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
cat("UMAP complete!\n")


# -----------------------------------------------------------------------------
# SECTION 8 — Clustering
# -----------------------------------------------------------------------------
# Goal: group cells that are transcriptionally similar together
#
# Two steps:
#
# FindNeighbors:
#   For each cell, find its K nearest neighbors in PCA space
#   (which other cells have most similar gene expression?)
#   Builds a graph where similar cells are connected
#
# FindClusters:
#   Uses Louvain algorithm to find groups (communities) in the graph
#   Groups of highly interconnected cells = clusters = likely same cell type
#
# resolution parameter:
#   Low (0.1-0.3)  = fewer, larger clusters
#   High (0.8-1.5) = more, smaller clusters
#   0.3 is good starting point for this dataset size

cat("\nFinding neighbors and clusters...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3, verbose = FALSE)

n_clusters <- length(levels(seurat_obj$seurat_clusters))
cat("Number of clusters found:", n_clusters, "\n")
cat("\nCells per cluster:\n")
print(table(seurat_obj$seurat_clusters))
# Expected: 15 clusters


# -----------------------------------------------------------------------------
# SECTION 9 — Cell Type Annotation
# -----------------------------------------------------------------------------
# Now we have clusters but they are just numbered (0, 1, 2...)
# We need to figure out what cell type each cluster represents
#
# Strategy: compare our clusters to the published cell type labels
# The authors already annotated cell types — stored in oupSample.cellType
# We cross-tabulate our clusters vs their labels
# Whichever cell type dominates a cluster = that cluster's identity

cat("\n--- Cross-tabulation: Our Clusters vs Published Cell Types ---\n")
cross_tab <- table(seurat_obj$seurat_clusters, seurat_obj$oupSample.cellType)
print(cross_tab)
# Reading this table:
# Cluster 0 → dominated by oligo = Oligodendrocyte
# Cluster 3 → dominated by astro = Astrocyte
# Cluster 8 → dominated by mg   = Microglia
# etc.

# Annotate clusters based on cross-tabulation
# RenameIdents changes the cluster labels in the Seurat object
seurat_obj <- RenameIdents(seurat_obj,
  "0"  = "Oligodendrocyte",
  "1"  = "Oligodendrocyte",
  "2"  = "Oligodendrocyte",
  "3"  = "Astrocyte",
  "4"  = "OPC",              # Oligodendrocyte Precursor Cells
  "5"  = "Neuron",
  "6"  = "Oligodendrocyte",
  "7"  = "Astrocyte",
  "8"  = "Microglia",
  "9"  = "Mixed",            # unclear mixed population
  "10" = "Mixed",
  "11" = "Neuron",
  "12" = "Endothelial",
  "13" = "Oligodendrocyte",
  "14" = "Neuron"
)

# Save annotation to metadata column for easy access later
seurat_obj$cell_type_annotated <- Idents(seurat_obj)

cat("\nCell counts per annotated cell type:\n")
print(table(seurat_obj$cell_type_annotated))
# Expected:
# Oligodendrocyte: 7368
# Astrocyte: 2022
# OPC: 1128
# Neuron: 1198
# Microglia: 541
# Mixed: 573
# Endothelial: 128


# -----------------------------------------------------------------------------
# SECTION 10 — Visualization
# -----------------------------------------------------------------------------

cat("\nGenerating UMAP plots...\n")

# Plot 1: numbered clusters
p1 <- DimPlot(seurat_obj,
              reduction = "umap",
              group.by  = "seurat_clusters",
              label     = TRUE,
              pt.size   = 0.3) +
  ggtitle("Seurat Clusters") +
  NoLegend()

# Plot 2: annotated cell types
p2 <- DimPlot(seurat_obj,
              reduction = "umap",
              group.by  = "cell_type_annotated",
              label     = TRUE,
              repel     = TRUE,
              pt.size   = 0.3) +
  ggtitle("Cell Type Annotation")

# Plot 3: AD vs Control
p3 <- DimPlot(seurat_obj,
              reduction = "umap",
              group.by  = "oupSample.batchCond",
              pt.size   = 0.3,
              cols      = c("AD" = "#E41A1C", "ct" = "#377EB8")) +
  ggtitle("AD vs Control")

# Combine and save
combined_umap <- p1 | p2 | p3
ggsave("results/day1/UMAP_overview.png",
       plot = combined_umap, width = 18, height = 6, dpi = 150)
cat("UMAP overview saved: results/day1/UMAP_overview.png\n")

# Plot 4: Annotated UMAP (larger, cleaner)
p_annotated <- DimPlot(seurat_obj,
                        reduction = "umap",
                        group.by  = "cell_type_annotated",
                        label     = TRUE,
                        repel     = TRUE,
                        pt.size   = 0.5) +
  ggtitle("AD Entorhinal Cortex — Cell Types",
          subtitle = "GSE138852 | Grubman et al. 2019 | n=12,958 nuclei") +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 14, face = "bold"))
ggsave("results/day1/UMAP_annotated.png",
       plot = p_annotated, width = 10, height = 8, dpi = 150)
cat("Annotated UMAP saved: results/day1/UMAP_annotated.png\n")


# -----------------------------------------------------------------------------
# SECTION 11 — AD vs Control Cell Type Composition
# -----------------------------------------------------------------------------
# Key biological question: are any cell types enriched or depleted in AD?
# This is called "differential abundance analysis"
#
# We calculate what % of AD cells are each type
# vs what % of Control cells are each type

cat("\n--- Cell Type Composition: AD vs Control ---\n")

# Raw counts
prop_table <- table(seurat_obj$cell_type_annotated,
                    seurat_obj$oupSample.batchCond)
cat("\nRaw cell counts:\n")
print(prop_table)

# Percentages
prop_pct <- prop.table(prop_table, margin = 2) * 100
cat("\nPercentages (% of each condition):\n")
print(round(prop_pct, 2))

# KEY FINDING:
# Astrocytes: 25.7% Control → 5.7% AD  (DEPLETED in AD)
# OPCs:       14.7% Control → 2.8% AD  (DEPLETED in AD)
# Oligos:     41.9% Control → 71.5% AD (ENRICHED in AD)
# This matches Grubman et al. 2019 findings!

# Visualize as stacked bar chart
prop_df <- as.data.frame(prop_pct)
colnames(prop_df) <- c("CellType", "Condition", "Percentage")
prop_df$Condition <- ifelse(prop_df$Condition == "ct", "Control", "AD")
prop_df$Condition <- factor(prop_df$Condition, levels = c("Control", "AD"))

p_prop <- ggplot(prop_df, aes(x = Condition, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title    = "Cell Type Composition: AD vs Control",
    subtitle = "Entorhinal Cortex — GSE138852 (Grubman et al. 2019, Nat Neurosci)",
    x        = "Condition",
    y        = "Percentage of nuclei (%)",
    fill     = "Cell Type"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey40")
  )
ggsave("results/day1/celltype_proportions.png",
       plot = p_prop, width = 8, height = 6, dpi = 150)
cat("Proportion plot saved: results/day1/celltype_proportions.png\n")


# -----------------------------------------------------------------------------
# SECTION 12 — Save Seurat Object
# -----------------------------------------------------------------------------
# Save the complete Seurat object so you don't have to rerun everything
# on Day 2 — just load this file and continue
# Like a checkpoint in a pipeline

cat("\nSaving Seurat object...\n")
saveRDS(seurat_obj, "data/alzheimers/seurat_day1.rds")
cat("Seurat object saved: data/alzheimers/seurat_day1.rds\n")

# How to reload on Day 2:
# seurat_obj <- readRDS("data/alzheimers/seurat_day1.rds")


# -----------------------------------------------------------------------------
# SECTION 13 — Session Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("=============================================================\n")
cat("DAY 1 COMPLETE — Summary\n")
cat("=============================================================\n")
cat("Dataset:          GSE138852 (Grubman et al. 2019, Nat Neurosci)\n")
cat("Tissue:           Entorhinal Cortex\n")
cat("Donors:           6 AD + 6 Control\n")
cat("Nuclei loaded:    13,214\n")
cat("Nuclei after QC:  12,958 (removed 256 low quality)\n")
cat("Genes analyzed:   10,850\n")
cat("Clusters found:   15\n")
cat("Cell types ID'd:  7 (Oligo, Astro, OPC, Neuron, Microglia, Endo, Mixed)\n")
cat("\nKey biological finding:\n")
cat("  Astrocytes: 25.7% (Control) → 5.7% (AD)  [DEPLETED]\n")
cat("  OPCs:       14.7% (Control) → 2.8% (AD)  [DEPLETED]\n")
cat("  Oligos:     41.9% (Control) → 71.5% (AD) [ENRICHED]\n")
cat("\nOutputs saved to: results/day1/\n")
cat("  QC_violinplot.png\n")
cat("  QC_scatter.png\n")
cat("  PCA_elbow.png\n")
cat("  PCA_plot.png\n")
cat("  UMAP_overview.png\n")
cat("  UMAP_annotated.png\n")
cat("  celltype_proportions.png\n")
cat("\nSeurat object: data/alzheimers/seurat_day1.rds\n")
cat("=============================================================\n")

# Print session info for reproducibility
cat("\nR Session Info:\n")
sessionInfo()
