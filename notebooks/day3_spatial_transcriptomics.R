# =============================================================================
# Day 3 — Spatial Transcriptomics: Mouse Brain Visium Analysis
# Dataset: 10x Visium Mouse Brain Anterior Section
# Tools: Seurat spatial functions
# Author: [Your Name]
#
# Day 1 answered: WHICH cell types change in AD? (composition)
# Day 2 answered: WHAT genes change in those cells? (expression)
# Day 3 answers:  WHERE are cell types located in brain tissue? (space)
#
# Key concept: Visium captures gene expression at spatial coordinates
# Each "spot" = ~10-20 cells at a specific X,Y position on tissue
# You get gene expression AND physical location simultaneously
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1 — Load Libraries
# -----------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/scratch/mallya/bysanir/scrna_project")
.libPaths("/scratch/mallya/bysanir/R/library")

# Increase memory limit for large matrix operations
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

dir.create("results/day3", recursive = TRUE, showWarnings = FALSE)


# -----------------------------------------------------------------------------
# SECTION 2 — Load Visium Data
# -----------------------------------------------------------------------------
# Visium data has two components:
#
# 1. Count matrix (like scRNA-seq):
#    filtered_feature_bc_matrix/
#    ├── matrix.mtx.gz    → sparse gene x spot matrix
#    ├── barcodes.tsv.gz  → spot barcodes
#    └── features.tsv.gz  → gene names
#
# 2. Spatial information (unique to Visium):
#    spatial/
#    ├── tissue_positions_list.csv → X,Y coordinates of each spot
#    ├── tissue_hires_image.png    → H&E stained tissue image
#    ├── tissue_lowres_image.png   → lower resolution image
#    └── scalefactors_json.json    → pixel scaling factors
#
# Key difference from scRNA-seq:
#   scRNA-seq: cells have gene expression but NO physical location
#   Visium:    spots have gene expression AND X,Y coordinates on tissue

cat("Reading Visium count matrix...\n")
counts <- Read10X(
  data.dir = "data/visium_mouse/filtered_feature_bc_matrix"
)
cat("Matrix dimensions (genes x spots):", dim(counts), "\n")
# Expected: 32285 genes x 2695 spots

# Create Seurat object with Spatial assay
brain <- CreateSeuratObject(
  counts  = counts,
  assay   = "Spatial",   # important: label as Spatial not RNA
  project = "VisiumMouseBrain"
)
cat("Seurat object created!\n")
print(brain)


# -----------------------------------------------------------------------------
# SECTION 3 — Add Spatial Coordinates
# -----------------------------------------------------------------------------
# tissue_positions_list.csv contains:
# Column 1: barcode         → spot ID (matches count matrix column names)
# Column 2: in_tissue       → 1 if spot is on tissue, 0 if off
# Column 3: array_row       → row position in Visium array grid
# Column 4: array_col       → column position in Visium array grid
# Column 5: pxl_row_in_fullres → Y pixel coordinate in full resolution image
# Column 6: pxl_col_in_fullres → X pixel coordinate in full resolution image
#
# The pixel coordinates are what you use to plot spots on the tissue image

cat("\nLoading spatial coordinates...\n")
coords <- read.csv(
  "data/visium_mouse/spatial/tissue_positions_list.csv",
  header = FALSE
)
colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres")

# Keep only spots that landed on tissue
coords <- coords[coords$in_tissue == 1, ]
rownames(coords) <- coords$barcode
cat("Tissue spots:", nrow(coords), "\n")

# Match spots between count matrix and coordinates
common_spots <- intersect(colnames(brain), rownames(coords))
cat("Common spots:", length(common_spots), "\n")
# If these numbers differ, some spots were filtered by min.cells threshold

brain  <- brain[, common_spots]
coords <- coords[common_spots, ]

# Add coordinates to Seurat metadata
# These will be used for all spatial plots
brain$array_row <- coords$array_row
brain$array_col <- coords$array_col
brain$pxl_row   <- coords$pxl_row_in_fullres
brain$pxl_col   <- coords$pxl_col_in_fullres

cat("Coordinates added to metadata!\n")
cat("Preview of metadata:\n")
print(brain@meta.data[1:3, ])


# -----------------------------------------------------------------------------
# SECTION 4 — Quality Control
# -----------------------------------------------------------------------------
# Same concept as Day 1 but numbers are very different:
#
# scRNA-seq (Day 1):        Visium spatial (Day 3):
# ─────────────────         ───────────────────────
# 1 cell per barcode        ~10-20 cells per spot
# Median ~917 UMIs          Median ~25,888 UMIs
# Median ~646 genes         Median ~6,228 genes
# %mito ~0.24%              %mito ~13%
#
# Higher counts because Visium captures multiple cells per spot
# Higher %mito because whole cells (not just nuclei) are captured
# Note: mouse mito genes = "mt-" (lowercase), human = "MT-" (uppercase)

brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^mt-")

cat("\nQC Metrics Summary:\n")
cat("\nnCount_Spatial (UMIs per spot):\n")
print(summary(brain$nCount_Spatial))

cat("\nnFeature_Spatial (genes per spot):\n")
print(summary(brain$nFeature_Spatial))

cat("\npercent.mt:\n")
print(summary(brain$percent.mt))

# Violin plots
p_qc <- VlnPlot(
  brain,
  features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"),
  pt.size  = 0.1,
  ncol     = 3
)
ggsave("results/day3/QC_spatial_violin.png",
       plot = p_qc, width = 12, height = 5, dpi = 150)
cat("QC violin plot saved!\n")

# Filter low quality spots
spots_before <- ncol(brain)
brain <- subset(brain,
  subset = nFeature_Spatial > 500  &   # remove spots with too few genes
           nCount_Spatial   > 1000 &   # remove spots with too few UMIs
           percent.mt       < 25       # remove spots with high mito content
)
cat("\nSpots before filtering:", spots_before, "\n")
cat("Spots after filtering: ", ncol(brain), "\n")
cat("Spots removed:         ", spots_before - ncol(brain), "\n")


# -----------------------------------------------------------------------------
# SECTION 5 — Normalization
# -----------------------------------------------------------------------------
# Same SCTransform as Day 1
# Corrects for differences in sequencing depth between spots
# Selects 3000 most variable genes

cat("\nRunning SCTransform...\n")
brain <- SCTransform(
  brain,
  assay           = "Spatial",
  vars.to.regress = "percent.mt",
  verbose         = TRUE
)
cat("SCTransform done!\n")
cat("Variable features:", length(VariableFeatures(brain)), "\n")


# -----------------------------------------------------------------------------
# SECTION 6 — Dimensionality Reduction and Clustering
# -----------------------------------------------------------------------------
# Same workflow as Day 1:
# PCA → UMAP → FindNeighbors → FindClusters
#
# But now each "cell" is actually a tissue spot
# Clusters = groups of spots with similar gene expression
# These often correspond to anatomical brain regions

cat("\nRunning PCA...\n")
brain <- RunPCA(brain, verbose = FALSE)

cat("Running UMAP...\n")
brain <- RunUMAP(brain, dims = 1:30, verbose = FALSE)

cat("Finding clusters...\n")
brain <- FindNeighbors(brain, dims = 1:30, verbose = FALSE)
brain <- FindClusters(brain, resolution = 0.5, verbose = FALSE)

cat("Clusters found:", length(levels(brain$seurat_clusters)), "\n")
print(table(brain$seurat_clusters))


# -----------------------------------------------------------------------------
# SECTION 7 — Spatial Visualization
# -----------------------------------------------------------------------------
# THIS is what makes spatial transcriptomics special
# Instead of UMAP, plot spots at their actual X,Y tissue positions
# Color by cluster, gene expression, or any metadata
#
# Reading a spatial plot:
#   Each dot = one Visium spot (~10-20 cells)
#   Position = actual physical location on brain tissue section
#   Color = whatever variable you're plotting
#
# If your annotation is correct:
#   White matter regions → high Plp1/Mbp expression
#   Olfactory bulb → high Omp/S100a5 expression
#   Hypothalamus → specific neuropeptide genes

# Plot 1 — clusters on tissue
p_spatial_clusters <- ggplot(brain@meta.data,
                              aes(x = pxl_col, y = -pxl_row,
                                  color = seurat_clusters)) +
  geom_point(size = 1.2, alpha = 0.8) +
  scale_color_manual(values = rainbow(10)) +
  labs(title    = "Spatial Clusters — Mouse Brain Anterior",
       subtitle = "10x Visium | Each dot = one spot (~10-20 cells)",
       x = "", y = "", color = "Cluster") +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "right")

ggsave("results/day3/spatial_clusters.png",
       plot = p_spatial_clusters, width = 10, height = 8, dpi = 150)
cat("Spatial cluster plot saved!\n")

# Get normalized expression matrix for gene plotting
expr_matrix <- GetAssayData(brain, assay = "SCT", layer = "data")

# Plot known brain marker genes spatially
# These should show anatomically meaningful patterns if data is good
brain_markers <- c(
  "Mbp",    # Myelin Basic Protein — white matter tracts
  "Gfap",   # Glial Fibrillary Acidic Protein — astrocytes
  "Cx3cr1", # Microglia — immune surveillance
  "Snap25", # Synaptosomal protein — neurons
  "Pdgfra"  # Platelet-derived growth factor — OPCs
)

# Keep only genes present in dataset
brain_markers <- brain_markers[brain_markers %in% rownames(brain)]
cat("\nAvailable canonical markers:", brain_markers, "\n")

# Plot each marker spatially
for (gene in brain_markers) {
  expr <- as.numeric(expr_matrix[gene, colnames(brain)])
  brain@meta.data[[paste0("expr_", gene)]] <- expr

  p_gene <- ggplot(brain@meta.data,
                   aes(x = pxl_col, y = -pxl_row,
                       color = .data[[paste0("expr_", gene)]])) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_color_gradient(low = "grey90", high = "darkred") +
    labs(title    = paste("Spatial expression:", gene),
         subtitle = "Mouse Brain Anterior — 10x Visium",
         x = "", y = "", color = "Expression\n(SCT normalized)") +
    theme_void() +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

  ggsave(paste0("results/day3/spatial_", gene, ".png"),
         width = 8, height = 7, dpi = 150)
  cat("Saved: spatial_", gene, ".png\n", sep = "")
}


# -----------------------------------------------------------------------------
# SECTION 8 — Spatially Variable Genes
# -----------------------------------------------------------------------------
# Regular variable genes = differ between cells/spots (ignores location)
# Spatially variable genes = expression is spatially organized
#
# A spatially variable gene forms PATTERNS on the tissue:
#   High in cortex, low in striatum
#   High in white matter tracts, low in grey matter
#   Gradients across brain regions
#
# We measure this using spatial correlation:
#   How much does a gene's expression correlate with X or Y position?
#   High correlation = spatially organized = spatially variable gene

cat("\nCalculating spatially variable genes...\n")

# Get coordinates matrix
coords_mat <- as.matrix(brain@meta.data[, c("pxl_col", "pxl_row")])

# Calculate variance per gene — start with top variable genes
gene_vars     <- apply(expr_matrix, 1, var)
top_var_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]

# For each gene, calculate correlation with spatial position
spatial_scores <- data.frame(
  gene          = top_var_genes,
  cor_x         = NA,
  cor_y         = NA,
  spatial_score = NA
)

for (i in seq_along(top_var_genes)) {
  gene  <- top_var_genes[i]
  expr  <- as.numeric(expr_matrix[gene, ])
  cor_x <- abs(cor(expr, coords_mat[, "pxl_col"], method = "spearman"))
  cor_y <- abs(cor(expr, coords_mat[, "pxl_row"], method = "spearman"))
  spatial_scores$cor_x[i]         <- cor_x
  spatial_scores$cor_y[i]         <- cor_y
  spatial_scores$spatial_score[i] <- max(cor_x, cor_y)
}

# Sort by spatial score
spatial_scores <- spatial_scores[order(-spatial_scores$spatial_score), ]
top_svgs       <- head(spatial_scores$gene, 12)

cat("Top 12 spatially variable genes:\n")
print(spatial_scores[1:12, c("gene", "spatial_score")])

# Plot top 6 SVGs individually
for (gene in top_svgs[1:6]) {
  expr <- as.numeric(expr_matrix[gene, colnames(brain)])
  brain@meta.data[[paste0("svg_", gene)]] <- expr

  p_svg <- ggplot(brain@meta.data,
                  aes(x = pxl_col, y = -pxl_row,
                      color = .data[[paste0("svg_", gene)]])) +
    geom_point(size = 1.5, alpha = 0.8) +
    scale_color_gradient(low = "grey90", high = "darkred") +
    labs(title    = paste("Spatially Variable Gene:", gene),
         subtitle = "Mouse Brain Anterior — 10x Visium",
         x = "", y = "", color = "Expression") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  ggsave(paste0("results/day3/svg_", gene, ".png"),
         width = 8, height = 7, dpi = 150)
  cat("Saved: svg_", gene, ".png\n", sep = "")
}

# Combined panel of top 6 SVGs
svg_plots <- list()
for (gene in top_svgs[1:6]) {
  svg_plots[[gene]] <- ggplot(brain@meta.data,
    aes(x = pxl_col, y = -pxl_row,
        color = .data[[paste0("svg_", gene)]])) +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_gradient(low = "grey90", high = "darkred") +
    labs(title = gene, x = "", y = "", color = "") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

p_svg_combined <- wrap_plots(svg_plots, ncol = 3)
ggsave("results/day3/svg_top6_combined.png",
       plot = p_svg_combined, width = 14, height = 9, dpi = 150)
cat("Combined SVG panel saved!\n")


# -----------------------------------------------------------------------------
# SECTION 9 — Find Spatial Cluster Markers and Annotate
# -----------------------------------------------------------------------------
# Find genes that define each spatial cluster
# Then use these to identify which brain region each cluster represents

cat("\nFinding spatial cluster markers...\n")
Idents(brain) <- "seurat_clusters"

spatial_markers <- FindAllMarkers(
  brain,
  assay           = "SCT",
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25,
  verbose         = FALSE
)

cat("Spatial markers found:", nrow(spatial_markers), "\n")

# Top 2 markers per cluster for annotation
cat("\nTop 2 markers per cluster:\n")
print(spatial_markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC) %>%
  select(cluster, gene, avg_log2FC) %>%
  arrange(cluster), n = 20)

# Save markers
write.csv(spatial_markers,
          "results/day3/spatial_cluster_markers.csv",
          row.names = FALSE)

# Annotate clusters based on marker genes
# Cluster → Brain region mapping:
# 0: Ido1, Cd4       → Immune/Meningeal
# 1: Ighm, Krt80     → Choroid Plexus
# 2: Tnnc1, Hkdc1    → Vascular
# 3: Plp1, Trf       → White Matter (oligodendrocytes)
# 4: Magel2, Arhgap36 → Hypothalamus
# 5: S100a5, Omp     → Olfactory Bulb
# 6: Myoc, Ogn       → Meninges
# 7: Trim54, Ucn     → Brainstem
# 8: Adamts19, Dnah11 → Choroid Plexus 2
# 9: Steap1, Smim22  → Specialized Neurons

brain <- RenameIdents(brain,
  "0" = "Immune_Meningeal",
  "1" = "Choroid_Plexus",
  "2" = "Vascular",
  "3" = "White_Matter",
  "4" = "Hypothalamus",
  "5" = "Olfactory_Bulb",
  "6" = "Meninges",
  "7" = "Brainstem",
  "8" = "Choroid_Plexus_2",
  "9" = "Specialized_Neuron"
)

brain$region_annotated <- Idents(brain)

cat("\nSpot counts per brain region:\n")
print(table(brain$region_annotated))


# -----------------------------------------------------------------------------
# SECTION 10 — Final Spatial Plots with Annotation
# -----------------------------------------------------------------------------

# Annotated regions plot
p_regions <- ggplot(brain@meta.data,
                    aes(x = pxl_col, y = -pxl_row,
                        color = region_annotated)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_brewer(palette = "Set3") +
  labs(title    = "Mouse Brain — Annotated Spatial Regions",
       subtitle = "10x Visium Anterior Section",
       x = "", y = "", color = "Region") +
  theme_void() +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle   = element_text(hjust = 0.5, color = "grey40"),
        legend.position = "right")

ggsave("results/day3/spatial_annotated_regions.png",
       plot = p_regions, width = 10, height = 8, dpi = 150)
cat("Annotated regions plot saved!\n")

# Validation: Plp1 expression should overlap White Matter cluster
expr_plp1 <- as.numeric(expr_matrix["Plp1", colnames(brain)])
brain@meta.data$Plp1_expr <- expr_plp1

p_plp1 <- ggplot(brain@meta.data,
                 aes(x = pxl_col, y = -pxl_row, color = Plp1_expr)) +
  geom_point(size = 1.5, alpha = 0.8) +
  scale_color_gradient(low = "grey90", high = "darkred") +
  labs(title    = "Plp1 — White Matter Validation",
       subtitle = "Should overlap with White_Matter cluster annotation",
       x = "", y = "", color = "Expression") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/day3/spatial_Plp1_validation.png",
       plot = p_plp1, width = 8, height = 7, dpi = 150)
cat("Plp1 validation plot saved!\n")


# -----------------------------------------------------------------------------
# SECTION 11 — Save and Summary
# -----------------------------------------------------------------------------
saveRDS(brain, "data/alzheimers/seurat_day3_brain.rds")
cat("Day 3 object saved!\n")

cat("\n")
cat("=============================================================\n")
cat("DAY 3 COMPLETE — Summary\n")
cat("=============================================================\n")
cat("Dataset:           10x Visium Mouse Brain Anterior\n")
cat("Spots loaded:      2,695\n")
cat("Spots after QC:    2,611\n")
cat("Genes analyzed:    32,285\n")
cat("Spatial clusters:  10\n")
cat("Brain regions:     10 annotated\n")

cat("\nTop spatially variable genes:\n")
print(spatial_scores[1:6, c("gene", "spatial_score")])

cat("\nKey regions mapped:\n")
cat("  White Matter   → Plp1, Trf (myelin)\n")
cat("  Olfactory Bulb → S100a5, Omp\n")
cat("  Hypothalamus   → Magel2, Arhgap36\n")
cat("  Choroid Plexus → Ighm\n")
cat("  Meninges       → Myoc, Ogn\n")

cat("\nOutputs:\n")
cat("  results/day3/QC_spatial_violin.png\n")
cat("  results/day3/spatial_clusters.png\n")
cat("  results/day3/spatial_Mbp/Gfap/Cx3cr1/Snap25/Pdgfra.png\n")
cat("  results/day3/svg_top6_combined.png\n")
cat("  results/day3/spatial_annotated_regions.png\n")
cat("  results/day3/spatial_Plp1_validation.png\n")
cat("  results/day3/spatial_cluster_markers.csv\n")
cat("  data/alzheimers/seurat_day3_brain.rds\n")
cat("=============================================================\n")

sessionInfo()
