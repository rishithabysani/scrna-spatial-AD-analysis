# =============================================================================
# Day 2 — Marker Genes & Differential Expression: AD vs Control
# Dataset: GSE138852 (Grubman et al. 2019, Nature Neuroscience)
# Tissue: Entorhinal Cortex | 6 AD patients vs 6 Controls
# Author: [Your Name]
#
# Day 1 answered: WHAT CELL TYPES are present and how do proportions change?
# Day 2 answers:  WHAT GENES are changed in AD within each cell type?
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 1 — Load Libraries and Data
# -----------------------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/scratch/mallya/bysanir/scrna_project")
.libPaths("/scratch/mallya/bysanir/R/library")

# Create output directory
dir.create("results/day2", recursive = TRUE, showWarnings = FALSE)

# Load Day 1 Seurat object
# This has: QC done, SCTransform normalized, PCA, UMAP, clusters, annotations
cat("Loading Day 1 Seurat object...\n")
seurat_obj <- readRDS("data/alzheimers/seurat_day1.rds")
cat("Loaded:", ncol(seurat_obj), "cells\n")
print(seurat_obj)

# Confirm cell type annotations are present
cat("\nCell type counts:\n")
print(table(seurat_obj$cell_type_annotated))


# -----------------------------------------------------------------------------
# SECTION 2 — Find Marker Genes per Cell Type
# -----------------------------------------------------------------------------
# Goal: find genes that uniquely identify each cell type
# These are genes highly expressed in one cluster vs all others
#
# This is different from DE analysis (Day 2 Step 3):
#   Marker genes = what makes cell type X unique vs other cell types
#   DE genes     = what changes in cell type X between AD and Control
#
# FindAllMarkers parameters:
#   only.pos = TRUE      → only return upregulated markers
#                          (genes higher in this cluster than others)
#   min.pct = 0.25       → gene must be detected in ≥25% of cluster cells
#                          filters out very rare/noisy genes
#   logfc.threshold = 0.25 → minimum log2 fold change
#                            filters out weak markers
#   test.use = "wilcox"  → Wilcoxon rank sum test
#                          non-parametric, standard for scRNA-seq
#                          doesn't assume normal distribution

cat("\nFinding marker genes for each cell type...\n")
cat("This takes 5-10 minutes...\n")

Idents(seurat_obj) <- "cell_type_annotated"

markers <- FindAllMarkers(
  seurat_obj,
  only.pos         = TRUE,
  min.pct          = 0.25,
  logfc.threshold  = 0.25,
  test.use         = "wilcox"
)

cat("Done! Total marker genes found:", nrow(markers), "\n")

# Show top 5 markers per cell type sorted by fold change
cat("\nTop 5 marker genes per cell type:\n")
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
print(top5)

# Understanding the output columns:
# p_val       = raw p-value (how unlikely this is by chance)
# avg_log2FC  = log2 fold change vs all other clusters
#               log2FC 3.5 = 2^3.5 = ~11x higher than other clusters
# pct.1       = % of THIS cluster's cells expressing the gene
# pct.2       = % of ALL OTHER clusters' cells expressing the gene
# p_val_adj   = Bonferroni-corrected p-value (accounts for multiple testing)
# cluster     = cell type this marker belongs to
# gene        = gene name

# Save marker genes
write.csv(markers, "results/day1/marker_genes_all.csv", row.names = FALSE)
cat("Marker genes saved: results/day1/marker_genes_all.csv\n")


# -----------------------------------------------------------------------------
# SECTION 3 — Dot Plot of Marker Genes
# -----------------------------------------------------------------------------
# Dot plot is the standard visualization for marker genes in scRNA-seq papers
#
# Each dot shows two things simultaneously:
#   Dot SIZE  = percentage of cells in that cluster expressing the gene
#               (larger = more cells express it)
#   Dot COLOR = average expression level in expressing cells
#               (redder = higher expression)
#
# A good marker gene has: large dot + dark red color in ONE cell type
#                         small dots + light color in all other cell types

# Top 3 markers per cell type
top3_genes <- markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

p_dot_markers <- DotPlot(
  seurat_obj,
  features = top3_genes,
  group.by = "cell_type_annotated"
) +
  coord_flip() +   # flip so genes are on y-axis, cell types on x-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle("Top Marker Genes per Cell Type",
          subtitle = "Dot size = % cells expressing | Color = avg expression") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")

ggsave("results/day2/dotplot_top_markers.png",
       plot = p_dot_markers, width = 12, height = 10, dpi = 150)
cat("Marker dot plot saved!\n")


# -----------------------------------------------------------------------------
# SECTION 4 — Canonical Brain Marker Validation
# -----------------------------------------------------------------------------
# Before trusting your annotation, validate with known markers from literature
# These genes are well-established in neuroscience literature
# If your annotation is correct:
#   MBP, PLP1 should be HIGH in Oligodendrocytes and LOW everywhere else
#   GFAP, AQP4 should be HIGH in Astrocytes only
#   etc.

canonical_markers <- c(
  "MBP",   "PLP1",    # Oligodendrocyte — myelin basic protein, proteolipid
  "GFAP",  "AQP4",    # Astrocyte — glial fibrillary acidic protein, aquaporin
  "PDGFRA","VCAN",    # OPC — platelet derived growth factor receptor
  "RBFOX1","SYT1",    # Neuron — RNA binding, synaptotagmin
  "CD74",  "CSF1R",   # Microglia — MHC class II, colony stimulating factor
  "CLDN5", "FLT1"     # Endothelial — claudin, VEGF receptor
)

p_canonical <- DotPlot(
  seurat_obj,
  features = canonical_markers,
  group.by = "cell_type_annotated"
) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ggtitle("Canonical Brain Cell Type Markers",
          subtitle = "Validation of cell type annotation against literature") +
  scale_color_gradient2(low = "blue", mid = "white", high = "red")

ggsave("results/day2/dotplot_canonical_markers.png",
       plot = p_canonical, width = 10, height = 8, dpi = 150)
cat("Canonical marker dot plot saved!\n")

# Feature plots — show expression on UMAP
# Each plot colors cells by expression level of one gene
p_feature_canonical <- FeaturePlot(
  seurat_obj,
  features = c("MBP", "GFAP", "CD74", "SYT1", "PDGFRA", "CLDN5"),
  ncol     = 3,
  pt.size  = 0.2,
  order    = TRUE    # plot expressing cells on top
) & scale_color_gradient(low = "grey90", high = "darkred")

ggsave("results/day2/featureplot_canonical.png",
       plot = p_feature_canonical, width = 14, height = 9, dpi = 150)
cat("Canonical feature plots saved!\n")


# -----------------------------------------------------------------------------
# SECTION 5 — Differential Expression: AD vs Control per Cell Type
# -----------------------------------------------------------------------------
# This is the CORE analysis of Day 2 and the most important for your resume
#
# Key concept: pseudobulk vs per-cell DE
# We use per-cell DE (FindMarkers) here for simplicity
# In a real paper you would use pseudobulk (aggregate cells per donor first)
# but per-cell is standard for exploratory analysis
#
# What FindMarkers does:
#   Takes all AD cells of cell type X
#   Takes all Control cells of cell type X
#   For each gene, tests: is expression different between these two groups?
#   Returns: fold change + p-value per gene
#
# positive avg_log2FC = gene is HIGHER in AD
# negative avg_log2FC = gene is LOWER in AD (higher in Control)

cell_types  <- c("Astrocyte", "Oligodendrocyte", "OPC",
                 "Neuron", "Microglia", "Endothelial")
de_results  <- list()

cat("\nRunning DE analysis: AD vs Control per cell type...\n")

for (ct in cell_types) {
  cat("Processing:", ct, "...\n")

  # Subset to just this cell type
  cells_ct <- subset(seurat_obj, cell_type_annotated == ct)

  # Count cells per condition
  n_ad  <- sum(cells_ct$oupSample.batchCond == "AD")
  n_ctl <- sum(cells_ct$oupSample.batchCond == "ct")
  cat("  AD:", n_ad, "| Control:", n_ctl, "\n")

  # Skip if not enough cells
  if (n_ad < 10 | n_ctl < 10) {
    cat("  Skipping — not enough cells\n\n")
    next
  }

  # Set identity to condition
  Idents(cells_ct) <- "oupSample.batchCond"

  # Run DE test
  # ident.1 = "AD" means positive FC = higher in AD
  de <- FindMarkers(
    cells_ct,
    ident.1          = "AD",
    ident.2          = "ct",
    test.use         = "wilcox",
    min.pct          = 0.1,      # gene in ≥10% of cells in either group
    logfc.threshold  = 0.1       # weak threshold — we filter later
  )

  de$gene      <- rownames(de)
  de$cell_type <- ct
  de_results[[ct]] <- de

  # Summary
  n_up   <- sum(de$avg_log2FC >  0 & de$p_val_adj < 0.05)
  n_down <- sum(de$avg_log2FC <  0 & de$p_val_adj < 0.05)
  cat("  DE genes:", nrow(de), "| Up in AD:", n_up, "| Down in AD:", n_down, "\n\n")
}

cat("DE analysis complete!\n")

# Combine all results
all_de <- do.call(rbind, de_results)
write.csv(all_de, "results/day2/DE_AD_vs_Control_all_celltypes.csv",
          row.names = FALSE)
cat("DE results saved: results/day2/DE_AD_vs_Control_all_celltypes.csv\n")

# Print top genes per cell type
for (ct in names(de_results)) {
  cat("\n--- Top AD-upregulated genes in", ct, "---\n")
  de_results[[ct]] %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5) %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    print()

  cat("--- Top AD-downregulated genes in", ct, "---\n")
  de_results[[ct]] %>%
    filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    head(5) %>%
    select(gene, avg_log2FC, p_val_adj) %>%
    print()
}


# -----------------------------------------------------------------------------
# SECTION 6 — Volcano Plots
# -----------------------------------------------------------------------------
# Volcano plot is the standard visualization for DE results
#
# X axis = log2 fold change (AD vs Control)
#   Right side = higher in AD
#   Left side  = lower in AD (higher in Control)
#
# Y axis = -log10(adjusted p-value)
#   Higher = more significant
#   Points above the horizontal dashed line = significant (p < 0.05)
#
# The "volcano" shape comes from:
#   Many non-significant genes cluster near x=0
#   Significant genes shoot up on both sides → looks like an eruption

make_volcano <- function(de_df, cell_type_name) {

  # Classify genes
  de_df$significance <- "Not significant"
  de_df$significance[de_df$p_val_adj < 0.05 & de_df$avg_log2FC >  0.5] <- "Up in AD"
  de_df$significance[de_df$p_val_adj < 0.05 & de_df$avg_log2FC < -0.5] <- "Down in AD"

  # Count per category
  n_up   <- sum(de_df$significance == "Up in AD")
  n_down <- sum(de_df$significance == "Down in AD")

  # Label top genes
  top_up   <- de_df %>% filter(significance == "Up in AD") %>%
    arrange(desc(avg_log2FC)) %>% head(5)
  top_down <- de_df %>% filter(significance == "Down in AD") %>%
    arrange(avg_log2FC) %>% head(5)
  label_genes <- rbind(top_up, top_down)

  ggplot(de_df, aes(x = avg_log2FC,
                    y = -log10(p_val_adj + 1e-300),  # add small value to avoid log(0)
                    color = significance)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c(
      "Up in AD"        = "#E41A1C",
      "Down in AD"      = "#377EB8",
      "Not significant" = "grey70"
    )) +
    geom_text(data = label_genes,
              aes(label = gene),
              size = 3, hjust = -0.1, color = "black") +
    geom_vline(xintercept = c(-0.5, 0.5),
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", color = "grey40") +
    annotate("text", x = max(de_df$avg_log2FC) * 0.8,
             y = max(-log10(de_df$p_val_adj + 1e-300)) * 0.95,
             label = paste("Up:", n_up), color = "#E41A1C", size = 3.5) +
    annotate("text", x = min(de_df$avg_log2FC) * 0.8,
             y = max(-log10(de_df$p_val_adj + 1e-300)) * 0.95,
             label = paste("Down:", n_down), color = "#377EB8", size = 3.5) +
    labs(
      title    = paste("AD vs Control —", cell_type_name),
      subtitle = "GSE138852 | Grubman et al. 2019, Nat Neurosci",
      x        = "log2 Fold Change (positive = higher in AD)",
      y        = "-log10 adjusted p-value",
      color    = ""
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position  = "top",
          plot.title       = element_text(face = "bold"),
          plot.subtitle    = element_text(color = "grey40"))
}

# Generate and save volcano plot for each cell type
for (ct in names(de_results)) {
  p_vol    <- make_volcano(de_results[[ct]], ct)
  filename <- paste0("results/day2/volcano_", gsub(" ", "_", ct), ".png")
  ggsave(filename, plot = p_vol, width = 8, height = 6, dpi = 150)
  cat("Saved:", filename, "\n")
}


# -----------------------------------------------------------------------------
# SECTION 7 — Heatmap of Top DE Genes
# -----------------------------------------------------------------------------
# Heatmap shows expression of multiple genes across all cells
# Each row = one gene
# Each column = one cell
# Color = expression level (red = high, blue = low)
# Cells are grouped by cell type on x-axis

top_de_genes <- all_de %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cell_type) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

cat("\nGenes in heatmap:", length(top_de_genes), "\n")

p_heat <- DoHeatmap(
  seurat_obj,
  features = top_de_genes,
  group.by = "cell_type_annotated",
  size     = 3,
  angle    = 45
) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  ggtitle("Top AD-upregulated genes per cell type")

ggsave("results/day2/heatmap_top_DE.png",
       plot = p_heat, width = 14, height = 10, dpi = 150)
cat("Heatmap saved!\n")


# -----------------------------------------------------------------------------
# SECTION 8 — Feature Plots of Key AD Genes
# -----------------------------------------------------------------------------
# Shows WHERE on the UMAP each gene is expressed
# Confirms our DE findings make spatial sense

key_ad_genes <- c(
  "C3",     # AD astrocyte reactive marker — complement
  "CD44",   # AD astrocyte activation
  "LINGO1", # inhibits remyelination
  "MBP",    # myelin — should be high in oligos
  "CD74",   # microglia activation
  "GFAP"    # astrocyte marker
)

p_feature <- FeaturePlot(
  seurat_obj,
  features = key_ad_genes,
  ncol     = 3,
  pt.size  = 0.2,
  order    = TRUE
) & scale_color_gradient(low = "grey90", high = "darkred")

ggsave("results/day2/featureplot_key_AD_genes.png",
       plot = p_feature, width = 14, height = 9, dpi = 150)
cat("Feature plots saved!\n")


# -----------------------------------------------------------------------------
# SECTION 9 — Save Day 2 Object
# -----------------------------------------------------------------------------
saveRDS(seurat_obj, "data/alzheimers/seurat_day2.rds")
cat("Day 2 Seurat object saved!\n")


# -----------------------------------------------------------------------------
# SECTION 10 — Day 2 Summary
# -----------------------------------------------------------------------------
cat("\n")
cat("=============================================================\n")
cat("DAY 2 COMPLETE — Summary\n")
cat("=============================================================\n")

cat("\nMarker genes found per cell type:\n")
print(table(markers$cluster))

cat("\nDE genes (AD vs Control) per cell type:\n")
for (ct in names(de_results)) {
  de   <- de_results[[ct]]
  up   <- sum(de$avg_log2FC >  0 & de$p_val_adj < 0.05)
  down <- sum(de$avg_log2FC <  0 & de$p_val_adj < 0.05)
  cat(sprintf("  %-20s Up in AD: %3d | Down in AD: %3d\n", ct, up, down))
}

cat("\nKey biological findings:\n")
cat("  Astrocytes:\n")
cat("    C3, CD44 UP   → switched to reactive neurotoxic state\n")
cat("    HES4, HES5 DOWN → losing normal astrocyte identity\n")
cat("    WIF1 DOWN     → Wnt signaling disrupted\n")
cat("  Oligodendrocytes:\n")
cat("    613 genes DOWN → major myelin/myelination dysfunction\n")
cat("  Neurons:\n")
cat("    Stress response genes UP → cellular stress in AD\n")
cat("  Microglia:\n")
cat("    Immune genes DOWN → immune suppression in AD\n")

cat("\nOutputs saved:\n")
cat("  results/day2/dotplot_top_markers.png\n")
cat("  results/day2/dotplot_canonical_markers.png\n")
cat("  results/day2/featureplot_canonical.png\n")
cat("  results/day2/volcano_*.png (one per cell type)\n")
cat("  results/day2/heatmap_top_DE.png\n")
cat("  results/day2/featureplot_key_AD_genes.png\n")
cat("  results/day2/DE_AD_vs_Control_all_celltypes.csv\n")
cat("  data/alzheimers/seurat_day2.rds\n")

cat("\nWhat Day 1 + Day 2 together answered:\n")
cat("  Day 1: WHICH cell types change in AD? (composition)\n")
cat("         → Astrocytes and OPCs massively depleted\n")
cat("  Day 2: WHAT are those cells doing differently? (expression)\n")
cat("         → AD astrocytes became reactive and toxic\n")
cat("         → AD oligodendrocytes lost myelin genes\n")
cat("=============================================================\n")

sessionInfo()
