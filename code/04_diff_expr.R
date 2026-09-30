#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Differential Gene Expression Analysis (Pre-TriKE to Post-TriKE samples)
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: September 20, 2026
#############################################################################################
#############################################################################################

# -------------------------------------------------------------------------------------------
# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(DESeq2)

# Functions
source("code/utils/functions.R")
source("code/utils/colormaps.R")

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Dirs
projDir <- getwd()
dataDir <- 'inputs/data/' # symlinked from data_delivery
objDir <- 'out/Objects/'
resDir <- 'out/Figures/DGEA/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))


#############################################################################################

# TriKE therapy has 3 components:
      # (1) CD16: targets NK/T cells
      # (2) IL15 (growth factor): stimulant for any cells expressing IL15R (B cells, T cells)
      # (3) Checkpoint molecule: tumor expression
# Delivered sub-q, distributes to pharmaco state, if it sees a tumor it pulls T cells/NK cells to the tumor

# We have 2 pretx samples, 2 posttx
# We are interested in what is happening in the T cell compartment postx (compared to pretx). Specifically:
      # Changes in GEX patterns reflecting activation
      # Ramping up adhesion molecule genes, proliferation genes
      # Chemokines & receptors, interleukin receptors to prime for IFN-g response
      # Increased IL2R, CCL4, CCL5
      # Pathways associated with tissue residence, migration, adhesion, etc.

# Some options for analysis:
      # Identify genes that are differentially expressed in the tcells between timepoints, for each dog individually
        # Each dog individually 
        # Pseudobulked together, such that each dog is 1 biological replicate -- this is limited since we only have 2 of each
      # GSEA for specific activation or adhesion pathways
      # AUCell?

#############################################################################################

# ------------------------------------------------------------------------------
# Start very broad - all cells in the sample
# ------------------------------------------------------------------------------

# Set Idents
Idents(s1) <- s1$orig.ident

# Join layers
s1 <- JoinLayers(s1)

# Get markers
markers_trike01 <- FindMarkers(object = s1, ident.1 = 'trike_01_2w', ident.2 = 'trike_01_pretx')
markers_trike02 <- FindMarkers(object = s1, ident.1 = 'trike_02_2w', ident.2 = 'trike_02_pretx')

# Plot
prettierVolc(markers_trike01, title = "TriKE 01: All cells", left_lab = "Upregulated PreTx", right_lab = 'Upregulated PostTx')
ggsave(file=paste0(resDir,'volcano_trike01_all_cells.png'), dpi=400, height=6, width=7)
prettierVolc(markers_trike02, title = "TriKE 02: All cells", left_lab = "Upregulated PreTx", right_lab = 'Upregulated PostTx')
ggsave(file=paste0(resDir,'volcano_trike02_all_cells.png'), dpi=400, height=6, width=7)


# ------------------------------------------------------------------------------
# Now let's break it down by cell types - using Ammons labels
# ------------------------------------------------------------------------------

# Create col of cell type + sample condition
s1$celltype.cond <- paste(s1$ammons.l1, s1$orig.ident, sep = "_")
unique(s1$celltype.cond)

# Set idents
Idents(s1) <- s1$celltype.cond
 
#### Run t-tests

# TriKE01
for (ct in unique(s1$ammons.l1)) {
  markers <- FindMarkers(s1, ident.1 = paste0(ct, "_trike_01_2w"), ident.2 = paste0(ct, "_trike_01_pretx"), min.pct = 0.1)
  
  # Plot
  p <- prettierVolc(markers, title = paste0("TriKE 01: ", ct), left_lab = "Upregulated PreTx", right_lab = "Upregulated PostTx")
  
  # Save
  ct_file <- gsub("[^A-Za-z0-9_-]", "_", ct)
  ggsave(paste0(resDir, "t-tests/volcano_trike01_", ct_file, ".png"), p, dpi = 400, height = 6, width = 7)
  
  # Save markers too
  write.table(markers, paste0(resDir, "t-tests/degs/markers_trike01_", ct_file, ".tsv"), sep = "\t", quote = FALSE)
}

# TriKE02
for (ct in unique(s1$ammons.l1)) {
  # Get markers
  markers <- FindMarkers(s1, ident.1 = paste0(ct, "_trike_02_2w"), ident.2 = paste0(ct, "_trike_02_pretx"), min.pct = 0.1)
  
  # Plot
  p <- prettierVolc(markers, title = paste0("TriKE 02: ", ct), left_lab = "Upregulated PreTx", right_lab = "Upregulated PostTx")

  # Save
  ct_file <- gsub("[^A-Za-z0-9_-]", "_", ct)
  ggsave(paste0(resDir, "t-tests/volcano_trike02_", ct_file, ".png"), p, dpi = 400, height = 6, width = 7)
  
  # Save markers list too
  write.table(markers, paste0(resDir, "t-tests/markers_trike02_", ct_file, ".tsv"), sep = "\t", quote = FALSE)
}


# ------------------------------------------------------------------------------
# Looks like there are some genes getting in that appear to be FPs - tweaking params to see if we can clean up results
# ------------------------------------------------------------------------------

# CD8/NK
markers_trike01_cd8 <- FindMarkers(object = s1, ident.1 = 'CD8/NK cell_trike_01_2w', ident.2 = 'CD8/NK cell_trike_01_pretx', logfc.threshold=0.6, min.pct=.25, min.diff.pct=0.1)
prettierVolc(markers_trike01_cd8, title = "TriKE 01: CD8/NK cells (min.pct=.25, min.diff.pct=0.1)", left_lab = "Upregulated PreTx", right_lab = 'Upregulated PostTx')
ggsave(file=paste0(resDir,'volcano_2_trike01_CD8_NK.png'), dpi=400, height=6, width=7)

# B cells (control?)
markers_trike01_bcell <- FindMarkers(object = s1, ident.1 = 'B cell_trike_01_2w', ident.2 = 'B cell_trike_01_pretx', logfc.threshold=0.6, min.pct=.25, min.diff.pct=0.1)
prettierVolc(markers_trike01_bcell, title = "TriKE 01: B cells (min.pct=.25, min.diff.pct=0.1)", left_lab = "Upregulated PreTx", right_lab = 'Upregulated PostTx')
ggsave(file=paste0(resDir,'volcano_2_trike01_B_cell.png'), dpi=400, height=6, width=7)

# CD4
markers_trike01_cd4 <- FindMarkers(object = s1, ident.1 = 'CD4 T cell_trike_01_2w', ident.2 = 'CD4 T cell_trike_01_pretx',  logfc.threshold=0.6, min.pct=.25, min.diff.pct=0.1)
prettierVolc(markers_trike01_cd4, title = "TriKE 01: CD4 cells (min.pct=.25, min.diff.pct=0.1)", left_lab = "Upregulated PreTx", right_lab = 'Upregulated PostTx')
ggsave(file=paste0(resDir,'volcano_2_trike01_CD4.png'), dpi=400, height=6, width=7)


# ------------------------------------------------------------------------------
# We are dealing with very high p-values and likely many false positives still
# 2 possible solutions:
    # (1) Using pseudobulked expression of cell types + tx group with DESeq2
    # (2) Subsetting cell type groups to 500 randomized cells to decrease statistical power
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# (1) Run DESeq2 with pseudobulked expression

# The Seurat wrapper for DESeq2 enforces a minimum of 3 replicates per condition, so we need to build the de object from scratch
# ------------------------------------------------------------------------------

# Pseudobulk expression
pseudo.s1 <- AggregateExpression(s1, assays = "RNA", return.seurat = T, group.by = c("orig.ident", "ammons.l1"))
head(Cells(pseudo.s1))

# Raw pseudobulk counts
counts <- GetAssayData(pseudo.s1, assay = "RNA", layer = "counts")

# Add metadata
pseudo.s1$condition <- ifelse(grepl("-pretx", pseudo.s1$orig.ident), "pretx", "2w")
pseudo.s1$dog <- sub("_(pretx|2w)$", "", pseudo.s1$orig.ident)

# Fix
pseudo.s1$dog <- sub("-(2w|pretx)_.*$", "", pseudo.s1$orig.ident)
table(pseudo.s1$dog, pseudo.s1$condition)

# Run paired DESeq2 for every cell type
      # We are asking whether, after accounting for differences in baseline expr levels across dogs (replicates), do we see
      # evidence of of gex changes Pre-TriKE to Post-TriKE?
ct <- "CD8/NK cell"

de_results <- lapply(unique(pseudo.s1$ammons.l1), function(ct) {
  # Grab samples per cell type 
  samples <- which(pseudo.s1$ammons.l1 == ct)
  # Grab pseudobulk counts for given cell type
  ct_counts <- counts[, samples, drop = FALSE]
  
  # Build metadata
  coldata <- data.frame(dog = factor(pseudo.s1$dog[samples]), condition = factor(pseudo.s1$condition[samples], levels = c("pretx", "2w")))
  rownames(coldata) <- colnames(ct_counts)
  
  # Remove genes with 0 counts
  ct_counts <- ct_counts[rowSums(ct_counts) > 0, , drop = FALSE]
  
  # Build DESeq2 obj -- this normalizes counts and fits a neg binomial model for every gene 
  dds <- DESeqDataSetFromMatrix(
    countData = round(ct_counts),
    colData = coldata,
    design = ~ dog + condition)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results + some reformatting
  res <- as.data.frame(results(dds, contrast = c("condition", "2w", "pretx")))
  res$gene <- rownames(res)
  res$cell_type <- ct
  
  
  # Plot top 10 genes as a heatmap
  # Top 10 genes by adjusted p-value
  top_genes <- res |>
    dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) |>
    dplyr::arrange(padj) |>
    dplyr::slice_head(n = 10) |>
    dplyr::pull(gene)
  
  # For grabbing more top genes to plot
  # top_genes <- res |>
  #   dplyr::filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) |>
  #   dplyr::arrange(padj) |>
  #   dplyr::slice_head(n = 20) |>
  #   dplyr::pull(gene)
  
  if (length(top_genes) > 3) {
    
    # Normalized counts from the same full-gene DESeq2 object
    mat <- counts(dds, normalized = TRUE)[top_genes, , drop = FALSE]
    mat <- t(scale(t(log2(mat + 1))))
    
    # Fixing sample IDs and reordering for viz
    colnames(mat) <- sub("_.*$", "", colnames(mat))
    colnames(mat) <- gsub("-", "_", colnames(mat))
    
    mat <- mat[, c("trike_01_pretx", "trike_01_2w", "trike_02_pretx", "trike_02_2w")]
    
    print('Plotting heatmap.....')
    # Unclustered samples
    pheatmap::pheatmap(
      mat,
      #color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      fontsize_row = 8,
      main = ct,
      filename = paste0(resDir, "heatmaps_sig_genes/heatmap_", gsub("[^A-Za-z0-9_-]", "_", ct), ".png"),
      width = 6,
      height = 6)
  
    # Clustering by sample
    pheatmap::pheatmap(
      mat,
      #color = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100),
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      fontsize_row = 8,
      main = ct,
      filename = paste0(resDir, "heatmaps_sig_genes/heatmap_clustered_", gsub("[^A-Za-z0-9_-]", "_", ct), ".png"),
      width = 6,
      height = 6)
  }
  
  res
}) |> dplyr::bind_rows() # Combining all cell-type results into a single df

# Save res
write.table(de_results, paste0(resDir,"deseq2_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Plotting results - looping over each cell type
for (group in unique(de_results$cell_type)) {
  # Extract deseq2 markers
  bulk_de_results <- de_results |>
    dplyr::filter(cell_type == group) |>
    dplyr::filter(!is.na(padj)) |>
    dplyr::arrange(padj)
  
  # For saving filename
  cell_type_clean <- gsub("[ /]", "_", group)
  
  # Plot
  prettierVolc(bulk_de_results, title = group, left_lab = "Upregulated Pre-TriKE", right_lab = "Upregulated Post-TriKE")
  ggsave(paste0(resDir, 'volcano_bulk_', cell_type_clean, '.png'), dpi=400, height=6, width=7)
}

bulk_de_results <- de_results |>
  dplyr::filter(cell_type == group) |>
  dplyr::filter(!is.na(padj)) |>
  dplyr::arrange(padj)


# Print session info
sessionInfo()
