#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Normalization & Clustering
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: August 5, 2026
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse)
library(Seurat)
library(patchwork)
library(scales)
library(SingleR)
library(scran)

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
resDir <- 'out/Figures/SingleR/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------

s1 <- readRDS(file=glue::glue("{objDir}merged_object.rds"))

# ------------------------------------------------------------------------------
# Run SingleR
# ------------------------------------------------------------------------------
# Import reference datasets
ref1 <- celldex::HumanPrimaryCellAtlasData() # Human primary scRNAseq atlas 
ref2 <- celldex::BlueprintEncodeData() # Human bulk RNAseq
ref3 <- celldex::ImmGenData() # Mouse immune bulk RNAseq

# Join layers
s1 <- JoinLayers(s1)

# Convert seurat object to sce object
integrated_sce <- as.SingleCellExperiment(s1)

# Object has already been log-normalized

## Convert mouse gene names to uppercase
rownames(ref3) <- toupper(rownames(ref3))

# Subset to common genes
common1 <- intersect(rownames(integrated_sce),rownames(ref1))
common2 <- intersect(rownames(integrated_sce),rownames(ref2)) 
common3 <- intersect(rownames(integrated_sce),rownames(ref3))

# Check common genes
length(common1) 
length(common2)
length(common3)

# Subset
integrated_sce1 <- integrated_sce[common1,]
integrated_sce2 <- integrated_sce[common2,]
integrated_sce3 <- integrated_sce[common3,]


# Run SingleR with different reference datasets -- outputs dataframes
hpca_results_main <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.main, assay.type.ref = "logcounts")
hpca_results_fine <- SingleR(test=integrated_sce1, ref=ref1, labels=ref1$label.fine, assay.type.ref = "logcounts")
blue_results <- SingleR(test=integrated_sce2, ref=ref2, labels=ref2$label.main, assay.type.ref = "logcounts")
immgen_results <- SingleR(test=integrated_sce3, ref=ref3, labels=ref3$label.main, assay.type.ref = "logcounts")

# Save results as csv
write.table(hpca_results_main, paste0(resDir,"SingleR_Results_HPCA_Main.tsv"), sep = "\t", quote = FALSE)
write.table(hpca_results_fine, paste0(resDir,"SingleR_Results_HPCA_Fine.tsv"), sep = "\t", quote = FALSE)
write.table(blue_results, paste0(resDir,"SingleR_Results_BlueEncode.tsv"), sep = "\t", quote = FALSE)
write.table(immgen_results, paste0(resDir,"SingleR_Results_ImmGen.tsv"), sep = "\t", quote = FALSE)


### Plotting Results ###
# Plot score heatmaps
png(paste0(resDir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Main.png"),height=2600,width=4000,res=400)
plotScoreHeatmap(hpca_results_main)
dev.off()

png(paste0(resDir,"SingleR_Heatmap_CellAssignmentScores_HPCA_Fine.png"),height=2600,width=5000, res=400)
plotScoreHeatmap(hpca_results_fine)
dev.off()

png(paste0(resDir,"SingleR_Heatmap_CellAssignmentScores_BlueEncode.png"),height=2600,width=4000, res=400)
plotScoreHeatmap(blue_results)
dev.off()

png(paste0(resDir,"SingleR_Heatmap_CellAssignmentScores_ImmGen.png"),height=2600,width=4000, res=400)
plotScoreHeatmap(immgen_results)
dev.off()


### Merge results with object ###
s1 <- AddMetaData(s1, hpca_results_main$labels, col.name = 'hpca_main_raw')
s1 <- AddMetaData(s1, hpca_results_main$pruned.labels, col.name = 'hpca_main_pruned')

s1 <- AddMetaData(s1, hpca_results_fine$labels, col.name = 'hpca_fine_raw')
s1 <- AddMetaData(s1, hpca_results_fine$pruned.labels, col.name = 'hpca_fine_pruned')

s1 <- AddMetaData(s1, blue_results$labels, col.name = 'blue_raw')
s1 <- AddMetaData(s1, blue_results$pruned.labels, col.name = 'blue_pruned')

s1 <- AddMetaData(s1, immgen_results$labels, col.name = 'immgen_raw')
s1 <- AddMetaData(s1, immgen_results$pruned.labels, col.name = 'immgen_pruned')

# ------------------------------------------------------------------------------
# Plot results onto UMAP
# ------------------------------------------------------------------------------
# Plot - HPCA main
Idents(s1) <- s1@meta.data$hpca_main_pruned
DimPlot(s1, label = F , repel = T, label.size = 3) + ggtitle('HPCA Main Annotations')
ggsave(paste0('hpca_main_umap.png'), height=6, width=9, dpi=400) 

# Plot - HPCA fine
Idents(s1) <- s1@meta.data$hpca_fine_pruned
DimPlot(s1, label = F , repel = T, label.size = 3) + ggtitle('HPCA Fine Annotations')
ggsave(paste0('hpca_fine_umap.png'), height=6, width=20, dpi=400) 

# Plot - BlueEncode
Idents(s1) <- "blue_pruned"
DimPlot(s1, reduction='umap') + ggtitle('BlueEncode Annotations')
ggsave(paste0('blueencode_umap.png'), height=7, width=8, dpi=400)

# Plot - ImmGen
Idents(s1) <- "immgen_pruned"
DimPlot(s1, reduction='umap') + ggtitle('ImmGen Annotations')
ggsave(paste0('immgen_umap.png'), height=5, width=8, dpi=400)

# Save progress
message("SingleR annotation finished. Saving progress...")
saveRDS(s1, paste0(objDir, sample_id, '.rds'))
message("Annotated object saved.")


