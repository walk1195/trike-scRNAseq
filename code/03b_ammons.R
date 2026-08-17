#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Label transfer with Ammons et al. canine PBMC cell types
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: August 17, 2026
#############################################################################################
#############################################################################################

# -------------------------------------------------------------------------------------------
# Ammons et al : scRNAseq atlas of circulating canine leukocytes
# https://pubmed.ncbi.nlm.nih.gov/37275879/

# Running automated annotation using this canine PBMC dataset as a reference for additional cross validation of cell types

# -------------------------------------------------------------------------------------------

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse)
library(Seurat)
library(RColorBrewer)

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
resDir <- 'out/Figures/Annotation/'
refDir <- '/scratch.global/walk1195/trike/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))

# Ammons reference
ref <- readRDS(glue::glue("{refDir}GSE225599_final_dataSet_H.rds"))

# ------------------------------------------------------------------------------
# Get ref variable features
# ------------------------------------------------------------------------------

# Reformatting the reference a bit
ref@assays$SCT <- NULL
ref@assays$integrated <- NULL

# Get variable features
ref <- FindVariableFeatures(ref, selection.method = "vst", nfeatures = 3000)

# ------------------------------------------------------------------------------
# Run label transfer
# ------------------------------------------------------------------------------

# Find anchors bw query and reference
anchors <- FindTransferAnchors(reference = ref, query = s1, dims = 1:30)

# Transfer labels using anchors
predictions.l1 <- TransferData(anchorset = anchors, refdata = ref$celltype.l1, dims = 1:30)
predictions.l2 <- TransferData(anchorset = anchors, refdata = ref$celltype.l2, dims = 1:30) 

# Add preds back to seurat obj

# L1
s1 <- AddMetaData(s1, metadata = predictions.l1$predicted.id, col.name='ammons.l1')
s1 <- AddMetaData(s1, metadata = predictions.l1$prediction.score.max, col.name='ammons.l1.score')

# L2
s1 <- AddMetaData(s1, metadata = predictions.l2$predicted.id, col.name='ammons.l2')
s1 <- AddMetaData(s1, metadata = predictions.l2$prediction.score.max, col.name='ammons.l2.score')


# Save results
write.table(predictions, file=glue::glue("{resDir}ammons_predictions_l1.tsv"), sep='\t', quote=F, row.names=T, col.names=T)
write.table(predictions.l2, file=glue::glue("{resDir}ammons_predictions_l2.tsv"), sep='\t', quote=F, row.names=T, col.names=T)

# ------------------------------------------------------------------------------
# Plot results
# ------------------------------------------------------------------------------
# L1 preds
prettierDimPlot(s1, group.by = 'ammons.l1', reduction='umap.harmony', cols=ammonsCols.l1, save_fig = T, resDir=resDir, box_labels = T,
                file_name = 'umap_ammons_l1.png')
prettierDimPlot(s1, group.by = 'ammons.l1', reduction='umap.harmony', cols=ammonsCols.l1, save_fig = T, resDir=resDir, box_labels = F, width=9,
                file_name = 'umap_ammons_l1_with_legend.png')


# L2 preds
prettierDimPlot(s1, group.by = 'ammons.l2', reduction='umap.harmony', cols=ammonsCols.l2, save_fig = T, resDir=resDir, box_labels=T,
                file_name = 'umap_ammons_l2.png')
prettierDimPlot(s1, group.by = 'ammons.l2', reduction='umap.harmony', cols=ammonsCols.l2, save_fig = T, resDir=resDir, box_labels = F, width=10,
                file_name = 'umap_ammons_l2_with_legend.png')

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
saveRDS(s1, file=glue::glue("{objDir}filtered_object.rds"))

##### Session info
sessionInfo()


