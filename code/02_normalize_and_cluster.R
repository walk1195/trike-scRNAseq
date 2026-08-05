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
library(EnhancedVolcano)

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
resDir <- 'out/Figures/Clustering/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------

# Sample IDs
sample_ids <- readLines(paste0(projDir,'/inputs/samples.txt'))

# Read in data
s1 <- readRDS(file=glue::glue("{objDir}merged_object.rds"))


# ------------------------------------------------------------------------------
# Log-normalize and scale
# ------------------------------------------------------------------------------

# Log normalization
s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = 10000)
# Variable Features
s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 3000)
# Plot
top10 <- head(VariableFeatures(s1), 10)
plot1 <- VariableFeaturePlot(s1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
ggsave('top_variable_features.png', plot=plot1 + plot2, width = 10, height = 5, dpi=400)

# Scale
s1 <- ScaleData(s1)

### PCA, UMAP, & clustering

# If 50 or fewer cells...
if (ncol(s1) <= 300) {
  s1 <- RunPCA(s1, npcs = 15)
  
  # Elbow plot
  ElbowPlot(s1, ndims = 15) # Elbow plot
  ggsave('elbow_plot.png', width = 6, height = 5, dpi=400)
  
  # Neighbors
  s1 <- FindNeighbors(s1, dims = 1:15)
  # UMAP
  s1 <- RunUMAP(s1, dims = 1:15) # only do this once
  
  # Clustering
  resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  for (res in resolutions){
    # Cluster
    s1 <- FindClusters(s1, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
    # Plot UMAP
    dir.create('umaps')
    DimPlot(s1, reduction = "umap", group.by = paste0('res',res), alpha = 0.7)
    ggsave(paste0('umaps/clusters_', res, '.png'), width = 6, height = 5, dpi=400)
  }
  
  # Clustree to choose resolution
  clusterings <- s1@meta.data
  clustree(clusterings, prefix = "res")
  ggsave(paste('clustree_diagram.png', sep=""), width = 7, height = 9, dpi=400)
  
  # If more than 50 cells...
} else {
  s1 <- RunPCA(s1)
  
  # Elbow plot
  ElbowPlot(s1, ndims = 40) # Elbow plot
  ggsave('elbow_plot.png', width = 6, height = 5, dpi=400)
  
  
  pc_iterations <- c(15,20,25,30)
  for (n_pcs in pc_iterations) {
    # Neighbors
    s1 <- FindNeighbors(s1, dims = 1:n_pcs)
    # UMAP
    s1 <- RunUMAP(s1, dims = 1:n_pcs) # only do this once
    # Plot
    DimPlot(s1, reduction = "umap", group.by = 'nFeature_RNA', alpha = 0.7) + ggtitle(glue::glue('{n_pcs} PCs'))
    ggsave(paste0(n_pcs,'_pcs.png'), width = 6, height = 5, dpi=400)
  }
  
  # Using 25 pcs as baseline
  s1 <- FindNeighbors(s1, dims = 1:25)
  s1 <- RunUMAP(s1, dims = 1:25) # only do this once
  
  # Clustering
  resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)
  for (res in resolutions){
    # Cluster
    s1 <- FindClusters(s1, resolution=res, cluster.name=paste0('res',res)) # Setting smallest resolution
    # Plot UMAP
    dir.create('umaps')
    DimPlot(s1, reduction = "umap", group.by = paste0('res',res), alpha = 0.7)
    ggsave(paste0('umaps/clusters_', res, '.png'), width = 6, height = 5, dpi=400)
  }
  
  # Clustree to choose resolution
  clusterings <- s1@meta.data
  clustree(clusterings, prefix = "res")
  ggsave(paste('clustree_diagram.png', sep=""), width = 7, height = 9, dpi=400)
  
}


# Save progress
message("Clustering finished. Saving progress...")
saveRDS(s1, paste0(objDir, sample_id, '.rds'))
message("Clustered object saved.")


########################################################################
#           QC Plots -- Post Clustering
########################################################################
message("Generating QC plots after clustering...")

setwd(resDir)

FeaturePlot(s1, features='nFeature_RNA')
ggsave(paste('umap_nFeature.png', sep=""), width = 6, height = 5, dpi=400)

FeaturePlot(s1, features='nCount_RNA')
ggsave(paste('umap_nCount.png', sep=""), width = 6, height = 5, dpi=400)

FeaturePlot(s1, features='percent.mt')
ggsave(paste('umap_mt.png', sep=""), width = 6, height = 5, dpi=400)



########################################################################
#           Marker Gene Plots -- Post Clustering
########################################################################
dir.create(paste0(resDir, '/Marker_genes'))
setwd(paste0(resDir, '/Marker_genes'))

# -----------------------------------------------------------
plot_features_safe <- function(obj, features, filename) {
  # Keep only genes present in the object
  valid_features <- features[features %in% rownames(obj)]
  # Warn if some are missing
  missing <- setdiff(features, valid_features)
  if (length(missing) > 0) {
    message("Missing genes: ", paste(missing, collapse = ", "))
  }
  # Skip if no valid genes
  if (length(valid_features) == 0) {
    message("No valid genes found. Skipping: ", filename)
    return(NULL)
  }
  # -----------------------------------------------------------
  
  # Plot and save
  p <- FeaturePlot(obj, features = valid_features)
  ggsave(filename, plot = p, width = 6, height = 5, dpi = 400)
}

bcell_markers <- c('CD79A', 'CD79B', 'PAX5', 'CD19')
plot_features_safe(s1, bcell_markers, "umap_bcell_markers.png")

tcell_markers <- c('CD3E', 'CD3D', 'CD4', 'CD8A')
plot_features_safe(s1, tcell_markers, "umap_tcell_markers.png")

myeloid <- c('CD163', 'CD68', 'S100A8', 'CSF1R')
plot_features_safe(s1, myeloid, "umap_myeloid_markers.png")

prolif <- c('MKI67', 'TOP2A', 'CENPF','BUB1')
plot_features_safe(s1, prolif, "umap_prolif_markers.png")






