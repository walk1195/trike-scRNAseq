#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Normalization & Clustering
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: August 5, 2026
#############################################################################################
#############################################################################################

# -------------------------------------------------------------------------------------------
# TODO:
# (1) Inspect cluster markers of suspicious cell clusters at various resolutions

# (2) Inspect result of changing n_pc parameter
# -------------------------------------------------------------------------------------------

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
# Library-size normalization
# ------------------------------------------------------------------------------

########################### Some notes about normalization ##################################

# https://www.reddit.com/r/bioinformatics/comments/c3jd6s/seurat_scrnaseq_normalization/

# Library-size norm. is by far the most common and simple method of normalizing scRNAseq
# Lib size = total sum of counts for all genes in a given cell
# Normalized counts for a given cell = feature_counts / lib.size * scale.factor
  ### where scale.factor is either
      ### 10e6 (cpm, recommended for bulk RNAseq), or
      ### 10k (seurat recommended for scRNAseq)

# CPM can greatly inflate the difference between a 0 and a 1 count in the log-space
# Say  you had 4 cells, each with 10k total counts, and for a particular gene, they have 0, 1, 2, and 3 counts respectively

  ## Using log(1+counts/million) this gives:
  
     # [0, 1, 2, 3] -> [0, 4.6, 5.3, 5.7]

  ## Using log(1+counts/10k) this gives:
  
     # [0, 1, 2, 3] -> [0, .7, 1.1, 1.4]

# Using cpm introduces a large difference between 0 to 1, and a much smaller difference between 1 count a 2 counts

# With cp10k, the 0-1 difference is only 1.7x the 1-2 difference.

#### In general, it's probably best to use a scale factor that is similar in order of magnitude to the # of counts you're getting.
     # For bulk or even full-length transcript single-cell data, you get millions of counts.
     # With UMI single-cell data you get 1000s to 10,000s of counts.

### Using the median raw count depth of the dataset

# Log normalization -- this is performed on cell-by-cell basis, so having independent counts layers doesn't matter
s1 <- NormalizeData(s1, normalization.method = "LogNormalize", scale.factor = median(s1$nCount_RNA)) # Scale by the median total UMI per cell = 8691 (first pass), 8366 (second pass)

# Join layers
#s1 <- JoinLayer(s1)
#s1[["RNA"]] <- split(s1[["RNA"]], f = s1$orig.ident)

############################################################################################
# Variable Features

# https://github.com/satijalab/seurat/issues/9808
# FindVariableFeatures() will find variable genes for each layer
# then find a union set of variable genes based on the ranking and the number of appearance (as an HVG) of them in each layer/sample.
# If some of the genes are specific to one layer, it will not be counted as a HVG in the final list.

s1 <- FindVariableFeatures(s1, selection.method = "vst", nfeatures = 3000)

s1 <- s1 %>%
  ScaleData() %>%
  RunPCA()

# Elbow plot
p <- ElbowPlot(s1, ndims = 50)
ggsave(file=glue::glue("{resDir}elbow_plot.png"), plot=p, dpi=400, height=5, width=6)

# Neighbors & UMAP
s1 <- s1 %>%
  FindNeighbors(dims = 1:30) %>%
  RunUMAP(dims = 1:30)


####### Integration ########
# Run integration
s1 <- IntegrateLayers(object = s1, method = HarmonyIntegration, orig.reduction = "pca",
                       new.reduction = "harmony", verbose = TRUE)

# Neighbors & UMAP
s1 <- s1 %>%
  FindNeighbors(dims = 1:30, reduction='harmony') %>%
  RunUMAP(dims = 1:30, reduction='harmony', reduction.name='umap.harmony')


# ------------------------------------------------------------------------------
# Plotting cluster results
# ------------------------------------------------------------------------------
resDir <- glue::glue("{resDir}thirdpass/")
dir.create(resDir)

############## UMAPs ###############

### By sample ID
prettierDimPlot(s1, group.by='orig.ident', cols=sampleCols, reduction='umap.harmony', save_fig=T, title=NULL, file_name='umap_by_sample_integrated.png', resDir=resDir)

prettierDimPlot(s1, group.by='orig.ident', cols=sampleCols, save_fig=F, title=NULL, file_name='umap_by_sample_unintegrated.png', resDir=resDir)


### By cell types
# Blue encode
prettierDimPlot(s1, group.by='blue_pruned', cols=encodeCols, reduction='umap.harmony', box_labels = T,
                save_fig=T, title='Blue Encode Annotations', file_name='umap_blue_encode_integrated.png', resDir=resDir)
# Blue encode (unintegrated)
prettierDimPlot(s1, group.by='blue_pruned', cols=encodeCols, box_labels = T,
                save_fig=T, title='Blue Encode Annotations', file_name='umap_blue_encode.png', resDir=resDir)

# Immgen
prettierDimPlot(s1, group.by='immgen_pruned', cols=immgenCols, reduction='umap.harmony', box_labels=T,
                save_fig=T, title='Immgen Annotations', file_name='umap_immgen_integrated.png', resDir=resDir)
# Immgen (unintegrated)
prettierDimPlot(s1, group.by='immgen_pruned', cols=immgenCols, box_labels=T,
                save_fig=T, title='Immgen Annotations', file_name='umap_immgen.png', resDir=resDir)

# HPCA main
prettierDimPlot(s1, group.by='hpca_main_pruned', cols=hpca_mainCols, reduction='umap.harmony', save_fig=T, box_labels=T,
                title='HPCA main Annotations', file_name='umap_hpca_main_integrated.png', resDir=resDir)
# HPCA main (unintegrated)
prettierDimPlot(s1, group.by='hpca_main_pruned', cols=hpca_mainCols, save_fig=T, box_labels=T,
                title='HPCA main Annotations', file_name='umap_hpca_main.png', resDir=resDir)

### By doublet calls
prettierDimPlot(s1, group.by='scDblFinder.class', cols=doubletCols, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='doublets_integrated.png')
prettierDimPlot(s1, group.by='scDblFinder.class', cols=doubletCols, save_fig=T, resDir=resDir, file_name='doublets.png')


############## FeaturePlots ###############
### QC metrics (integrated)
prettierFeatPlot(s1, features='nFeature_RNA', reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='nFeature_integrated.png')
prettierFeatPlot(s1, features='nCount_RNA', reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='nCount_integrated.png')
prettierFeatPlot(s1, features='percent.mt', reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='MT_integrated.png')

### QC metrics (unintegrated)
prettierFeatPlot(s1, features='nFeature_RNA', save_fig=T, resDir=resDir, file_name='nFeature.png')
prettierFeatPlot(s1, features='nCount_RNA', save_fig=T, resDir=resDir, file_name='nCount.png')
prettierFeatPlot(s1, features='percent.mt', save_fig=T, resDir=resDir, file_name='MT.png')

### Cell type markers
resDir <- glue::glue("{resDir}marker_genes/")
dir.create(resDir)

# Markers of interest
activation <- c('CD69', 'IL2RA', 'TNFRSF9')
effector_cytokines <- c('IFNG', 'TNF', 'IL2')
cytotoxicity <- c('GZMB', 'GZMA', 'NKG7', 'PRF1', 'CTSW', 'CCL4', 'CCL5')
effector_mem <- c('CD44', 'CD27', 'IL7R', 'CD28')
prolif <- c('MKI67', 'TOP2A', 'STMN1')
signaling <- c('CD3D', 'CD3E')
# Plot
prettierFeatPlot(s1, features=activation, ncol=2, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='activation_genes.png')
prettierFeatPlot(s1, features=effector_cytokines, ncol=2, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='effector_cyto_genes.png')
prettierFeatPlot(s1, features=cytotoxicity, ncol=3, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='cytotoxicity_genes.png', height=10, width=12)
prettierFeatPlot(s1, features=effector_mem, ncol=2, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='effector_mem_genes.png')
prettierFeatPlot(s1, features=prolif, ncol=2, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='prolif_genes.png')
prettierFeatPlot(s1, features=signaling, ncol=2, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='tcell_genes.png', height=4, width=8)

# ------------------------------------------------------------------------------
# Inspecting cells < 1000 feats
# ------------------------------------------------------------------------------
s1$low_feat <- with(s1@meta.data, ifelse(nFeature_RNA < 1000, TRUE, FALSE))
prettierDimPlot(s1, group.by='low_feat', cols=doubletCols, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='umap_low_feature_cells.png')

prettierDimPlot(s1, group.by='orig.ident', split.by=T, cols=sampleCols, reduction='umap.harmony', save_fig=T, resDir=resDir, file_name='umap_split_by_sample.png', width=15)


# ------------------------------------------------------------------------------
# Clustering at different resolutions
# ------------------------------------------------------------------------------
resDir <- glue::glue("{resDir}umaps/")
dir.create(resDir)

# Grab custom cols to replace seurat default
random_cols <- generate_color_vector(20, seed=10)

resolutions = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
for (res in resolutions){
  # Cluster
  s1 <- FindClusters(s1, resolution=res, cluster.name=paste0('res',res))
  # Plot UMAP
  prettierDimPlot(s1, cols=random_cols, group.by=paste0('res',res), reduction = "umap.harmony", save_fig=T, resDir=resDir, file_name=paste0('res', res, '.png'))
}

# Clustree to choose resolution
clusterings <- s1@meta.data
clustree(clusterings, prefix = "res")
ggsave(paste0(resDir, 'clustree_diagram.png'), width = 7, height = 9, dpi=400)


# ------------------------------------------------------------------------------
# TODO: Calculating cluster markers for a few resolutions 
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# TODO: Testing out a few different pc values
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# TODO: inspect cpms of single genes
# ------------------------------------------------------------------------------
#### Inspecting single gene cpms ----- TODO: we'll come back to this later to create a function that can grab cpms and plot them ###
s1.joined <- JoinLayers(s1)
counts <- GetAssayData(s1.joined, assay = "RNA", layer = "counts")

# expr_df <- FetchData(s1, vars = c(gene, "orig.ident"))
# 
# ggplot(expr_df, aes(x = orig.ident, y = .data[[gene]])) +
#   geom_boxplot(outlier.size = 0.3) +
#   labs(x = NULL, y = paste0(gene, " expression")) +
#   theme_classic()

gene <- 'KLF2'
cpm <- counts[gene, ] / Matrix::colSums(counts) * 1e6

cpm_df <- data.frame(sample = s1$orig.ident, cpm = as.numeric(cpm)) |>
  dplyr::filter(cpm > 0)


ggplot(cpm_df, aes(x = sample, y = cpm)) +
  geom_boxplot(outlier.size = 0.3) +
  labs(x = NULL, y = paste0(gene, " CPM")) +
  theme_classic()




##### Save
# firstpass object
saveRDS(s1, file=glue::glue("{objDir}merged_object.rds"))

# thirdpass object
saveRDS(s1, file=glue::glue("{objDir}filtered_object.rds"))


##### Session info
sessionInfo()




