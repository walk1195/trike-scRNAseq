#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Finalizing cell typing for downstream analysis
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: September 22, 2026
#############################################################################################
#############################################################################################

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

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))

#############################################################################################

# We can start to form initial broad annotations based on some very distinct clusters
cell_types <- c(
  "0" = "CD4+ T cells",
  "1" = "Monocytes",
  "2" = "B cells",
  "3" = "CD8+/NK T cells",
  "4" = "DN T cells?",
  "5" = "Eosinophils",
  "6" = "CD8+/NK T cells",
  "7" = "Neutrophils",
  "8" = "DCs",
  "9" = "Plasma cells",
  "10" = "CD8+/NK T cells",
  "11" = "gd T cells")

s1$cell_types_broad <- unname(cell_types[as.character(s1$res0.06)])

DimPlot(s1, group.by='cell_types_broad', reduction='umap.harmony')
DimPlot(s1, group.by='ammons.l1', reduction='umap.harmony')

# Dotplot of main markers
cell_markers <- c('JCHAIN', 'CD79B', 'CD3E', 'CD4', 'CD8A', 'GZMA', 'GATA3', 'S100A12', 'S100A9', 'DLA-DRA',
                  'CCL5', 'CCL4', 'IL2R', 'IL15RA')

DotPlot(s1, features = cell_markers,
        group.by = "cell_types_broad", cluster.idents=T,
        dot.scale = 6) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(axis.text=element_text(size=8)) +
  RotatedAxis()


#############################################################################################
# A given gene across all cell groups
VlnPlot(s1, features='CD8A', group.by='cell_types_broad')
VlnPlot(s1, features='CCL5', group.by='cell_types_broad')
VlnPlot(s1, features='CCL4', group.by='cell_types_broad')
VlnPlot(s1, features='IL2RA', group.by='cell_types_broad')
VlnPlot(s1, features='IL2RB', group.by='cell_types_broad')

VlnPlot(s1, features='CD79B', group.by='cell_types_broad')
VlnPlot(s1, features='CD4', group.by='cell_types_broad')
VlnPlot(s1, features='S100A8', group.by='cell_types_broad')
VlnPlot(s1, features='S100A12', group.by='cell_types_broad')


# Adding in a split by sample, specific for cell types
Idents(s1) <- 'cell_types_broad'

VlnPlot(s1, features='IL2RB', group.by='cell_types_broad', split.by = 'orig.ident', cols=sampleCols, idents=c('CD8+/NK T cells','DN T cells'))
VlnPlot(s1, features='IL2RA', group.by='cell_types_broad', split.by = 'orig.ident', cols=sampleCols)
VlnPlot(s1, features='CCL5', group.by='cell_types_broad', split.by = 'orig.ident', cols=sampleCols, idents='CD4+ T cells')


df <- FetchData(s1, vars = c("IL2RB", "cell_types_broad", "orig.ident"))
ggplot(df[df$cell_types_broad %in% c("CD8+/NK T cells", "DN T cells"), ],
       aes(cell_types_broad, IL2RB, fill = orig.ident)) +
  geom_boxplot() +
  scale_fill_manual(values = sampleCols) +
  labs(y = "IL2RB CPM", fill = "Sample") +
  theme_classic()

##### Save
saveRDS(s1, file=glue::glue("{objDir}filtered_object.rds"))

##### Session info
sessionInfo()


