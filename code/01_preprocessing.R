#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Initial Preprocessing & QC
### Samples: Vallera 1, 2, 3, 4
### Author: Grace Walker
### Date: July 27, 2026
#############################################################################################
#############################################################################################

# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse)
library(Seurat)
library(patchwork)

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
resDir <- 'out/Figures/Preprocessing/'

# List of samples
sample_ids <- readLines(paste0(projDir,'/inputs/cellranger_ids.txt'))

# New IDs
new_sample_ids <- readLines(paste0(projDir,'/inputs/samples.txt'))

# ------------------------------------------------------------------------------
# Generate Seurat objects
# ------------------------------------------------------------------------------
objects <- c()
for (i in 1:length(sample_ids)) {
  
  sample_id <- sample_ids[[i]]
  
  print(glue::glue('Generating seurat object for {sample_id}'))
  
  umgc_id <- paste0('cellranger_', sample_id, '_GEX_FL/')
  fullDir <- paste0(dataDir, umgc_id, 'filtered_feature_bc_matrix/')
  counts <- Read10X(data.dir = fullDir)
  
  new_sample_id <- new_sample_ids[[i]]
  
  s1 <- CreateSeuratObject(counts = counts, project = new_sample_id, min.cells = 3, min.features = 200)
  
  # Calculate MT percent
  s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
  
  objects[[new_sample_id]] <- s1
}


# ------------------------------------------------------------------------------
# Generate QC plots for each sample
# ------------------------------------------------------------------------------

for (i in 1:length(objects)) {
  # Get object
  s1 <- objects[[i]]
  sample_id <- names(objects)[i]
  # Set curr res dir
  curr_dir <- paste0(resDir, sample_id, '/')
  # Create subdir
  dir.create(curr_dir, recursive = T)
  # Plot  
  generate_qc_plots(s1, sample_id=sample_id, resDir=curr_dir, feature_min=350, count_min=1000, mt_threshold = 15, sample_colors=sampleCols)
}

# ------------------------------------------------------------------------------
# Generate QC plots across all samples
# ------------------------------------------------------------------------------

# Setting custom order of samples
custom_order <- c('trike_01_pretx', 'trike_01_2w', 'trike_02_pretx', 'trike_02_2w')

# Plot
all_sample_qc_plots(objects, order=custom_order, custom_cols=sampleCols, cell_count=T, nFeat=T, nCount=T, mt=T, density=T)

# Testing just density plot
all_sample_qc_plots(objects, order=custom_order, custom_cols=sampleCols, density=T)

# ------------------------------------------------------------------------------
# Filter each dataset
# ------------------------------------------------------------------------------
feature_min = 350
count_min = 1000
mt_threshold = 15

for (i in 1:length(objects)) {
  # Get obj
  s1 <- objects[[i]]
  sample_id <- names(objects)[i]
  
  # Prefilter count
  prefilter_count = length(Cells(s1))
  
  # Filter obj
  s1@meta.data$keep <- with(s1@meta.data, ifelse(nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt_threshold, TRUE, FALSE))
  s1 <- subset(s1, subset = keep == TRUE)
  
  # Post filter count
  postfilter_count = length(Cells(s1))
  
  total = prefilter_count - postfilter_count
  print(glue::glue("Cells removed from {sample_id} : {total}"))
  
  # Save
  saveRDS(s1, glue::glue("{objDir}{sample_id}_filtered.rds"))
}

# ------------------------------------------------------------------------------
# Mark doublets (not removing yet)
# ------------------------------------------------------------------------------

# Grab no. of cells
cellNum <- dim(s1)[2]
percentDoub <- ((8e-6)*cellNum)+0.0005 # linear equation I calculated from 10x genomics table
nExp_poi <- round(percentDoub*nrow(s1@meta.data)) 
s1 <- doubletFinder(s1, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
colnames(s1@meta.data)[length(colnames(s1@meta.data))] <- "DF.classifications"
Idents(s1) <- "DF.classifications"
DF.name <- colnames(s1@meta.data)[grepl("DF.classification", colnames(s1@meta.data))]


# ------------------------------------------------------------------------------
# Merge samples into single obj
# ------------------------------------------------------------------------------

# Perform merge
s1.merged <- merge(x = objects[[1]], y = objects[-1], add.cell.ids = names(objects))

# Save
saveRDS(s1.merged, file=glue::glue("{objDir}merged_object.rds"))


########### Session info ###########
sessionInfo()
