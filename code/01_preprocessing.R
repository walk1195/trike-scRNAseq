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

# Dirs
projDir <- getwd()
dataDir <- 'inputs/data/' # symlinked from data_delivery
objDir <- 'out/Objects/'
resDir <- 'out/Figures/Preprocessing/'

# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Functions
source("code/utils/functions.R")
source("code/utils/colormaps.R")

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
  objects[[new_sample_id]] <- s1
}

all_obj <- objects # store for later
objects <- all_obj
# ------------------------------------------------------------------------------
# Generate QC plots for each sample
# ------------------------------------------------------------------------------

# Testing on single sample
objects <- objects[1]

for (i in 1:length(objects)) {
  # Get object
  s1 <- objects[[i]]
  sample_id <- names(objects)[i]
  # Set curr res dir
  curr_dir <- paste0(resDir, sample_id, '/')
  # Create subdir
  dir.create(curr_dir, recursive = T)
  # Plot  
  generate_qc_plots(s1, sample_id=sample_id, resDir=curr_dir, feature_min=300, count_min=500, mt_threshold = 20, sample_colors=sampleCols)
}


###### Filter each dataset #######
s1 <- objects[[1]]

prefilter_count = length(Cells(s1))
s1@meta.data$keep <- with(s1@meta.data, ifelse(nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt_threshold, TRUE, FALSE))
s1 <- subset(s1, subset = keep == TRUE)
# Count
postfilter_count = length(Cells(s1))
total = prefilter_count - postfilter_count
print(paste("Cells removed:", total))

# Save
saveRDS(s1, paste0(objDir,"/ORv03_PreTx_processed.rds"))
saveRDS(s1, paste0(objDir,"/ORv03_D10_processed.rds"))


# Save session info
sessionInfo()
