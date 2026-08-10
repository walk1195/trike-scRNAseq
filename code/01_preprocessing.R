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
library(scDblFinder)

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
# To read in individual samples if iterating back over QC
objects <- c()
for (obj in list.files(objDir)) {
  s1 <- readRDS(file=glue::glue("{objDir}{obj}"))
  sample_id <- gsub('_filtered.rds', '', obj)
  objects[[sample_id]] <- s1
}
objects$merged_object.rds <- NULL
# ------------------------------------------------------------------------------

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

#### Calculate multiplet rate ####

# For 10x data, scDblFinder automatically estimates the doublet rate parameter (dbr) on a 0.08% per 1000 cells basis
# This is not up to date with GEM-X technology, so setting dbr manually:

# 10x Chromium GEM-X Single Cell 3' Reagent Kits v4 User Guide:
  # https://cdn.10xgenomics.com/image/upload/v1725314293/support-documents/CG000731_ChromiumGEM-X_SingleCell3v4_UserGuide_RevB.pdf
  # pg 19-20 shows the multiplet rate table
  # for 20k cell recovery, expected doublet rate is 8% 

# Call doublet finder function (per sample)
for (sample_id in names(objects)) {
  # Get obj
  s1 <- objects[[sample_id]]
  # Call function
  s1 <- run_doublet_finder(s1, multiplet_rate = 0.08)
  
  # Plot results
  curr_dir <- paste0(resDir, sample_id, '/')
  
  p <- VlnPlot(s1, features=c('nFeature_RNA', 'nCount_RNA'), split.by='scDblFinder.class', pt.size = 0)
  ggsave(paste0(curr_dir,'doublet_finder_vln_plots.png'), plot=p, dpi=400, height=6, width=9)
  
  # Update list
  objects[[sample_id]] <- s1
}

# ------------------------------------------------------------------------------
# Merge samples into single obj
# ------------------------------------------------------------------------------

# Perform merge
s1.merged <- merge(x = objects[[1]], y = objects[-1], add.cell.ids = names(objects))

########### Save progress ###########

saveRDS(s1.merged, file=glue::glue("{objDir}merged_object.rds"))

########### Session info ###########
sessionInfo()
