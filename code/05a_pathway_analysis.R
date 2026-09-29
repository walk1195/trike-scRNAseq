#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: GSEA and ORA Pathway Analysis
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: September 28, 2026
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
resDir <- 'out/Figures/Pathway_Analysis/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# ------------------------------------------------------------------------------
# Read in data
# ------------------------------------------------------------------------------
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))






