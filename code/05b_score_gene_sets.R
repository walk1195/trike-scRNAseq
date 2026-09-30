#!/usr/bin/env Rscript

#############################################################################################
############################################################################################
### TriKE scRNAseq: Gene Set Scoring and Individual Gene Plots
### Samples: n=4 (see README)
### Author: Grace Walker
### Date: September 30, 2026
#############################################################################################
#############################################################################################

# -------------------------------------------------------------------------------------------
# Set up environment
# -------------------------------------------------------------------------------------------

# Libs
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(clusterProfiler)
library(org.Cf.eg.db)
library(enrichplot)
library(ggpubr)


# Params
options(future.globals.maxSize = 4 * 1024^3)  # Set to 4 GB
gc() # free up memory

# Dirs
projDir <- getwd()
dataDir <- 'inputs/data/' # symlinked from data_delivery
objDir <- 'out/Objects/'
resDir <- 'out/Figures/Gene_scoring/'

if (!dir.exists(resDir)) {
  dir.create(resDir, recursive=TRUE)
}

# Data 
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))

# Functions
source("code/utils/functions.R")
source("code/utils/colormaps.R")


# ------------------------------------------------------------------------------
# Scoring Gene Sets
# ------------------------------------------------------------------------------
# Checking if viral gene set is present in a certain subset of cells
genes <- AnnotationDbi::select(org.Cf.eg.db, keys = "GO:0045069", keytype = "GOALL", columns = c("SYMBOL", "GOALL"))
genes <- AnnotationDbi::select(org.Cf.eg.db, keys = "GO:0032609", keytype = "GOALL", columns = c("SYMBOL", "GOALL"))

viral_gs <- list(genes$SYMBOL)
ifn_gs <- list(genes$SYMBOL)

# Calculate module score
s1 <- AddModuleScore(s1, features=ifn_gs, name='Interferon_Production')

# Plot
FeaturePlot(s1, features='Interferon_Production1', reduction='umap.harmony')

# Compare by timepoint
s1$timepoint <- ifelse(grepl("_pretx$", s1$orig.ident), "pretx", ifelse(grepl("_2w$", s1$orig.ident), "2w", NA))
s1$timepoint <- factor(s1$timepoint, levels = c("pretx", "2w"))




# ------------------------------------------------------------------------------
# Plotting CPMs of individual genes
# ------------------------------------------------------------------------------

# Join layers
s1.joined <- JoinLayers(s1)

# Genes of interest
genes_list <- c("TNF", "IFNG", "CSF2", "CCL4", "CCL5", "IL2RA", "MKI67", "TOP2A")

# Extract counts
counts <- GetAssayData(s1.joined, assay = "RNA", layer = "counts")
cpm <- t(t(counts) / colSums(counts)) * 1e6

# For each gene, generate boxplots of log expr between timepoints (by cell type)
for (gene in genes_list) {
  # Prepare df for plotting
  plot_df <- data.frame(
    CPM = cpm[gene, ],
    ammons = s1$ammons.l1,
    encode = s1$blue_raw,
    timepoint = s1$timepoint) |>
    dplyr::mutate(logCPM = log2(CPM + 1))
  
  # By Ammons labels
  ggplot(plot_df, aes(x = timepoint, y = logCPM, fill = timepoint)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ammons) +
    scale_x_discrete(labels = c("pretx" = "PreTx", "2w" = "2w")) +
    scale_fill_manual(values = c("pretx" = "lightblue", "2w" = "darkblue"), labels = c("pretx" = "PreTx", "2w" = "2w")) +
    theme_classic() +
    labs(title = gene, x = NULL, y = "log2(CPM + 1)", fill = "Timepoint") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 11),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(size=12),
          legend.title = element_text(face = "bold", size = 13),
          legend.text = element_text(size = 14))
  ggsave(paste0(resDir, "ammons_boxplot_log_", gene, ".png"), dpi=400, height=7, width=9)
  
  # By encode labels
  ggplot(plot_df, aes(x = timepoint, y = logCPM, fill = timepoint)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~encode) +
    scale_x_discrete(labels = c("pretx" = "PreTx", "2w" = "2w")) +
    scale_fill_manual(values = c("pretx" = "lightblue", "2w" = "darkblue"), labels = c("pretx" = "PreTx", "2w" = "2w")) +
    theme_classic() +
    labs(title = gene, x = NULL, y = "log2(CPM + 1)", fill = "Timepoint") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
          strip.text = element_text(face = "bold", size = 11),
          axis.title = element_text(face = "bold", size = 13),
          axis.text = element_text(size=12),
          legend.title = element_text(face = "bold", size = 13),
          legend.text = element_text(size = 14))
  ggsave(paste0(resDir, "encode_boxplot_log_", gene, ".png"), dpi=400, height=8, width=9)
}

# Now plot for each dog individually
for (gene in genes_list) {
  plot_df <- data.frame(CPM=cpm[gene, ], ammons=s1$ammons.l1, encode=s1$blue_raw, orig.ident=s1$orig.ident) |>
    dplyr::mutate(
      logCPM=log2(CPM+1),
      dog=sub("_(pretx|2w)$","",orig.ident),
      timepoint=factor(sub("^.*_(pretx|2w)$","\\1",orig.ident), levels=c("pretx","2w")))
  
  for (group in c("ammons", "encode")) {
    for (d in c("trike_01", "trike_02")) {
      p <- ggplot(dplyr::filter(plot_df, dog == d), aes(timepoint, logCPM, fill=timepoint)) +
        geom_boxplot(outlier.size=.5) +
        facet_wrap(as.formula(paste("~", group))) +
        scale_x_discrete(labels=c(pretx="PreTx", `2w`="2w")) +
        scale_fill_manual(
          values=c(pretx="lightblue", `2w`="darkblue"),
          labels=c(pretx="PreTx", `2w`="2w")) +
        theme_classic() +
        labs(title=glue::glue("{gene} ({d})"), x=NULL, y="log2(CPM + 1)", fill="Timepoint") +
        theme(plot.title=element_text(hjust=.5, face="bold", size=16),
          strip.text=element_text(face="bold", size=11),
          axis.title=element_text(face="bold", size=13),
          axis.text=element_text(size=12),
          legend.title=element_text(face="bold", size=13),
          legend.text=element_text(size=14))
      ggsave(paste0(resDir, group, "_boxplot_log_", gene, "_", sub("_","",d), ".png"), p, dpi=400, height=7, width=9)
    }
  }
}






