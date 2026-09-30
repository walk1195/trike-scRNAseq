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

# Barplot of contribution by sample
plot_df <- s1@meta.data |>
  dplyr::count(ammons.l1, orig.ident) |>
  group_by(ammons.l1) |>
  mutate(prop = n / sum(n))

ggplot(plot_df, aes(x = ammons.l1, y = prop, fill = orig.ident)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = sampleCols) +
  labs(x = NULL, y = "Relative contribution", fill = "Sample") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plot_df <- as.data.frame(table(s1$ammons.l1, s1$orig.ident))
colnames(plot_df) <- c("ammons.l1", "orig.ident", "n")
plot_df$prop <- ave(plot_df$n, plot_df$orig.ident, FUN = function(x) x / sum(x))
plot_df$orig.ident <- factor(plot_df$orig.ident, levels = c("trike_01_pretx", "trike_01_2w", "trike_02_pretx", "trike_02_2w"))

ggplot(plot_df, aes(x = orig.ident, y = prop, fill = ammons.l1)) +
  geom_col() +
  theme_classic() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = ammonsCols.l1) +
  theme(axis.text = element_text(size=12), axis.title=element_text(size=14, face='bold'), legend.text = element_text(size=12),
        axis.text.x = element_text(angle=45, size=10, hjust=1)) +
  labs(x = NULL, y = "Cell type composition", fill = NULL)
ggsave(paste0(resDir, 'ammons_l1_proportions_barplot.png'), dpi=400, height=6, width=8)


#### paired
plot_df <- as.data.frame(table(s1$ammons.l1, s1$orig.ident))
colnames(plot_df) <- c("ammons.l1", "orig.ident", "n")
plot_df$prop <- ave(plot_df$n, plot_df$orig.ident, FUN = function(x) x / sum(x))

# Get timpoint col
plot_df$timepoint <- factor(sub("^.*_(pretx|2w)$", "\\1", plot_df$orig.ident), levels = c("pretx", "2w"))
plot_df$timepoint <- factor(plot_df$timepoint, levels = c("pretx", "2w"), labels = c("PreTx", "2w"))

# Fixing CD34+ label that's too big
plot_df$ammons.l1 <- ifelse(as.character(plot_df$ammons.l1) == "CD34+ Unclassified",
  "CD34+\nUnclassified", as.character(plot_df$ammons.l1))

# Plot
ggplot(plot_df, aes(x = timepoint, y = prop, group = dog, color = dog)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~ammons.l1, ncol = 4) +
  scale_y_continuous(labels = scales::percent) +
  scale_color_manual(values = c("trike_01" = "royalblue3", "trike_02" = "darkorange1")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 12)) +
  labs(x = NULL, y = "Cell type proportion", color = "Dog")
ggsave(paste0(resDir,'ammons_proportions_paired_dotplot.png'), dpi=400, height=5, width=7)

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------
saveRDS(s1, file=glue::glue("{objDir}filtered_object.rds"))

##### Session info
sessionInfo()


