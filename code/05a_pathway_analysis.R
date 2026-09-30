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
library(clusterProfiler)
library(org.Cf.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggpubr)


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

# Data 
s1 <- readRDS(file=glue::glue("{objDir}filtered_object.rds"))

# Functions
source("code/utils/functions.R")
source("code/utils/colormaps.R")

# ------------------------------------------------------------------------------
# Read in DGEA results
# ------------------------------------------------------------------------------
de_results <- read.delim("out/Figures/DGEA/deseq2_results.tsv")

# Subset to just CD8/NK cell results
bulk_results_cd8 <- de_results |> dplyr::filter(cell_type == "CD8/NK cell")

# ------------------------------------------------------------------------------
# Generate ranked genes list
      # Ranking genes by -log10(pval) x direction of FC (sign)
      # An alternative ranking metric is the DESeq2 `stat` metric -- but MSI recommends the former
# ------------------------------------------------------------------------------
gene_list <- bulk_results_cd8 |>
  dplyr::filter(!is.na(pvalue), !is.na(log2FoldChange), pvalue > 0) |>
  dplyr::mutate(score = -log10(pvalue) * sign(log2FoldChange)) |>
  dplyr::select(gene, score)

# setting gene names
gene_list <- gene_list$score |> setNames(gene_list$gene)
gene_list <- sort(gene_list, decreasing = TRUE)

# ------------------------------------------------------------------------------
# GSEA: GO terms

# Some parameters to note:
      # minGSSize - minimum # genes a given pathway must have from the ranked gene list to be tested
            # This can be set to avoid small gene sets where enrichment is driven by only a handful of genes (e.g. 5)
      # maxGSSize - max # genes allowed in a GO term
      # pvalueCutoff - threshold to filter results

# using 'fgsea' for GSEA analysis, please cite Korotkevich et al (2019).
# ------------------------------------------------------------------------------

### BP - Biological Processes
          # There are ties in the preranked stats (14.86% of the list).
          # The order of those tied genes will be arbitrary, which may produce unexpected results.
gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Cf.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed=124) # for reproducibility

# Extract results
results_go_bp <- gsea_go_bp@result |>
  dplyr::arrange(p.adjust)

# Read back in 
results_go_bp_old <- read.delim(paste0(resDir, "GO/gsea_go_bp_results.tsv"))
               
# Grabbing top pathways
top_pathways <- results_go_bp|>
  dplyr::slice_head(n = 10) |>
  dplyr::pull(ID)

## Enrichplots --- to show enrichment of the top pathway hits

# Single pathway
for (i in top_pathways) {
  # Plot
  p <- gseaplot2(gsea_go_bp, geneSetID=i, title=gsea_go_bp@result$Description[match(i, gsea_go_bp@result$ID)])
  # Save
  file_name <- gsub(" ", "_", gsea_go_bp@result$Description[match(i, gsea_go_bp@result$ID)])
  ggsave(paste0(resDir, "GO/pathway_", file_name, ".png"), p, dpi=400, height=5, width=7)
}

# Top 3 pathways
p <- gseaplot2(gsea_go_bp, geneSetID = top_pathways[1:3], color=c('green', '#1BB6AFFF', '#EE6100FF'))
ggsave(paste0(resDir, "GO/enrichment_plot_top_3_pathways.png"), p, dpi=400, height=5, width=7)


## Barchart
results_go_bp_log <- results_go_bp |>
  dplyr::filter(ID != "GO:0044419") |>
  dplyr::mutate(neglog10p = -log10(p.adjust))

ggbarplot(results_go_bp_log[1:10,], 
          x = "Description", 
          y = "NES",
          fill = "neglog10p", 
          color = "black",
          sort.val = "asc", 
          rotate = TRUE,
          ggtheme = theme_classic()) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x=NULL) +
  theme(axis.text = element_text(face='bold')) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)") 
#  scale_fill_distiller(palette = "Blues", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)") 
  #scale_y_continuous(limits = c(-2, NA))
#ggsave(paste0(resDir, 'GO/BP/barchart_bp_top10.png'), dpi=400, height=5, width=8)
ggsave(paste0(resDir, 'GO/BP/barchart_bp_top10_reds.png'), dpi=400, height=5, width=8)

# Plotting all the other results (~380) in sets of 50
ggbarplot(results_go_bp_log[100:150,], 
          x = "Description", 
          y = "NES",
          fill = "neglog10p", 
          color = "black",
          sort.val = "asc", 
          rotate = TRUE,
          ggtheme = theme_classic()) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x=NULL) +
  theme(axis.text = element_text(face='bold')) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)") 
ggsave(paste0(resDir, 'GO/BP/barchart_bp_top50.png'), dpi=400, height=5, width=8)


# It looks like there are some interesting pathways that are not necessarily the top 50, but still significantly enriched
        # e.g. cell killing, type II interferon production, etc.
# Utilizing the function simplify() to collapse redundant pathways
gsea_go_bp_simple <- clusterProfiler::simplify(gsea_go_bp, cutoff = 0.5, by = "p.adjust", select_fun = min) # 32 pathways
gsea_go_bp_simple <- clusterProfiler::simplify(gsea_go_bp, cutoff = 0.6, by = "p.adjust", select_fun = min) # 67 pathways
gsea_go_bp_simple <- clusterProfiler::simplify(gsea_go_bp, cutoff = 0.7, by = "p.adjust", select_fun = min) # 147 pathways

# Check results
dim(gsea_go_bp_simple)

# This lowered our results from 358 to 67
# log-transforming p-vals
results_go_bp_log_simple <- gsea_go_bp_simple |>
  dplyr::filter(ID != "GO:0044419") |>
  dplyr::mutate(neglog10p = -log10(p.adjust))

# And plot top hits now
highlight_terms <- c("integrin-mediated signaling pathway", "leukocyte activation", 
                     "cell adhesion", "regulation of apoptotic process")

# Replicate the same ascending NES sort that sort.val = "asc" applies internally
plot_df <- results_go_bp_log_simple[1:20,] |> dplyr::arrange(NES)
label_cols <- ifelse(plot_df$Description %in% highlight_terms, "red", "black")
ggbarplot(results_go_bp_log_simple[1:20,], 
          x = "Description", 
          y = "NES",
          fill = "neglog10p", 
          color = "black",
          sort.val = "asc", 
          rotate = TRUE,
          ggtheme = theme_classic()) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x = NULL) +
  theme(axis.text.y = element_text(face = 'bold', size = 10, color = label_cols),
        axis.text.x = element_text(face = 'bold', size = 10)) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)")
ggsave(paste0(resDir, 'GO/BP/barchart_bp_top20_simplified_0.6.png'), dpi=400, height=5, width=9)


# Plotting enrichplots for key pathways
pathways_of_interest <- c('natural killer cell mediated cytotoxicity', 'cell killing', 'type II interferon production')
pathways_of_interest <- c('cell adhesion')

for (i in pathways_of_interest) {
  id <- gsea_go_bp@result$ID[match(i, gsea_go_bp@result$Description)]
  p <- gseaplot2(gsea_go_bp, geneSetID = id, title = i)
  ggsave(paste0(resDir, "GO/BP/pathway_test", gsub(" ", "_", i), ".png"),
         p, dpi=400, height=5, width=7)
}

### MF - Molecular Function
gsea_go_mf <- gseGO(
  geneList = gene_list,
  OrgDb = org.Cf.eg.db,
  keyType = "SYMBOL",
  ont = "MF",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed=T) # for reproducibility


# Extract results
results_go_mf <- gsea_go_mf@result |>
  dplyr::arrange(p.adjust)

# Grabbing top pathways
top_pathways <- results_go_mf|>
  dplyr::slice_head(n = 10) |>
  dplyr::pull(ID)

## Enrichplots --- to show enrichment of the top pathway hits

# Single pathway
for (i in top_pathways) {
  # Plot
  p <- gseaplot2(gsea_go_mf, geneSetID=i, title=gsea_go_mf@result$Description[match(i, gsea_go_mf@result$ID)])
  # Save
  file_name <- gsub("[ \\/]", "_", gsea_go_mf@result$Description[match(i, gsea_go_mf@result$ID)])
  ggsave(paste0(resDir, "GO/MF/pathway_", file_name, ".png"), p, dpi=400, height=5, width=7)
}

## Barchart
results_go_mf_filtered <- results_go_mf |>
  dplyr::mutate(neglog10p = -log10(p.adjust))

ggbarplot(results_go_mf_filtered[1:10,], 
          x = "Description", 
          y = "NES",
          fill = "neglog10p", 
          color = "black",
          sort.val = "asc", 
          rotate = TRUE,
          ggtheme = theme_classic()) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x=NULL) +
  theme(axis.text = element_text(face='bold')) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)")
  #scale_fill_distiller(palette = "Blues", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)") 
ggsave(paste0(resDir, 'GO/MF/barchart_mf_top10.png'), dpi=400, height=5, width=9)
ggsave(paste0(resDir, 'GO/MF/barchart_mf_top10_reds.png'), dpi=400, height=5, width=9)


### CC - Cellular Compartment
gsea_go_cc <- gseGO(
  geneList = gene_list,
  OrgDb = org.Cf.eg.db,
  keyType = "SYMBOL",
  ont = "CC",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed=T) # for reproducibility

# Extract results
results_go_cc <- gsea_go_cc@result |>
  dplyr::arrange(p.adjust)

# Grabbing top pathways
top_pathways <- results_go_cc|>
  dplyr::slice_head(n = 10) |>
  dplyr::pull(ID)

## Enrichplots --- to show enrichment of the top pathway hits

# Single pathway
for (i in top_pathways) {
  # Plot
  p <- gseaplot2(gsea_go_cc, geneSetID=i, title=gsea_go_cc@result$Description[match(i, gsea_go_cc@result$ID)])
  # Save
  file_name <- gsub(" ", "_", gsea_go_cc@result$Description[match(i, gsea_go_cc@result$ID)])
  ggsave(paste0(resDir, "GO/CC/pathway_", file_name, ".png"), p, dpi=400, height=5, width=7)
}

## Barchart
results_go_cc_filtered <- results_go_cc |>
  dplyr::mutate(neglog10p = -log10(p.adjust))

ggbarplot(results_go_cc_filtered[1:10,], 
          x = "Description", 
          y = "NES",
          fill = "neglog10p", 
          color = "black",
          sort.val = "asc", 
          rotate = TRUE,
          ggtheme = theme_classic()) +
  geom_hline(yintercept = 0, linetype = "solid") +
  labs(x=NULL) +
  theme(axis.text = element_text(face='bold')) +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = guide_colorbar(order = 2), name = "-log10(p.adjust)") 
ggsave(paste0(resDir, 'GO/CC/barchart_cc_top10_red.png'), dpi=400, height=5, width=8)
ggsave(paste0(resDir, 'GO/CC/barchart_cc_top10.png'), dpi=400, height=5, width=8)


### All GO terms at once
gsea_go_all <- gseGO(
  geneList = gene_list,
  OrgDb = org.Cf.eg.db,
  keyType = "SYMBOL",
  ont = "ALL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.01,
  seed=T)

### Results
# Sort by NES, then group
results_go_all <- gsea_go_all@result |>
  dplyr::filter(!is.na(NES), !is.na(p.adjust)) |>
  dplyr::filter(ID != "GO:0044419") |> # removing this go term - vague name
  dplyr::group_by(ONTOLOGY) |>
  dplyr::slice_min(p.adjust, n = 10) |>
  dplyr::ungroup() |>
  dplyr::mutate(nes_sign = ifelse(NES > 0, "pos", "neg")) |>
  dplyr::arrange(dplyr::desc(nes_sign), ONTOLOGY, p.adjust) |>
  dplyr::mutate(Description = factor(Description, levels = rev(Description))) |>
  dplyr::select(-nes_sign)

# Or sort just by group
results_go_all <- gsea_go_all@result |>
  dplyr::filter(!is.na(NES), !is.na(p.adjust)) |>
  dplyr::group_by(ONTOLOGY) |>
  dplyr::slice_min(p.adjust, n = 10) |>
  dplyr::arrange(ONTOLOGY, p.adjust) |>
  dplyr::ungroup() |>
  dplyr::mutate(Description = factor(Description, levels = rev(Description)))

cols <- c(BP = "orange", MF = "seagreen3", CC = "deeppink2")
# Legend on top
ggplot(results_go_all, aes(Description, NES, fill = -log10(p.adjust))) +
  geom_col(color = "white") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_fill_distiller(palette = "Blues", direction = 1, guide = guide_colorbar(order = 2)) +
  scale_y_continuous(limits = range(results_go_all$NES)) +
  new_scale_fill() +
  geom_tile(aes(fill = ONTOLOGY), color = NA, width = 0, height = 0) +
  scale_fill_manual(values = cols, name = "GO Category", guide = guide_legend(order = 1)) +  theme_pubr() +
  theme(
    legend.position = "top",
    legend.location = "plot",
    legend.justification = "right",
    legend.box.just = "right",
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 5, 0),
    axis.text.y = element_text(face = "bold", color = cols[results_go_all$ONTOLOGY[match(levels(results_go_all$Description), results_go_all$Description)]])) +
  labs(y = "NES", fill = "-log10(p.adjust)", x=NULL)
ggsave(paste0(resDir, 'barplot_GO_terms_all.png'), dpi=400, height=8, width=11)
ggsave(paste0(resDir, 'barplot_GO_terms_all_by_group.png'), dpi=400, height=8, width=11)

# legend on side
ggplot(results_go_all, aes(Description, NES, fill = -log10(p.adjust))) +
  geom_col(color = "white") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  scale_y_continuous(limits = range(results_go_all$NES)) +
  new_scale_fill() +
  geom_tile(aes(fill = ONTOLOGY), color = NA, width = 0, height = 0) +
  scale_fill_manual(values = cols, name = "GO Category") +
  theme_classic() +
  theme(axis.text.y = element_text(face = "bold", color = cols[results_go_all$ONTOLOGY[match(levels(results_go_all$Description), results_go_all$Description)]])) +
  labs(y = "NES", fill = "-log10(p.adjust)", x=NULL)
ggsave(paste0(resDir, 'barplot_GO_terms_all_side_legend.png'), dpi=400, height=6, width=10) 
ggsave(paste0(resDir, 'barplot_GO_terms_all_side_legend_by_group.png'), dpi=400, height=6, width=10) 

# Save results dfs
write.table(results_go_all, paste0(resDir, "GO/gsea_go_results_all.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_go_bp, paste0(resDir, "GO/gsea_go_bp_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_go_cc, paste0(resDir, "GO/gsea_go_cc_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_go_mf, paste0(resDir, "GO/gsea_go_mf_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# GSEA: KEGG terms
# ------------------------------------------------------------------------------

### First we have to convert from gene SYMBOL to ENTREZ ID
gene_df <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Cf.eg.db) # ~1.94% of genes were dropped

# Reformatting
gene_list_kegg <- gene_list[gene_df$SYMBOL] # our gene list was already sorted by log(p.value) x sign(FC)
names(gene_list_kegg) <- gene_df$ENTREZID 

# Run GSEA
gsea_kegg <- gseKEGG(
  geneList=gene_list_kegg,
  organism="cfa",
  minGSSize=10,
  maxGSSize=500,
  pvalueCutoff=0.05,
  seed=TRUE)

# Plot results
results_kegg <- gsea_kegg@result |>
  dplyr::filter(!is.na(NES), !is.na(p.adjust))

kegg_top_20 <- results_kegg |>
  dplyr::slice_max(abs(NES), n = 50) |>
  dplyr::mutate(neglog10p = -log10(p.adjust), Description = reorder(Description, NES))

ggplot(kegg_top_20, aes(Description, NES, fill = neglog10p)) +
  geom_col(color = "black") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = range(results_kegg$NES)) +
  scale_fill_distiller(palette = "YlOrRd", direction = -1) +
  theme_pubr() +
  theme(axis.text.y = element_text(face = "bold")) +
  labs(x = NULL, y = "NES", fill = "-log10(p.adjust)")
ggsave(paste0(resDir, 'KEGG/barplot_KEGG_terms_top20.png'), dpi=400, height=7, width=9) 
ggsave(paste0(resDir, 'KEGG/barplot_KEGG_terms_top50.png'), dpi=400, height=14, width=9) 

# Save results df
write.table(results_kegg, paste0(resDir, "KEGG/gsea_kegg_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)


# ------------------------------------------------------------------------------
# GSEA: Hallmark
# ------------------------------------------------------------------------------

# Grab hallmark gene sets
hallmark <- msigdbr(species = "Canis lupus familiaris", collection = "H")

hallmark_sets <- hallmark |>
  dplyr::select(gs_name, gene_symbol)

# Run hallmark
gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_sets,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed = TRUE)

# Plot results
results_hallmark <- gsea_hallmark@result |>
  dplyr::filter(!is.na(NES), !is.na(p.adjust))

hallmark_top_20 <- results_hallmark |>
  dplyr::slice_max(abs(NES), n = 20) |>
  dplyr::mutate(neglog10p = -log10(p.adjust), Description = reorder(Description, NES))

hallmark_top_10 <- results_hallmark |>
  dplyr::slice_max(abs(NES), n = 10) |>
  dplyr::mutate(neglog10p = -log10(p.adjust), Description = reorder(Description, NES))

ggplot(hallmark_top_20, aes(Description, NES, fill = neglog10p)) +
  geom_col(color = "black") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = range(results_kegg$NES)) +
  scale_fill_distiller(palette = "YlOrRd", direction = -1) +
  theme_pubr() +
  theme(axis.text.y = element_text(face = "bold")) +
  labs(x = NULL, y = "NES", fill = "-log10(p.adjust)")
ggsave(paste0(resDir, 'Hallmark/barplot_hallmark_top20.png'), dpi=400, height=7, width=9) 

# Enrichplots
for (i in hallmark_top_10$ID) {
  # Plot
  p <- gseaplot2(gsea_hallmark, geneSetID=i, title=gsea_hallmark@result$Description[match(i, gsea_hallmark@result$ID)])
  # Save
  file_name <- gsub(" ", "_", gsea_hallmark@result$Description[match(i, gsea_hallmark@result$ID)])
  ggsave(paste0(resDir, "Hallmark/pathway_", file_name, ".png"), p, dpi=400, height=5, width=7)
}

# Save results df
write.table(results_hallmark, paste0(resDir, "Hallmark/gsea_hallmark_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)



# ------------------------------------------------------------------------------
# GSEA: REACTOME
# ------------------------------------------------------------------------------

# Download REACTOME files
# Pathways
reactome_pathways <- read.delim("https://reactome.org/download/current/ReactomePathways.txt", header=FALSE, col.names=c("ID","Description","Species"))

# Mapping file
        # Using the 'All levels' version to test for enrichment at various hierarchy levels. Details about the different mapping files at
        # https://reactome.org/download-data
reactome_map_ensembl <- read.delim("https://reactome.org/download/current/Ensembl2Reactome_All_Levels.txt", header=FALSE, col.names=c("ENSEMBL","ID","URL","Event","Evidence","Species"))

# Subset to canine only
reactome_pathways_cfa <- reactome_pathways |>
  dplyr::filter(Species == "Canis familiaris")
reactome_map_cfa_ensembl <- reactome_map_ensembl |>
  dplyr::filter(Species == "Canis familiaris")

# Look's like these are ENSEMBL ID format -- converting gene list from SYMBOL to ENSEMBL 
gene_df <- bitr(names(gene_list), fromType="SYMBOL", toType="ENSEMBL", OrgDb=org.Cf.eg.db) # ~12.44% of genes were dropped
gene_list_reactome <- gene_list[gene_df$SYMBOL]
names(gene_list_reactome) <- gene_df$ENSEMBL # This is already sorted from a previous run

# Build reactome gene sets -- 
reactome_sets <- reactome_map_cfa_ensembl |>
  dplyr::inner_join(reactome_pathways_cfa, by="ID") |>
  dplyr::select(ID, Description, ENSEMBL)

# Convert to gene symbols
    # Note that ~376 ENSEMBL IDs map to multiple gene symbols. We are ignoring this for now. I need to check how this conversion handled those duplicates???
reactome_sets_symbol <- reactome_sets |>
  dplyr::inner_join(gene_df, by = "ENSEMBL") |>
  dplyr::select(ID, Description, SYMBOL) |>
  dplyr::distinct()

# Run GSEA
gsea_reactome <- GSEA(
  geneList = gene_list,
  TERM2GENE = reactome_sets_symbol |> dplyr::select(ID, SYMBOL),
  TERM2NAME = reactome_sets_symbol |> dplyr::select(ID, Description),
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  seed = 124)

# Results
results_reactome <- gsea_reactome@result |>
  dplyr::filter(!is.na(NES), !is.na(p.adjust))

reactome_top_20 <- results_reactome |>
 # dplyr::filter(NES > 0) |>
  dplyr::slice_max(NES, n = 20) |>
  dplyr::mutate(
    neglog10p = -log10(p.adjust),
    Description = reorder(Description, NES))

highlight_terms <- c("Interferon Signaling", "Cytokine Signaling in Immune system", "MAPK targets/ Nuclear events mediated by MAP kinases",
                     "Signaling by Interleukins")  # <- replace with your terms
label_cols <- ifelse(levels(factor(reactome_top_20$Description)) %in% highlight_terms, "red", "black")

ggplot(reactome_top_20, aes(Description, NES, fill = neglog10p)) +
  geom_col(color = "black") +
  geom_hline(yintercept = 0) +
  coord_flip() +
  scale_y_continuous(limits = range(results_reactome$NES)) +
  scale_fill_distiller(palette = "YlOrRd", direction = -1) +
  theme_pubr() +
  theme(
    legend.position = "top",
    legend.location = "plot",
    legend.justification = "right",
    axis.text.y = element_text(face = "bold", color = label_cols)
  ) +
  labs(x = NULL, y = "NES", fill = "-log10(p.adjust)")
ggsave(paste0(resDir, 'Reactome/barplot_reactome_top20.png'), dpi=400, height=7, width=10) 





