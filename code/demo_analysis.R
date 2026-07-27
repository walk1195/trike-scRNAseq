#!/usr/bin/env Rscript

# Setup paths
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- dirname(normalizePath(script_path))
proj_dir <- dirname(script_dir)
out_dir <- file.path(proj_dir, "code_out", "demo_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

# ---------------------------------------------------------------------
# Basic ggplot2
# ---------------------------------------------------------------------

# Load packages
library(tidyverse)

# Generate sample data
set.seed(42)
sample_data <- tibble(
    sample_id = paste0("Sample_", 1:20),
    condition = rep(c("Control", "Treatment"), each = 10),
    expression = c(rnorm(10, mean = 5, sd = 1.5), 
                   rnorm(10, mean = 8, sd = 1.5))
)

# Save data
write_csv(sample_data, "sample_data.csv")

# Create and save plot
p <- ggplot(sample_data, aes(x = condition, y = expression, fill = condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = "Gene Expression: Control vs Treatment",
         x = "Condition", y = "Expression Level") +
    scale_fill_brewer(palette = "Set2")

ggsave("expression_boxplot.png", p, width = 8, height = 6, dpi = 300)

# ---------------------------------------------------------------------
# Seurat analysis section
# ---------------------------------------------------------------------

# Configure reticulate to use conda Python (for Leiden algorithm)
library(reticulate)
reticulate_python <- Sys.getenv("RETICULATE_PYTHON")
if (reticulate_python != "") {
    use_python(reticulate_python, required = TRUE)
}

library(Seurat)

# Load bundled pbmc_small dataset
data("pbmc_small")

# Update to Seurat v5 format and convert assay to v5
pbmc_small <- UpdateSeuratObject(pbmc_small)
pbmc_small[["RNA"]] <- as(pbmc_small[["RNA"]], Class = "Assay5")

pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10)
# Clustering with Leiden algorithm
pbmc_small <- FindClusters(pbmc_small, algorithm = 4, resolution = 0.5)
pbmc_small <- RunUMAP(pbmc_small, dims = 1:10)

# Create UMAP plot and save as PDF
pdf("umap_leiden_clusters.pdf", width = 8, height = 6)
p_umap <- DimPlot(pbmc_small, reduction = "umap", label = TRUE) + 
    ggplot2::labs(title = "UMAP of Leiden Clusters")
print(p_umap)
dev.off()

# Summary of clusters
cat("Clustering complete. Cluster distribution:\n")
print(table(Idents(pbmc_small)))

cat("Done! Results in:", out_dir, "\n")
