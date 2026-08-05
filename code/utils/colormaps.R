#!/usr/bin/Rscript

################################################################################
# Setting consistent color maps for plots
################################################################################

# ------------------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------------------
#samples <- readLines("inputs/samples.txt")

set.seed(123)
#sampleCols <- setNames(grDevices::hcl.colors(length(samples), palette = "Dynamic"), samples)

sampleCols <- c(
  "trike_01_pretx" = "#DB9D85",
  "trike_01_2w"    = "#86B875",
  "trike_02_pretx" = "#4CB9CC",
  "trike_02_2w"    = "#4A78A8"
)

# ------------------------------------------------------------------------------
# Automated annotations
# ------------------------------------------------------------------------------
# Immgen
immgenCols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#17BECF",
                   "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C",
                   "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")


# Blue encode
encodeCols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#17BECF",
                 "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C",
                 "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#BC80BD", "#8DD3C7")


# Scibet
scibetCols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#17BECF",
               "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C",
               "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#BC80BD", "#8DD3C7", "#006400")

