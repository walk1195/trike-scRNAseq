#!/usr/bin/Rscript

################################################################################
# Setting consistent color maps for plotting
################################################################################

# ------------------------------------------------------------------------------
# Function to generate vector of n colors
# ------------------------------------------------------------------------------
# Adapted from https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r

generate_color_vector <- function(n, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  all_colors = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
  color_vec <- sample(all_colors, size=n)
  return(color_vec)
}

# ------------------------------------------------------------------------------
# Samples
# ------------------------------------------------------------------------------
sampleCols <- c(
  "trike_01_pretx" = "#DB9D85",
  "trike_01_2w"    = "#86B875",
  "trike_02_pretx" = "#4CB9CC",
  "trike_02_2w"    = "#4A78A8")

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

# HPCA main
hpca_mainCols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#17BECF",
                  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C",
                  "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#BC80BD")

# Scibet
scibetCols <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#17BECF",
               "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#FF7F00", "#FDBF6F", "#E31A1C",
               "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#BC80BD", "#8DD3C7", "#006400")

# Ammons PBMCs
ammonsCols.l1 <- c("B cell"="#1B9E77","CD34+ Unclassified"="#7570B3","CD4 T cell"="#FF7F00","CD8/NK cell"="tomato","Cycling T cell"="orchid3","DC"="#FB9A99",
                   "DN T cell"="#B2DF8A","gd T cell"="#1F78B4","Granulocyte"="#33A02C","Monocyte"="#FDBF6F","Neutrophil"="#E31A1C","Plasma cell"="#B15928")

ammonsCols.l2 <- c("B cell"="#1B9E77","Basophil"="#D95F02","CD34+ Unclassified"="#7570B3","CD4+ Naive"="#E7298A","CD4+ T reg"="#66A61E","CD4+ TCM"="#E6AB02",
                   "CD4+ TEM"="#A6761D","CD8 T cell"="#17BECF","CD8+ Memory"="#A6CEE3","DC"="#1F78B4","DN T cell"="#B2DF8A","Eosinophil"="#33A02C",
                   "gd T cell"="#FB9A99","M-MDSC"="#FF7F00","Monocyte"="#FDBF6F","Neutrophil"="#E31A1C","NK cell"="#CAB2D6","NK T cell"="#6A3D9A",
                   "Other T cell"="#FFFF99","Plasma cell"="#B15928","PMN-MDSC"="cadetblue1")


# ------------------------------------------------------------------------------
# Other categories
# ------------------------------------------------------------------------------

doubletCols <- c('#FF7F00','dodgerblue2')
