#!/usr/bin/Rscript

# Load libs
library(Seurat)
library(tidyverse)
library(scales)
library(EnhancedVolcano)

# ==============================================================================
# QC Plots - per sample
# ==============================================================================
generate_qc_plots <- function(seu.obj, sample_id, resDir, feature_min, count_min, mt_threshold, sample_colors) {
  
  # Calculate MT percent
  seu.obj[["percent.mt"]] <- PercentageFeatureSet(seu.obj, pattern = "^MT-")
  
  # QC plots
  VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), layer="counts", ncol = 3)
  ggsave(paste0(resDir,'qc_violin_plots.png'), dpi=400, width=10, height=6)
  
  # Scatter plots
  plot1 <- FeatureScatter(seu.obj, feature1 = "percent.mt", feature2 = "nFeature_RNA") + ggtitle('MT expression x nFeature_RNA')
  plot2 <- FeatureScatter(seu.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle('nCount_RNA x nFeature_RNA')
  scatter_plots <- plot1 + plot2
  ggsave(paste0(resDir,'qc_scatter_plots.png'), plot=scatter_plots, dpi=400, width=10, height=6)
  
  # Density plots
  p1 <- seu.obj@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + 
    scale_x_log10() + 
    scale_fill_manual(values = sample_colors) +
    scale_color_manual(values = sample_colors) +
    geom_vline(xintercept = feature_min,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position = "none") +
    ggtitle("nFeature")
  p2 <- seu.obj@meta.data %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill= orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() + 
    scale_x_log10() + 
    scale_fill_manual(values = sample_colors) +
    scale_color_manual(values = sample_colors) +
    geom_vline(xintercept = count_min,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    theme(legend.position = "none") +
    ggtitle("nCount")
  p3 <-seu.obj@meta.data %>% 
    ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10(labels = label_comma()) + 
    scale_fill_manual(values = sample_colors) +
    scale_color_manual(values = sample_colors) +
    geom_vline(xintercept = mt_threshold,color="red",linetype="dotted") +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("MT % Expression") +
    guides(color = guide_legend(title = "Sample"), 
           fill = guide_legend(title = "Sample"))
  density_plots <- p1 + p2 + p3
  ggsave(paste0(resDir,'qc_density_plots.png'), plot=density_plots, dpi=400, width=15, height=6)
}

# ==============================================================================
# QC plots - all samples
# ==============================================================================

all_sample_qc_plots <- function(obj_list, order, custom_cols, cell_count=F, nFeat=F, nCount=F, mt=F, density=F) {
  # Initialize list to store plots
  plots <- list()
  # Grab qc data
  qc_df <- do.call(rbind, lapply(obj_list, \(x) x@meta.data))
  qc_df$sample <- rep(names(obj_list), sapply(obj_list, ncol))
  
  # Set sample order
  if (!is.null(order)) {qc_df$sample <- factor(qc_df$sample, levels = order)}
  
  ### (1) Raw cell counts barplot
  if (cell_count) {
    cell_counts <- as.data.frame(table(qc_df$sample))
    colnames(cell_counts) <- c("sample", "count")
    p.cell_counts <- ggplot(cell_counts, aes(x=sample, y=count, fill=sample)) +
      geom_col(color='black') +
      scale_fill_manual(values = custom_cols) +
      theme_classic() +
      labs(title='Cell Counts (unfiltered)', x="Sample", y='# of Cells') +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
            axis.title = element_text(face='bold'),
            plot.title = element_text(hjust=0.5, face='bold'),
            legend.position = 'none')
    plots$cell_count <- p.cell_counts
    ggsave(paste0(resDir,'all_samples_counts_unfiltered.png'), plot=p.cell_counts, dpi=400, width=5, height=6)
  }
  
  ### (2) nFeature_RNA violin plot
  if (nFeat) {
    p.nFeat <- ggplot(qc_df, aes(x=sample, y=nFeature_RNA, fill=sample)) +
      geom_violin(scale="width") +
      scale_fill_manual(values = custom_cols) +
      theme_classic() +
      labs(title='nFeature_RNA', x="Sample", y="# Features") +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.title = element_text(face='bold'),
            plot.title = element_text(hjust=0.5, face='bold'))
    plots$nFeat <- p.nFeat
    ggsave(paste0(resDir,'all_samples_nFeature.png'), plot=p.nFeat, dpi=400, width=6, height=6)
  }
  
  ### (3) nCount_RNA violin plot
  if (nCount) {
    p.nCount <- ggplot(qc_df, aes(x=sample, y=nCount_RNA, fill=sample)) +
      geom_violin(scale="width") +
      scale_fill_manual(values = custom_cols) +
      theme_classic() +
      labs(title='nCount_RNA', x="Sample", y="# UMIs") +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.title = element_text(face='bold'),
            plot.title = element_text(hjust=0.5, face='bold'))
    plots$nCount <- p.nCount
    ggsave(paste0(resDir,'all_samples_nCount.png'), plot=p.nCount, dpi=400, width=6, height=6)
  }

  ### (4) MT Expression violin plot
  if (mt) {
    p.mt <- ggplot(qc_df, aes(x=sample, y=percent.mt, fill=sample)) +
      geom_violin(scale="width") +
      scale_fill_manual(values = custom_cols) +
      theme_classic() +
      labs(title='nCount_RNA', x="Sample", y="# UMIs") +
      theme(axis.text.x = element_text(angle=45, hjust=1),
            axis.title = element_text(face='bold'),
            plot.title = element_text(hjust=0.5, face='bold'))
    plots$mt <- p.mt
    ggsave(paste0(resDir,'all_samples_MT.png'), plot=p.mt, dpi=400, width=6, height=6)
  }
  
  ### (5) Density plots (all samples, all metrics)
  if (density) {
    p.density1 <- ggplot(qc_df, aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident, group=orig.ident)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() + 
      scale_x_log10() + 
      scale_fill_manual(values = custom_cols) +
      scale_color_manual(values = custom_cols) +
      geom_vline(xintercept = 350,color="red",linetype="dotted") + # Include this if plotting threshold line
      theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = 'none') +
      ggtitle("nFeature Distribution")
    p.density2 <- ggplot(qc_df, aes(x=nCount_RNA, color=orig.ident, fill=orig.ident, group=orig.ident)) + 
      geom_density(alpha = 0.2) + 
      theme_classic() + 
      scale_x_log10() + 
      scale_fill_manual(values = custom_cols) +
      scale_color_manual(values = custom_cols) +
      geom_vline(xintercept = 1000,color="red",linetype="dotted") + # Include this if plotting threshold line
      theme(plot.title = element_text(hjust=0.5, face="bold"), legend.position = 'none') +
      ggtitle("nCount Distribution")
    p.density3 <- ggplot(qc_df, aes(x=percent.mt, color=orig.ident, fill=orig.ident)) + 
      geom_density(alpha = 0.2) +
      theme_classic() +
      scale_x_log10(labels = label_comma()) +
      scale_fill_manual(values = custom_cols) +
      scale_color_manual(values = custom_cols) +
      geom_vline(xintercept = 15,color="red",linetype="dotted") +
      theme(plot.title = element_text(hjust=0.5, face="bold")) +
      ggtitle("MT % Expression") +
      guides(color = guide_legend(title = "Sample"),
             fill = guide_legend(title = "Sample"))
    density_plots <- p.density1 + p.density2 + p.density3
    plots$density <- density_plots
    ggsave(paste0(resDir,'all_samples_density_plots.png'), plot=density_plots, dpi=400, width=15, height=6)
  }
    return(plots)
}

# ==============================================================================
# Filter seurat object
# ==============================================================================

filter_out_cells <- function(seu.obj, feature_min=NULL, count_min=NULL, mt_threshold=NULL) {
  
  print('Filtering out low quality cells...')
  
  prefilter_count = length(Cells(seu.obj))
  
  seu.obj@meta.data$keep <- with(seu.obj@meta.data, ifelse(nFeature_RNA > feature_min & nCount_RNA > count_min & percent.mt < mt_threshold, TRUE, FALSE))
  seu.obj <- subset(seu.obj, subset = keep == TRUE)

  postfilter_count = length(Cells(seu.obj))
  
  total = prefilter_count - postfilter_count
  print(glue::glue("Cells removed from {sample_id} : {total}"))
  
  return(seu.obj)
}

# ==============================================================================
# Grab CPMs
# ==============================================================================

# TODO: insert function to extract cpm matrix from either (1) seurat obj or (2) cellranger output



# ==============================================================================
# Prettier UMAPs & VolcanoPlots
# ==============================================================================

############ formatUMAP ############
# Function adapted from Ammons repo https://github.com/dyammons/canine_osteosarcoma_atlas/blob/main/analysisCode/customFunctions.R

formatUMAP <- function(plot = NULL, smallAxes = F) {
  
  plot <- plot + labs(x = "UMAP1", y = "UMAP2") +
    theme(axis.text = element_blank(), 
          axis.ticks = element_blank(),
          axis.title = element_text(size= 20),
          plot.title = element_blank(),
          title = element_text(size= 20),
          axis.line = element_blank(),
          panel.border = element_rect(color = "black",
                                      fill = NA,
                                      size = 2)
    )
  
  if(smallAxes){
    
    axes <- ggplot(data.frame(x = 0:1, y = 0:1), aes(x, y)) +
      geom_blank() +
      coord_fixed() +
      labs(x = "UMAP1", y = "UMAP2") +
      theme(
        axis.line = element_line(
          colour = "black",
          arrow = arrow(angle = 30,
                        length = unit(0.1, "inches"),
                        ends = "last",
                        type = "closed")
        ),
        axis.title.y = element_text(colour = "black", size = 20),
        axis.title.x = element_text(colour = "black", size = 20, margin = margin(t = 8)),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    plot <- plot + theme(axis.title = element_blank(),
                         panel.border = element_blank(),
                         plot.margin = unit(c(-7, -7, -7, -7), "pt"))
    
    plot <- plot + inset_element(axes,
                                 left = 0,
                                 bottom = 0,
                                 right = 0.25,
                                 top = 0.25,
                                 align_to = "full")
  }
  
  return(plot)
}


############ prettyVolc ############                                          
prettyVolc <- function(plot = NULL, rightLab = NULL, leftLab = NULL, rightCol = "red", leftCol = "blue", arrowz = T
){
  
  p <- plot + scale_x_symmetric(mid = 0) + theme(legend.position = c(0.10, 0.9),
                                                 legend.background = element_blank(),
                                                 legend.key = element_blank(),
                                                 axis.title = element_text(size = 16),
                                                 axis.text = element_text(size = 12),
                                                 panel.grid.major = element_blank(),
                                                 panel.grid.minor = element_blank(),
                                                 panel.border = element_blank(),
                                                 panel.background = element_blank(),
                                                 axis.line = element_line(color="black"),
                                                 plot.title = element_blank()
  ) + 
    {if(arrowz){
      annotate("segment", x = 0.58*1.5, 
               y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
               xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[1], 
               yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
               lineend = "round", linejoin = "bevel", linetype ="solid", colour = rightCol,
               size = 1, arrow = arrow(length = unit(0.1, "inches"))
      ) 
    }} +
    {if(arrowz){
      annotate("segment", x = -0.58*1.5, 
               y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
               xend = c(max(abs(plot$data$log2FoldChange)),-max(abs(plot$data$log2FoldChange)))[2],
               yend = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.06, 
               lineend = "round", linejoin = "bevel", linetype ="solid", colour = leftCol,
               size = 1, arrow = arrow(length = unit(0.1, "inches"))
      )
    }} + 
    {if(!is.null(rightLab)){
      annotate(geom = "text", x = (max(abs(plot$data$log2FoldChange))-0.58*1.5)/2+0.58*1.5, 
               y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
               label = rightLab,
               hjust = 0.5,
               size = 5)
    }} + 
    {if(!is.null(leftLab)){
      annotate(geom = "text", x = -(max(abs(plot$data$log2FoldChange))-0.58*1.5)/2-0.58*1.5, 
               y = ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range[2]*1.09,
               label = leftLab,
               hjust = 0.5,
               size = 5)
    }} 
  
  return(p)
}






# Testing this function from the Ammons repo
########### prettyFeats ############
prettyFeats <- function(seu.obj = NULL, nrow = 3, ncol = NULL, features = "", color = "black", order = FALSE, titles = NULL, noLegend = F, bottomLeg = F, min.cutoff = NA, pt.size = NULL, title.size = 18, legJust = "bottom",showAxis = F,smallAxis=F, legInLine = F, returnPlots = F
) {
  
  DefaultAssay(seu.obj) <- "RNA"
  features <- features[features %in% c(unlist(seu.obj@assays$RNA@counts@Dimnames[1]),unlist(colnames(seu.obj@meta.data)))]
  
  if(is.null(ncol)){
    ncol = ceiling(sqrt(length(features)))
  }
  
  if(is.null(titles)){
    titles <- features #- add if statement
  }
  
  
  #strip the plots of axis and modify titles and legend -- store as large list
  plots <- Map(function(x,y,z) FeaturePlot(seu.obj,features = x, pt.size = pt.size, order = order, min.cutoff = min.cutoff) + labs(x = "UMAP1", y = "UMAP2") +
                 theme(axis.text= element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title = element_blank(), 
                       axis.line = element_blank(),
                       title = element_text(size= title.size, colour = y),
                       legend.position = "none"
                 ) + 
                 scale_color_gradient(breaks = pretty_breaks(n = 3), limits = c(NA, NA), low = "lightgrey", high = "darkblue") + 
                 ggtitle(z), x = features, y = color, z = titles) 
  
  
  asses <- ggplot() + labs(x = "UMAP1", y = "UMAP2") + 
    theme(axis.line = element_line(colour = "black", 
                                   arrow = arrow(angle = 30, length = unit(0.1, "inches"),
                                                 ends = "last", type = "closed"),
    ),
    axis.title.y = element_text(colour = "black", size = 20),
    axis.title.x = element_text(colour = "black", size = 20),
    panel.border = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
    )
  if(!noLegend){
    leg <- FeaturePlot(seu.obj,features = features[1], pt.size = 0.1) + 
      theme(legend.position = 'bottom',
            legend.direction = 'vertical',
            legend.justification = "center",
            
      ) + 
      scale_color_gradient(breaks = pretty_breaks(n = 1), labels = c("low", "high"), limits = c(0,1), low = "lightgrey", high = "darkblue") + 
      guides(color = guide_colourbar(barwidth = 1)) 
    
    if(bottomLeg){
      leg <- leg + theme(legend.direction = 'horizontal') + guides(color = guide_colourbar(barwidth = 8))
    }
    
    legg <- get_legend(leg)
  }
  
  #nrow <- ceiling(length(plots)/ncol) - add if statement
  patch <- area()
  
  counter=0
  for (i in 1:nrow) {
    for (x in 1:ncol) {
      counter = counter+1
      if (counter <= length(plots)) {
        patch <- append(patch, area(t = i, l = x, b = i, r = x))
      }
    }
  }
  
  if(!noLegend){
    if(!bottomLeg & !legInLine){
      legPos <- ifelse(legJust == "bottom",ceiling(length(features)/ncol),1)
      patch <- append(patch, area(t = legPos, l = ncol+1, b = legPos, r = ncol+1))
    }else if(legInLine){
      patch <- append(patch, area(t = nrow, l = ncol, b = nrow, r = ncol))
    }else{
      patch <- append(patch, area(t = ceiling(length(features)/ncol)+1, l = ncol, b = ceiling(length(features)/ncol)+1, r = ncol))
    }
  }else{
    if(smallAxis){
      patch <- append(patch, area(t = nrow, l = 1, b = nrow, r = 1))
    }else{
      patch <- append(patch, area(t = 1, l = 1, b = nrow, r = ncol))
    }
  }
  
  
  if(returnPlots){
    return(plots)
  }else{
    p <- Reduce( `+`, plots ) +  {if(noLegend & showAxis){asses}} +
      {if(!noLegend){legg}} + plot_layout(guides = "collect") +
      {if(!noLegend & bottomLeg){plot_layout(design = patch, heights = c(rep.int(1, nrow),0.2))}else if(!noLegend & !bottomLeg){plot_layout(design = patch, widths = c(rep.int(1, ncol),0.2))}else{plot_layout(design = patch, widths = rep.int(1, ncol))}}
    
    return(p)
  }
}









