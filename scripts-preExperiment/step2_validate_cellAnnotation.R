##################################################
## Project: SP pre-experiments 
## Script purpose: validate the cell type annotation
## Date: 2020-08-27
## Author: Yuting Liu
##################################################

## Section: set env
##################################################
setwd("/lustre/user/liclab/liuyt/SP/pre-scRNA")
library(Seurat)
library(dplyr)
library("viridis") 
library(DataCombine)


## Section: load the seurat object
##################################################
load('data/preB.seurat.pca30.umap.tsne.res.1.5.withAnnotation.RData')
#load('data/PreA.seurat.pca20.res.1.FindAllMarkers.genes.RData')

## Section: defining stacked vlnplot functions
##################################################
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.15, 0, -0.15, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(0.5)), 
          plot.margin = plot.margin ) +
    scale_fill_manual(values = c('Mic'=rgb(225,107,167, maxColorValue = 255),'Per'=rgb(196,115,171,maxColorValue = 255),
                                 'End'=rgb(144,136,192,maxColorValue = 255),'OPC'=rgb(91,153,210, maxColorValue = 255),
                                 'InCGE'=rgb(73,178,73,maxColorValue = 255),'InMGE'=rgb(47,178,81,maxColorValue = 255),
                                 'ExDp2'=rgb(168,126,182,maxColorValue = 255),'ExDp1'=rgb(174,159,51,maxColorValue = 255),
                                 'ExM-U'=rgb(135,170,63, maxColorValue = 255),'ExM'=rgb(226,131,37,maxColorValue = 255),
                                 'ExN'=rgb(240,117,109,maxColorValue = 255),'IP'=rgb(205,147,41,maxColorValue = 255),
                                 'PgG2M'=rgb(44,173,227,maxColorValue = 255),'PgS'=rgb(34,185,174,maxColorValue = 255),
                                 'oRG'=rgb(49,181,132,maxColorValue = 255),'vRG'=rgb(0,185,211, maxColorValue = 255))) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.15, 0, -0.15, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90), axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
 # ymaxs<- purrr::map_dbl(plot_list, extract_max)
#  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
 #                           scale_y_continuous(breaks = c(y)) + 
  #                          expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

## Section: vlnplot
##################################################
levels(pre) <- c('InCGE','InMGE','ExDp2','ExDp1','ExM-U','ExM','ExN','IP','PgG2M','oRG','OPC','End','Per','Mic')
gn.ls <- c('PTPRC','P2RY12','FN1','OLIG2','HOPX','PAX6','EOMES','PPP1R17',
           'NRP1','NEUROD6','SATB2','LPL','TBR1','GAD1')

pdf('results/QC_and_mergeIndividual/preB_res.1.5_MakerGene_stacked_Vlnplot.pdf', width = 6, height = 14)
StackedVlnPlot(pre, features = gn.ls)
dev.off()

## Section: heatmap
##################################################
pre.markers <- FindAllMarkers(pre, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- pre.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


g <- DoHeatmap(pre, features = top10$gene, angle = 90, group.bar = T, raster = F, size = 2 ) + scale_fill_viridis(option = "D") 
pdf('results/QC_and_mergeIndividual/preB_res.1.5_top5MakerGene_heatmap.pdf', width = 6, height = 6.4)
g
dev.off()

## Section: featureplot
##################################################
pdf('results/QC_and_mergeIndividual/preB_res.1.5_CanonicalMakerGene_FeaturePlot.pdf', width =15, height = 10 )
FeaturePlot(pre, features = c('PTPRC','P2RY12','FN1',
                              'HOPX','PAX6','OLIG2','EOMES','PPP1R17',
                              'GAD1','LHX6','NRP1','SATB2','NEUROD6',
                              'LPL','TBR1', 'ST18','EPHA7','NR4A2'),
            pt.size = .2, label.size = 4, reduction = 'tsne', ncol=5)
dev.off()
