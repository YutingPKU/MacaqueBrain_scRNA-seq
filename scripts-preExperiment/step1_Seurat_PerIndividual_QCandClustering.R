##################################################
## Project: SP pre-experiment 
## Script purpose: QC and clustering on individual data
## Date: 2020-08-26
## Author: Yuting Liu
##################################################

## Section: set env
##################################################
setwd("/lustre/user/liclab/liuyt/SP/pre-scRNA")
library(Seurat)
library(dplyr)

## Section: setup the seurat object
##################################################
pre <- Read10X(data.dir = 'cellranger-res/Pre-B/outs/filtered_feature_bc_matrix/')
pre <- CreateSeuratObject(counts = pre, project = '11002C', min.cells = 3)
pre

## Section: QC and selecting cells
##################################################
pre[["percent.mt"]] <- PercentageFeatureSet(pre, 
                        features = c("MTARC2","MTFR1L","MTERF1","MTFR2","MTRF1L","MTRES1",
                                     "MTO1","MTCH1","MTFMT","MTFR1","MTERF3","MTERF2","MTPAP",
                                     "MTERF4","MTCH2",'MTIF2',"MTG2","MTIF3","MTRF1","MTCL1"))
VlnPlot(pre, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pre, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#pre <- subset(pre, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA < 15000) #preA
pre <- subset(pre, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 11000) #preB
pre

## Section: normalizing
##################################################
pre <- NormalizeData(pre, normalization.method = "LogNormalize", scale.factor = 10000)

## Section: find variable genes
##################################################
pre <- FindVariableFeatures(pre, selection.method = "mvp", mean.cutoff = c(0.01, 5), 
                            dispersion.cutoff  = c(1,Inf))

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pre), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pre)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Section: scaling 
##################################################
all.genes <- rownames(pre)
pre <- ScaleData(pre, features = all.genes)

## Section: dimensional reduction
##################################################
pre <- RunPCA(pre, features = VariableFeatures(object = pre))
VizDimLoadings(pre, dims = 1:2, reduction = "pca")
DimPlot(pre, reduction = "pca")
DimHeatmap(pre, dims = 1:20, cells = 500, balanced = TRUE)

pre <- JackStraw(pre, num.replicate = 100, dims = 40)
pre <- ScoreJackStraw(pre, dims = 1:40)
JackStrawPlot(pre, dims = 1:40)
ElbowPlot(pre, ndims = 40)

## Section: cluster
##################################################
pre <- FindNeighbors(pre, dims = 1:30)
pre <- FindClusters(pre, resolution = c(0.5, 0.6, 0.8, 1, 1.2, 1.5,2) )

pre <- RunUMAP(pre, dims = 1:30)
pre <- RunTSNE(pre, dims = 1:30)

Idents(object = pre) <- 'RNA_snn_res.1'
DimPlot(pre, reduction = "umap", label = T)
DimPlot(pre, reduction = "tsne", label = T)

saveRDS(pre, file = "data/preB.seurat.pca30.umap.tsne.rds")

Idents(object = pre) <- 'RNA_snn_res.1'

pdf('results/QC_and_mergeIndividual/preB_tSNE_res.1.5.Cluster.pdf', width = 4.5, height = 4)
DimPlot(pre, reduction = "tsne", label = T, label.size = 4, pt.size = 0.4)+NoLegend()
dev.off()

pdf('results/QC_and_mergeIndividual/preB_UMAP_res.1.5.Cluster.pdf', width = 4.5, height = 4)
DimPlot(pre, reduction = "umap", label = T, label.size = 4, pt.size = 0.4)+NoLegend()
dev.off()

## Section: find markers
##################################################
pre <- readRDS('data/preB.seurat.pca30.umap.tsne.rds')
Idents(object = pre) <- 'RNA_snn_res.1.5'

pre.markers <- FindAllMarkers(pre, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pre.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

save(pre.markers, file = 'data/preB.seurat.pca30.res.1.5.FindAllMarkers.genes.RData')



id = c('CENPF','TOP2A','NUSAP1','ASPM','PRC1','CENPE') #PgG2M
id = c('BCAN','PTN','APOD','EPN2') #OPC
id = c('PPP1R17','EOMES','CORO1C') # IP
id = c('PTN','VIM','SFRP1','SLC1A3','HOPX') #oRG
id = c('CALB2') #InCGE
id = c('SST','MAF','PDZRN4') #InMGE

VlnPlot(pre, features = id, ncol = 3)
FeaturePlot(pre, features =id)
top10 <- pre.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pre, features = top10$gene) + NoLegend()

## Section: add cell type annotation
##################################################
#new.cluster.ids <- c('ExM-U','ExN','ExM','InMGE','InCGE','ExDp1','oRG','ExM','ExM','ExM-U','IP','ExDp2','PgG2M','InMGE','OPC')
new.cluster.ids <- c('InMGE','InCGE','ExM-U','InCGE','ExM-U','ExN','ExM','IP','PgG2M','oRG','ExDp1',
                     'ExM-U','InMGE','InCGE','OPC','ExN','ExN','ExDp2','ExN','End','Per','Mic')
names(new.cluster.ids) <- levels(pre)
pre <- RenameIdents(pre, new.cluster.ids)
save(pre, file = 'data/preB.seurat.pca20.umap.tsne.res.1.5.withAnnotation.RData')

pdf('results/QC_and_mergeIndividual/preB_tSNE_res.1.5.CellType.pdf', width = 4.5, height = 4)
DimPlot(pre, reduction = "tsne", label = TRUE , pt.size = 0.4, label.size = 4, 
        cols= c('ExM-U' = rgb(61,104,233, maxColorValue = 255),'ExN' = rgb(0,253,254, maxColorValue = 255),
                'ExM' = rgb(120,202,122, maxColorValue = 255), 'InMGE' = rgb(254,214,2,maxColorValue = 255),
                'InCGE' = rgb(253,162,4, maxColorValue = 255), 'ExDp1' = rgb(216,110,210, maxColorValue = 255),
                'oRG' = rgb(169,127,254,maxColorValue = 255), 'IP' = rgb(82,25,137,maxColorValue = 255),
                'ExDp2' = rgb(121,248,5, maxColorValue = 255),'PgG2M' = rgb(27,138,31, maxColorValue = 255),
                'OPC' = rgb(253,69,0,maxColorValue = 255),'Mic' = rgb(156,47,239,maxColorValue = 255),
                'End' = rgb(92,69,137, maxColorValue = 255),'Per' = rgb(138,0,0,maxColorValue = 255)))+NoLegend()
dev.off()



pdf('results/QC_and_mergeIndividual/preB_UMAP_res.1.5.CellType.pdf', width = 4.5, height = 4)
DimPlot(pre, reduction = "umap", label = TRUE , pt.size = 0.4, label.size = 4, 
        cols= c('ExM-U' = rgb(61,104,233, maxColorValue = 255),'ExN' = rgb(0,253,254, maxColorValue = 255),
                'ExM' = rgb(120,202,122, maxColorValue = 255), 'InMGE' = rgb(254,214,2,maxColorValue = 255),
                'InCGE' = rgb(253,162,4, maxColorValue = 255), 'ExDp1' = rgb(216,110,210, maxColorValue = 255),
                'oRG' = rgb(169,127,254,maxColorValue = 255), 'IP' = rgb(82,25,137,maxColorValue = 255),
                'ExDp2' = rgb(121,248,5, maxColorValue = 255),'PgG2M' = rgb(27,138,31, maxColorValue = 255),
                'OPC' = rgb(253,69,0,maxColorValue = 255),'Mic' = rgb(156,47,239,maxColorValue = 255),
                'End' = rgb(92,69,137, maxColorValue = 255),'Per' = rgb(138,0,0,maxColorValue = 255)))+NoLegend()

dev.off()

