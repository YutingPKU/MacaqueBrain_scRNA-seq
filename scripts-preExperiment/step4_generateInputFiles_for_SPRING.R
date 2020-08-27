##################################################
## Project: pre-experiment
## Script purpose: generate the matrix, gene list, cell clustering input files for SPRING
## Date: 2020-08-27
## Author: Yuting Liu
##################################################

## Section: set env
##################################################
setwd("/lustre/user/liclab/liuyt/SP/pre-scRNA/")
library(Matrix)

## Section: load the seurat object
##################################################
load('data/preB.seurat.pca30.umap.tsne.res.1.5.withAnnotation.RData')

## Section: get gene expression matrix, gene list and cluster info
##################################################
mt <- pre@assays$RNA@counts
writeMM(obj = mt, file="data/SRPING/preB_matrix.mtx")

# save genes and cells names
write(x = rownames(mt), file = "data/SRPING/preB_genes.tsv")
write(x = colnames(mt), file = "data/SRPING/preB_barcodes.tsv")

id <- as.character(pre@active.ident)
id <- c('type', id)
cl <- as.character(pre@meta.data$RNA_snn_res.1.5)
cl <- c('cluster',cl)
df <- data.frame(t(data.frame(cbind(cl,id))))
write.table(df, file = 'data/SRPING/preB_clustermeta.csv', col.names = F, row.names = F, quote = F, sep = ",")
