##################################################
## Project: SP pre-experiments 
## Script purpose: comparing cell type between individuals by cheking the overlap of marker genes
## Date: 2020-08-27
## Author: Yuting Liu
##################################################

## Section: set env
##################################################
setwd("/lustre/user/liclab/liuyt/SP/pre-scRNA")
library(Seurat)



## Section: load the seurat object
##################################################
load('data/PreA.seurat.pca20.res.1.withAnnotation.FindAllMarkers.genes.RData')
a.ls <- split(rownames(pre.markers), f = pre.markers$cluster)
load('data/preB.seurat.pca30.res.1.5.WithAnnotation.FindAllMarkers.genes.RData')
b.ls <- split(rownames(pre.markers), f = pre.markers$cluster)
names(a.ls)
names(b.ls)
b.ls <- b.ls[names(a.ls)]

## Section: fish'er exact test to check the overlap btw two sets
##################################################

getP <- function(v1,v2){
  n1 <- length(v1)
  n2 <- length(v2)
  n <- length(intersect(v1,v2))
  Convictions <- matrix(c(n1,n,(15000-n1),(n2-n)), nrow = 2,
                        dimnames =
                          list(c("Dizygotic", "Monozygotic"),
                               c("Convicted", "Not convicted")))
  re <- fisher.test(Convictions)
  return(-log10(re$p.value))
}

df <- data.frame(matrix(0, nrow = 11, ncol = 11))
rownames(df) <- paste0("A_", names(a.ls))
colnames(df) <- paste0('B_', names(b.ls))

for (i in 1:length(a.ls)) {
  for (j in 1:length(b.ls)) {
    p <- getP(unlist(a.ls[i]), unlist(b.ls[j]))
    df[i,j] <- p
  }
}

pdf('results/QC_and_mergeIndividual/preA.preB.cellCluster.markergenes.overlap.significance.pdf',width = 3.5, height = 3)
re <- pheatmap(df, cluster_cols = F, cluster_rows = F, color = colorRampPalette(c('white','red'))(40), border_color = NA)
plot.new()
re <- pheatmap(df, cluster_cols = F, cluster_rows = F, color = colorRampPalette(c('white','red'))(40))
dev.off()
