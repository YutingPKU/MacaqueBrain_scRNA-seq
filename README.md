# MacaqueBrain_scRNA-seq
## Aug 20,2020
pre-experiment for two individuals: 11002B(A) and 11002C(B)

`scripts-preExperiment`

1. Aligning and counting reads using CellRanger per individuals.
2. QC and cell type annotation using Seurat V3 per individuals.
3. Comparing cell composition between two individuals:
  - overlap of cell type marker genes;
  - comparing cell trajectory using SPRING;
