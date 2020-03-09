library(SingleCellExperiment)
library(Matrix)
library(tidyverse)

load("./data/10x/Sce_Dataset2.RData")

pbmc <- Matrix(counts(sce))
# Remove genes with all zero counts.
pbmc <- pbmc[rowSums(pbmc) > 0, ]

saveRDS(pbmc, "./data/10x/pbmc.rds")
