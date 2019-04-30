library(Matrix)
set.seed(666)

# The drop-seq dataset in Montoro et al. can be downloaded here:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354
trachea <- read.table("./data/GSE103354_Trachea_droplet_UMIcounts.txt")
trachea <- as.matrix(trachea)
trachea <- Matrix(trachea)

# I only retain the 158 tuft cells. The idea is that the population should
#   be homogeneous enough for a smallish number of factors to be able to
#   adequately describe the dataset.
cell.type <- sapply(sapply(colnames(trachea), strsplit, '_'), `[[`, 3)
tuft <- trachea[, cell.type == "Tuft"]

# Remove genes that have non-zero counts in 3 or fewer cells.
gene.cts <- rowSums(tuft > 0)
tuft <- tuft[gene.cts > 3, ]

# I include 500 highly variable genes.
gene.means <- apply(tuft, 1, mean)
gene.vars <- apply(tuft, 1, var)
gene.disps <- gene.vars / gene.means
top.genes <- order(gene.disps, decreasing = TRUE)[1:500]

# The selection process in biased in favor of highly expressed genes. To
#   make the overall sparsity more similar to the sparsity of the original
#   dataset, I also include 1500 randomly selected genes.
remaining.genes <- setdiff(1:nrow(tuft), top.genes)
rand.genes <- sample(remaining.genes, 1500)

tuft <- tuft[c(top.genes, rand.genes), ]

# The resulting dataset is 158 x 2000, with 77% of counts equal to zero.
saveRDS(tuft, "./data/tuft.rds")
