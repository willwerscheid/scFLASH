# Droplet dataset from Montoro et al.

library(Matrix)

zz <- download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE103354&format=file&file=GSE103354%5FTrachea%5Fdroplet%5FUMIcounts%2Etxt%2Egz",
                    destfile = "droplet.txt.gz")

droplet <- data.table::fread("droplet.txt.gz")

# fread puts the row names (genes) in the first column.
genes   <- droplet$V1
droplet <- droplet[, -1]

# Store as a sparse matrix.
droplet <- Matrix(as.matrix(droplet))
rownames(droplet) <- genes

saveRDS(droplet, "./data/droplet.rds")

file.remove("droplet.txt.gz")
