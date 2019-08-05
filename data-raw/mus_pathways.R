# Mus musculus REACTOME pathways

library(Matrix)

zz <- download.file("https://reactome.org/download/current/Ensembl2Reactome.txt",
                    destfile = "pathways.txt")

pathways <- data.table::fread("pathways.txt", header = FALSE, stringsAsFactors = TRUE)

pathways <- subset(pathways, V6 == "Mus musculus")

# Store as a sparse matrix of 1s and 0s.
pw.mat <- reshape2::acast(pathways, V1 ~ V4, length)
pw.mat <- Matrix(pw.mat)

# Convert ENSEMBL IDs to gene symbols.
ensembl.ids <- rownames(pw.mat)
ensembl <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
bm <- biomaRt::getBM(filters = "ensembl_gene_id",
                     attributes = c("mgi_symbol", "ensembl_gene_id"),
                     values = ensembl.ids,
                     mart = ensembl)
rownames(pw.mat) <- bm$mgi_symbol[match(ensembl.ids, bm$ensembl_gene_id)]

# Deal with duplicate gene symbols.
dups <- names(which(table(rownames(pw.mat)) > 1))
for (dup in dups) {
  idx <- which(rownames(pw.mat) == dup)
  pw.mat[idx[1], ] <- pmin(colSums(pw.mat[idx, ]), 1)
  pw.mat <- pw.mat[-idx[-1], ]
}

saveRDS(pw.mat, "./data/mus_pathways.rds")

zz <- file.remove("pathways.txt")
