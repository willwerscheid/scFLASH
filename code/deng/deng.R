# remotes::install_github("kkdey/singleCellRNASeqMouseDeng2014")

library(flashier)
library(singleCellRNASeqMouseDeng2014)

# Load data. ------------------------------------------------------------------

counts <- exprs(Deng2014MouseESC)
meta_data <- pData(Deng2014MouseESC)
gene_names <- rownames(counts)

preprocess <- function(dat, min.nzcts = 10) {
  size.factors <- colSums(dat)
  size.factors <- size.factors / mean(size.factors)
  gene_cts <- rowSums(dat > 0)
  dat <- dat[gene_cts >= min.nzcts, ]

  lunpc <- max(1 / min(size.factors) - 1 / max(size.factors), 1)
  fl.dat <- log1p(t(t(dat) / size.factors) / lunpc)

  return(list(
    dat = dat,
    fl.dat = fl.dat,
    size.factors = size.factors,
    excluded.genes = gene_cts < min.nzcts)
  )
}
Deng <- preprocess(counts)

# Run flashier. ---------------------------------------------------------------

# fl <- flash.init(Deng$fl.dat, var.type = 1) %>%
#   flash.add.greedy(
#     Kmax = 20,
#     ebnm.fn = ebnm::ebnm_point_normal
#   ) %>%
#   flash.backfit(verbose = 3) %>%
#   flash.nullcheck()

snmf.fl <- flash.init(Deng$fl.dat, var.type = 1) %>%
  flash.add.greedy(
    Kmax = 25,
    ebnm.fn = c(ebnm::ebnm_normal_scale_mixture, ebnm::ebnm_point_exponential),
    init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
  ) %>%
  flash.backfit(verbose = 3)

saveRDS(snmf.fl, "./data/deng_fl.rds")
