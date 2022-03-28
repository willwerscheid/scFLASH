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

disp.min.sd <- function(new, old, k) {
  return(sqrt(min(1 / ff.tau(new))))
}

do.fit <- function(K) {
  fl <- flash.init(
      Deng$fl.dat,
      S = min(Deng$fl.dat[Deng$fl.dat > 0]) / sqrt(ncol(Deng$fl.dat)),
      var.type = 1
    ) %>%
    flash.add.greedy(
      Kmax = K,
      ebnm.fn = c(ebnm::ebnm_normal_scale_mixture, ebnm::ebnm_point_exponential),
      init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
    ) %>%
    flash.set.verbose(
      disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
      colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
      colwidths = c(18, 14, 14, 14)
    ) %>%
    flash.backfit(verbose = 3)

  saveRDS(fl, paste0("./output/deng/deng_fl", K, ".rds"))
  return(fl)
}

fl6  <- do.fit(6)
fl10 <- do.fit(10)
fl25 <- do.fit(25)
