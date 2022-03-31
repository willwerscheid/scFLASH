# remotes::install_github("kkdey/singleCellRNASeqMouseDeng2014")

library(flashier)
library(Matrix)

# Load data. ------------------------------------------------------------------

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
    pc = lunpc,
    excluded.genes = gene_cts < min.nzcts)
  )
}

pbmc <- readRDS("../flashier-chapter/data/pbmc.rds")
pbmc <- pbmc[, colSums(pbmc) < 15000]
pbmc <- preprocess(pbmc)


# Run flashier. ---------------------------------------------------------------

n <- ncol(pbmc$fl.dat)
min.sd <- sd(log1p(rpois(1e7, 1/n) / pbmc$pc / median(pbmc$size.factors)))
disp.min.sd <- function(new, old, k) {
  return(sqrt(min(1 / ff.tau(new))))
}

K <- 20
t0 <- Sys.time()
snmf.fl <- flash.init(pbmc$fl.dat, S = min.sd, var.type = 1) %>%
  flash.init.factors(
    list(matrix(rowMeans(pbmc$fl.dat), ncol = 1), matrix(1, nrow = ncol(pbmc$fl.dat), ncol = 1)),
    ebnm.fn = c(ebnm::ebnm_point_normal, ebnm::ebnm_point_exponential)
  ) %>%
  flash.fix.factors(kset = 1, mode = 2) %>%
  flash.backfit() %>%
  flash.add.greedy(
    Kmax = K - 1,
    ebnm.fn = c(ebnm::ebnm_point_normal, ebnm::ebnm_point_exponential),
    init.fn = function(f) init.fn.default(f, dim.signs = c(0, 1))
  ) %>%
  flash.set.verbose(
    disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
    colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
    colwidths = c(18, 14, 14, 14)
  ) %>%
  flash.backfit(verbose = 3)
t1 <- Sys.time()
snmf.t <- t1 - t0
saveRDS(list(fl = snmf.fl, t = snmf.t), "./output/pbmc/smnf.rds")


t0 <- Sys.time()
greedy.fl <- flash.init(pbmc$fl.dat, S = min.sd, var.type = 1) %>%
  flash.init.factors(
    list(matrix(rowMeans(pbmc$fl.dat), ncol = 1), matrix(1, nrow = ncol(pbmc$fl.dat), ncol = 1)),
    ebnm.fn = ebnm::ebnm_point_exponential
  ) %>%
  flash.fix.factors(kset = 1, mode = 2) %>%
  flash.backfit() %>%
  flash.add.greedy(
    Kmax = K - 1,
    ebnm.fn = ebnm::ebnm_point_exponential,
    init.fn = function(f) init.fn.default(f, dim.signs = c(1, 1))
  )
t1 <- Sys.time()
greedy.t <- t1 - t0
saveRDS(list(fl = greedy.fl, t = greedy.t), "./output/pbmc/greedy.rds")


t0 <- Sys.time()
bf.fl <- greedy.fl %>%
  flash.set.verbose(
    disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
    colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
    colwidths = c(18, 14, 14, 14)
  ) %>%
  flash.backfit(verbose = 3)
t1 <- Sys.time()
bf.t <- t1 - t0
saveRDS(list(fl = bf.fl, t = bf.t), "./output/pbmc/bf.rds")


t0 <- Sys.time()
nnmf.res <- NNLM::nnmf(
  as.matrix(pbmc$fl.dat),
  init = list(H0 = matrix(1, nrow = 1, ncol = ncol(pbmc$fl.dat))),
  k = K - 1,
  verbose = 0
)
nnlm.fl <- flash.init(pbmc$fl.dat, S = min.sd, var.type = 1) %>%
  flash.set.verbose(0) %>%
  flash.init.factors(
    list(nnmf.res$W[, c(K, 1:(K-1))], t(nnmf.res$H[c(K, 1:(K-1)), ])),
    ebnm.fn = ebnm::ebnm_point_exponential
  ) %>%
  flash.fix.factors(kset = 1, mode = 2)
t1 <- Sys.time()
nnlm.t <- t1 - t0
saveRDS(list(fl = nnlm.fl, t = nnlm.t), "./output/pbmc/nnlm.rds")


t0 <- Sys.time()
nnlmbf.fl <- nnlm.fl %>%
  flash.set.verbose(
    disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
    colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
    colwidths = c(18, 14, 14, 14)
  ) %>%
  flash.backfit(verbose = 3)
t1 <- Sys.time()
nnlmbf.t <- t1 - t0
saveRDS(list(fl = nnlmbf.fl, t = nnlmbf.t), "./output/pbmc/nnlmbf.rds")


t0 <- Sys.time()
nzpe.fl <- flash.init(pbmc$fl.dat, S = min.sd, var.type = 1) %>%
  flash.set.verbose(0) %>%
  flash.init.factors(
    list(nnmf.res$W[, c(K, 1:(K-1))], t(nnmf.res$H[c(K, 1:(K-1)), ])),
    ebnm.fn = c(
      ebnm::ebnm_point_exponential,
      as.ebnm.fn(prior_family = "point_exponential", mode = "estimate")
    )
  ) %>%
  flash.fix.factors(kset = 1, mode = 2) %>%
  flash.set.verbose(
    disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
    colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
    colwidths = c(18, 14, 14, 14)
  ) %>%
  flash.backfit(verbose = 3)
t1 <- Sys.time()
nzpe.t <- t1 - t0
saveRDS(list(fl = nzpe.fl, t = nzpe.t), "./output/pbmc/nzpe.rds")

t0 <- Sys.time()
shifts <- sapply(nzpe.fl$F.ghat[-1], function(g) g$shift[1])
LL <- nzpe.fl$L.pm
LL[, 1] <- nzpe.fl$L.pm %*% c(1, shifts)
FF <- nzpe.fl$F.pm
FF[, -1] <- FF[, -1] - rep(shifts, each = nrow(FF))
nzpe.fl2 <- flash.init(pbmc$fl.dat, S = min.sd, var.type = 1) %>%
  flash.set.verbose(0) %>%
  flash.init.factors(
    list(LL[, 1, drop = FALSE], FF[, 1, drop = FALSE]),
    ebnm.fn = c(ebnm::ebnm_point_laplace, ebnm::ebnm_point_exponential)
  ) %>%
  flash.fix.factors(kset = 1, mode = 2) %>%
  flash.init.factors(
    list(LL[, -1], FF[, -1]),
    ebnm.fn = ebnm::ebnm_point_exponential
  ) %>%
  flash.set.verbose(
    disp.fns = c(display.elbo, display.elbo.diff, display.F.max.chg, disp.min.sd),
    colnames = c("ELBO", "Diff", "Max Chg (F)", "Min SD"),
    colwidths = c(18, 14, 14, 14)
  ) %>%
  flash.backfit(verbose = 3)
t1 <- Sys.time()
nzpe.t2 <- t1 - t0
saveRDS(list(fl = nzpe.fl2, t = nzpe.t2), "./output/pbmc/nzpe2.rds")
