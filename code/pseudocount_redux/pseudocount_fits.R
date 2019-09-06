source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")

do.fit <- function(size.factors, pseudocount, greedy.Kmax = 20) {
  cat("FITTING PSEUDOCOUNT:", as.character(pseudocount), "\n")

  processed <- do.call(preprocess,
                       c(list(data = droplet,
                              min.nzcts = 10,
                              max.libsize = 30000,
                              size.factors = pseudocount * size.factors,
                              prescale = "cells")))

  fl.g <- flashier(processed$data, var.type = 1, greedy.Kmax = greedy.Kmax,
                   prior.family = prior.point.normal(), verbose.lvl = 1)

  fl.g$flash.fit$given.tau.dim <- 1
  fl.g$flash.fit$given.tau <- 1 / min(fl.g$residuals.sd)^2

  fl <- flashier(init = fl.g, fit = "backfit.only", tol = 10, maxiter = 500,
                 verbose.lvl = 3, output.lvl = 1)

  cat("Calculating p-values...\n")
  p.vals <- calc.p.vals(droplet, processed, fl, var.type = 1)
  cat("Done.\n")

  processed$data <- NULL
  return(c(processed, list(p.vals = p.vals, fl = fl)))
}

pseudocounts <- c(1/100, 1/16, 1/4, 1/2, 1, 2, 4, 16, 100)

sf.scran <- scran::computeSumFactors(droplet)
sf.scran <- sf.scran / median(sf.scran)

res <- list()
for (pc in pseudocounts) {
  pc.str <- as.character(pc)
  res[[pc.str]] <- do.fit(sf.scran, pc)
}

saveRDS(res, "../../output/pseudocount/pseudocount_fits.rds")
