source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")

do.fit <- function(size.factors, title, greedy.Kmax = 20) {
  cat("FITTING:", title, "\n")

  processed <- do.call(preprocess,
                       c(list(data = droplet,
                              min.nzcts = 10,
                              max.libsize = 30000,
                              size.factors = size.factors,
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

sf.scran <- scran::computeSumFactors(droplet)
sf.scran <- sf.scran / median(sf.scran)

res <- list()
res$noscale <- do.fit(rep(1, ncol(droplet)), "No size factors")
res$libsize <- do.fit(NULL, "Library-size normalization")
res$scran <- do.fit(sf.scran, "scran size factors")

saveRDS(res, "../../output/size_factors/sizefactor_fits.rds")


# SCnorm decides that all scale factors are 1...?!

# t.scnorm <- system.time({
#   scnorm.sf <- SCnorm::SCnorm(droplet,
#                               Conditions = rep(1, ncol(droplet)),
#                               FilterCellNum = 10,
#                               NCores = 4,
#                               PrintProgressPlots = TRUE,
#                               reportSF = TRUE)
# })
# scnorm.sfsf <- SCnorm::results(scnorm.sf, type = "ScaleFactors")
