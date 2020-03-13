library(tidyverse)
library(flashier)

source("../utils.R")

pbmc <- readRDS("../../data/10x/pbmc.rds")
processed <- preprocess(pbmc,
                        min.nzcts = 10,
                        max.libsize = Inf,
                        prescale = "cells")

greedy.Kmax <- 20

do.fit <- function(prior.scale, greedy.Kmax = 20) {
  cat("FITTING: prior.scale =", prior.scale, "\n")

  t <- system.time({
    fl <- flash(processed$data, var.type = 1, greedy.Kmax = greedy.Kmax,
                regularize.var = prior.scale, verbose.lvl = 1) %>%
      flash.backfit(maxiter = 100, verbose.lvl = 3)
  })

  cat("Calculating p-values...\n")
  p.vals <- calc.p.vals(pbmc, processed, fl, var.type = 1)
  cat("Done.\n")

  processed$data <- NULL
  return(c(processed, list(elapsed.time = t[3], p.vals = p.vals, fl = fl)))
}

res <- list()
res$scale1 <- do.fit(prior.scale = 1)
res$scale10 <- do.fit(prior.scale = 10)
res$scale100 <- do.fit(prior.scale = 100)
res$scale1000 <- do.fit(prior.scale = 1000)
res$scale10000 <- do.fit(prior.scale = 10000)
res$scale100000 <- do.fit(prior.scale = 100000)

saveRDS(res, "../../output/var_reg/varreg_fits.rds")
