source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")

do.fit <- function(var.type, prescale = "none", greedy.Kmax = 20) {
  cat("FITTING: var.type =", var.type, "; prescale =", prescale, "\n")

  t <- system.time({
    processed <- do.call(preprocess,
                         c(list(data = droplet,
                                min.nzcts = 10,
                                max.libsize = 30000,
                                prescale = prescale)))

    fl <- do.call(flashier,
                  c(list(data = processed$data,
                         var.type = var.type,
                         greedy.Kmax = greedy.Kmax,
                         output.lvl = 1)))
  })

  cat("Calculating p-values...\n")
  p.vals <- calc.p.vals(droplet, processed, fl, var.type)
  cat("Done.\n")

  processed$data <- NULL
  return(c(processed, list(elapsed.time = t[3], p.vals = p.vals, fl = fl)))
}

res <- list()
res$const <- do.fit(var.type = 0)
res$bygene <- do.fit(var.type = 1)
res$bycell <- do.fit(var.type = 2)
res$kron <- do.fit(var.type = c(1, 2))
res$scale.cells <- do.fit(var.type = 1, prescale = "cells")
res$scale.genes <- do.fit(var.type = 2, prescale = "genes")
res$scale.both <- do.fit(var.type = 1, prescale = "both")

saveRDS(res, "../../output/var_type/vartype_fits.rds")
