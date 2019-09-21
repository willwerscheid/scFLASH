source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")


# 1. Fit more factors than necessary (greedily). ------------------------------

sf.scran <- scran::computeSumFactors(droplet)
sf.scran <- sf.scran / median(sf.scran)

processed <- preprocess(data = droplet,
                        min.nzcts = 10,
                        max.libsize = 30000,
                        size.factors = 0.5 * sf.scran,
                        prescale = "cells")

fl <- flashier(processed$data,
               var.type = 1,
               greedy.Kmax = 80,
               prior.family = prior.normal.scale.mix(),
               verbose.lvl = 3)


# 2. Inspect factor plots to decide on Kmax. ----------------------------------

res <- list(fl = fl,
            cell.prescaling.factors = processed$cell.prescaling.factors)
cell.types <- factor(sapply(strsplit(colnames(droplet), "_"), `[`, 3))
cell.types <- cell.types[-processed$dropped.cells]

plot.factors(res, cell.types, 1:20)
plot.factors(res, cell.types, 21:40)
plot.factors(res, cell.types, 41:60)
plot.factors(res, cell.types, 61:80) # nothing interesting here

# Factor 43 is a tuft factor I'd like to include to see whether it can
#   distinguish between tuft-1 and tuft-2 populations. So I'll go with K = 50.


# 3. Fit a backfitted flash object. -------------------------------------------

fl <- flashier(processed$data,
               var.type = 1,
               greedy.Kmax = 50,
               prior.family = prior.normal.scale.mix(),
               verbose.lvl = 1)

fl$flash.fit$given.tau.dim <- 1
fl$flash.fit$given.tau <- 1 / min(fl$residuals.sd)^2

fl <- flashier(init = fl, fit = "backfit.only", tol = 10, verbose.lvl = 3)

# The whole thing takes about 12 hours (413 backfitting iterations).
saveRDS(fl, "../../output/lowrank/full_fit50.rds")


# 4. Try using low-rank data representation (rank 100). -----------------------

set.seed(666)
lr.data <- rsvd::rsvd(processed$data, k = 100)

lr.fl <- flashier(lr.data,
                  var.type = 1,
                  greedy.Kmax = 50,
                  prior.family = prior.normal.scale.mix(),
                  verbose.lvl = 1)

lr.fl <- flashier(init = lr.fl, fit = "backfit.only", tol = 10, verbose.lvl = 3)

# Takes about 4 hours (160 backfitting iterations)
saveRDS(lr.fl, "../../output/lowrank/lowrank_fit50.rds")

