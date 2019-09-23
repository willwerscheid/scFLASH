source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")

sf.scran <- scran::computeSumFactors(droplet)
sf.scran <- sf.scran / median(sf.scran)

processed <- preprocess(data = droplet,
                        min.nzcts = 10,
                        max.libsize = 30000,
                        size.factors = 0.5 * sf.scran,
                        prescale = "cells")

# 1. Fit more factors than necessary (greedily). ------------------------------

# fl <- flashier(processed$data,
#                var.type = 1,
#                greedy.Kmax = 80,
#                prior.family = c(prior.normal.scale.mix(), prior.nonnegative()),
#                verbose.lvl = 3)

# 2. Inspect factor plots to decide on Kmax. ----------------------------------

# res <- list(fl = fl,
#             cell.prescaling.factors = processed$cell.prescaling.factors)
# cell.types <- factor(sapply(strsplit(colnames(droplet), "_"), `[`, 3))
# cell.types <- cell.types[-processed$dropped.cells]
#
# plot.factors(res, cell.types, 1:20)
# plot.factors(res, cell.types, 21:40)
# plot.factors(res, cell.types, 41:60)
# plot.factors(res, cell.types, 61:80) # nothing interesting here

# 3. Fit a backfitted flash object. -------------------------------------------

sink("tmp_g.txt")
t.g <- system.time({
  fl <- flashier(processed$data,
                 var.type = 1,
                 greedy.Kmax = 30,
                 prior.family = c(prior.normal.scale.mix(), prior.nonnegative()),
                 verbose.lvl = -1)
})
sink()

fl$flash.fit$given.tau.dim <- 1
fl$flash.fit$given.tau <- 1 / min(fl$residuals.sd)^2

sink("tmp_b.txt")
t.b <- system.time({
  fl <- flashier(init = fl, fit = "backfit.only", tol = 10, verbose.lvl = -1)
})
sink()

output.g <- data.table::fread("tmp_g.txt")
output.b <- data.table::fread("tmp_b.txt")

file.remove("tmp_g.txt")
file.remove("tmp_b.txt")

# Only keep the final iteration for each greedy factor.
output.g <- subset(output.g, c((output.g$Factor[1:(nrow(output.g) - 1)]
                                != output.g$Factor[2:nrow(output.g)]), TRUE))

fl$flash.fit <- NULL
fl$sampler <- NULL

res <- list(t = list(greedy = t.g[3], backfit = t.b[3]),
            cell.prescaling.factors = processed$cell.prescaling.factors,
            fl = fl,
            output = rbind(output.g, output.b))

saveRDS(res, "../../output/final_montoro/final_fit.rds")
