source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")
processed <- preprocess(droplet,
                        min.nzcts = 10,
                        max.libsize = 30000,
                        prescale = "cells")

greedy.Kmax <- 20

sink("../../output/backfit/greedy_output.txt")
t.g <- system.time({
  fl.g <- flashier(processed$data, var.type = 1, greedy.Kmax = greedy.Kmax,
                   verbose.lvl = -1)
})
sink()

# A hack to prevent tau from going to infinity. Sets a lower bound on the
#   gene-wise residual sds.
fl.g$flash.fit$given.tau.dim <- 1
fl.g$flash.fit$given.tau <- 1 / min(fl.g$residuals.sd)^2

sink("../../output/backfit/backfit_output.txt")
t.b <- system.time({
  fl.b <- flashier(init = fl.g, fit = "backfit.only", maxiter = 500,
                   verbose.lvl = -1)
})
sink()

# Remove large objects from results before saving to file.
fl.g$flash.fit <- NULL
fl.b$flash.fit <- NULL
fl.g$sampler <- NULL
fl.b$sampler <- NULL

res <- list(greedy = list(elapsed.time = t.g[3], fl = fl.g),
            backfit = list(elapsed.time = t.b[3], fl = fl.b))

saveRDS(res, "../../output/backfit/backfit_fits.rds")
