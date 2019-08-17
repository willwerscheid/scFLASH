source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")
processed <- preprocess.droplet(droplet)

greedy.Kmax <- 20

fl.kron1 <- flashier(processed$data, var.type = c(1, 2), greedy.Kmax = 1)
gene.sds <- fl.kron1$residuals.sd[[1]]
cell.sds <- fl.kron1$residuals.sd[[2]]

scaled.data <- processed$data / gene.sds
scaled.data <- t(t(scaled.data) / cell.sds)

sink("../../output/backfit/greedy_output.txt")
t.g <- system.time({
  fl.g <- flashier(scaled.data, var.type = 1, greedy.Kmax = greedy.Kmax,
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

res <- list(fl = list(greedy = fl.g, backfit = fl.b),
            t = list(greedy = t.g[3], backfit = t.b[3]))
saveRDS(res, "../../output/backfit/flashier_fits.rds")
