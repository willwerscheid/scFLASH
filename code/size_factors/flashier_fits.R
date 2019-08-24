source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")

sf.scran <- scran::computeSumFactors(droplet)
sf.scran <- sf.scran / median(sf.scran)

libsize <- preprocess.droplet(droplet)
noscale <- preprocess.droplet(droplet, scale.factors = rep(1, ncol(droplet)))
scran <- preprocess.droplet(droplet, scale.factors = sf.scran)

do.fit <- function(data, greedy.Kmax) {
  fl.kron1 <- flashier(data, var.type = c(1, 2), greedy.Kmax = 1)
  gene.sds <- fl.kron1$residuals.sd[[1]]
  cell.sds <- fl.kron1$residuals.sd[[2]]

  scaled.data <- data / gene.sds
  scaled.data <- t(t(scaled.data) / cell.sds)

  fl.g <- flashier(scaled.data, var.type = 1, greedy.Kmax = greedy.Kmax,
                   prior.family = prior.point.normal(), verbose.lvl = 1)

  fl.g$flash.fit$given.tau.dim <- 1
  fl.g$flash.fit$given.tau <- 1 / min(fl.g$residuals.sd)^2

  fl.b <- flashier(init = fl.g, fit = "backfit.only", tol = 10, maxiter = 500,
                   verbose.lvl = 3)

  fl.b$flash.fit <- NULL
  fl.b$sampler <- NULL

  return(list(fl = fl.b, gene.sds = gene.sds, cell.sds = cell.sds))
}

greedy.Kmax <- 20

res <- list()

res$noscale <- do.fit(noscale$data, greedy.Kmax)
res$libsize <- do.fit(libsize$data, greedy.Kmax)
res$scran <- do.fit(scran$data, greedy.Kmax)

res$noscale$sf <- noscale$scale.factors
res$libsize$sf <- libsize$scale.factors
res$scran$sf <- scran$scale.factors

res$noscale$elbo.adj <- noscale$elbo.adj
res$libsize$elbo.adj <- libsize$elbo.adj
res$scran$elbo.adj <- scran$elbo.adj

calc.p.for.scaled <- function(fl, data, sf, gene.sds, cell.sds) {
  return(calc.p.vals(data, sf,
                     mu = fitted(fl) * gene.sds * rep(cell.sds, each = nrow(data)),
                     sigma = outer(fl$residuals.sd * gene.sds, cell.sds)))
}

res$noscale$p.vals <- calc.p.for.scaled(res$noscale$fl,
                                        noscale$data,
                                        res$noscale$sf,
                                        res$noscale$gene.sds,
                                        res$noscale$cell.sds)
res$libsize$p.vals <- calc.p.for.scaled(res$libsize$fl,
                                        libsize$data,
                                        res$libsize$sf,
                                        res$libsize$gene.sds,
                                        res$libsize$cell.sds)
res$scran$p.vals <- calc.p.for.scaled(res$scran$fl,
                                      scran$data,
                                      res$scran$sf,
                                      res$scran$gene.sds,
                                      res$scran$cell.sds)

saveRDS(res, "../../output/size_factors/flashier_fits.rds")

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
