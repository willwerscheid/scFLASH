source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")
processed <- preprocess.droplet(droplet)

greedy.Kmax <- 20

cat("Fitting constant variance structure...\n")
t.const <- system.time({
  fl.const <- flashier(processed$data, var.type = 0, greedy.Kmax = greedy.Kmax,
                       output.lvl = 1)
})

cat("Fitting genewise variance structure...\n")
t.bygene <- system.time({
  fl.bygene <- flashier(processed$data, var.type = 1, greedy.Kmax = greedy.Kmax,
                        output.lvl = 1)
})

cat("Fitting cellwise variance structure...\n")
t.bycell <- system.time({
  fl.bycell <- flashier(processed$data, var.type = 2, greedy.Kmax = greedy.Kmax,
                        output.lvl = 1)
})

cat("Fitting Kronecker variance structure...\n")
t.kron <- system.time({
  fl.kron <- flashier(processed$data, var.type = c(1, 2), greedy.Kmax = greedy.Kmax,
                      output.lvl = 1)
})

cat("Fitting approximate Kroncker variance structure (scale cells)...\n")
t.scale.cells <- system.time({
  fl.kron1 <- flashier(processed$data, var.type = c(1, 2), greedy.Kmax = 1)
  cell.sds <- fl.kron1$residuals.sd[[2]]
  scaled.data <- t(t(processed$data) / cell.sds)
  fl.scale.cells <- flashier(scaled.data, var.type = 1, greedy.Kmax = greedy.Kmax,
                             output.lvl = 1)
})

cat("Fitting approximate Kronecker variance structure (scale genes)...\n")
t.scale.genes <- system.time({
  fl.kron1 <- flashier(processed$data, var.type = c(1, 2), greedy.Kmax = 1)
  gene.sds <- fl.kron1$residuals.sd[[1]]
  scaled.data <- processed$data / gene.sds
  fl.scale.genes <- flashier(scaled.data, var.type = 2, greedy.Kmax = greedy.Kmax,
                             output.lvl = 1)
})

cat("Fitting approximate Kronecker variance structure (scale both)...\n")
t.scale.both <- system.time({
  fl.kron1 <- flashier(processed$data, var.type = c(1, 2), greedy.Kmax = 1)
  gene.sds <- fl.kron1$residuals.sd[[1]]
  cell.sds <- fl.kron1$residuals.sd[[2]]
  scaled.data <- processed$data / gene.sds
  scaled.data <- t(t(scaled.data) / cell.sds)
  fl.scale.both <- flashier(scaled.data, var.type = 1, greedy.Kmax = greedy.Kmax,
                             output.lvl = 1)
})

fl.const$flash.fit <- NULL
fl.bygene$flash.fit <- NULL
fl.bycell$flash.fit <- NULL
fl.kron$flash.fit <- NULL
fl.scale.cells$flash.fit <- NULL
fl.scale.genes$flash.fit <- NULL
fl.scale.both$flash.fit <- NULL

res <- list(fl = list(const = fl.const,
                      bygene = fl.bygene,
                      bycell = fl.bycell,
                      kron = fl.kron,
                      scale.cells = fl.scale.cells,
                      scale.genes = fl.scale.genes,
                      scale.both = fl.scale.both),
            t = list(const = t.const[3],
                     bygene = t.bygene[3],
                     bycell = t.bycell[3],
                     kron = t.kron[3],
                     scale.cells = t.scale.cells[3],
                     scale.genes = t.scale.genes[3],
                     scale.both = t.scale.both[3]),
            kron1 = list(gene.sds = gene.sds,
                         cell.sds = cell.sds))

saveRDS(res, "../../output/var_type/flashier_fits.rds")

cat("Calculating p-values (constant variance)...\n")
const.p.vals <- calc.p.vals(processed$data,
                            processed$scale.factors,
                            mu = fitted(res$fl$const),
                            sigma = res$fl$const$residuals.sd)

cat("Calculating p-values (gene-wise variance)...\n")
bygene.p.vals <- calc.p.vals(processed$data,
                             processed$scale.factors,
                             mu = fitted(res$fl$bygene),
                             sigma = res$fl$bygene$residuals.sd)

cat("Calculating p-values (cell-wise variance)...\n")
bycell.p.vals <- calc.p.vals(processed$data,
                             processed$scale.factors,
                             mu = fitted(res$fl$bycell),
                             sigma = rep(res$fl$bycell$residuals.sd,
                                         each = nrow(processed$data)))

cat("Calculating p-values (Kronecker variance)...\n")
kron.p.vals <- calc.p.vals(processed$data,
                           processed$scale.factors,
                           mu = fitted(res$fl$kron),
                           sigma = outer(res$fl$kron$residuals.sd[[1]],
                                         res$fl$kron$residuals.sd[[2]]))

cat("Calculating p-values (pre-scaling cells)...\n")
scale.cells.p.vals <- calc.p.vals(processed$data,
                                  processed$scale.factors,
                                  mu = (fitted(res$fl$scale.cells)
                                        * rep(res$kron1$cell.sds,
                                              each = nrow(processed$data))),
                                  sigma = outer(res$fl$scale.cells$residuals.sd,
                                                res$kron1$cell.sds))

cat("Calculating p-values (pre-scaling genes)...\n")
scale.genes.p.vals <- calc.p.vals(processed$data,
                                  processed$scale.factors,
                                  mu = fitted(res$fl$scale.genes) * res$kron1$gene.sds,
                                  sigma = outer(res$kron1$gene.sds,
                                                res$fl$scale.genes$residuals.sd))

cat("Calculating p-values (pre-scaling both)...\n")
scale.both.p.vals <- calc.p.vals(processed$data,
                                 processed$scale.factors,
                                 mu = (fitted(res$fl$scale.both) * res$kron1$gene.sds
                                       * rep(res$kron1$cell.sds,
                                             each = nrow(processed$data))),
                                 sigma = outer((res$fl$scale.both$residuals.sd
                                                * res$kron1$gene.sds),
                                               res$kron1$cell.sds))

res$pvals <- list(const = const.p.vals,
                  bygene = bygene.p.vals,
                  bycell = bycell.p.vals,
                  kron = kron.p.vals,
                  scale.cells = scale.cells.p.vals,
                  scale.genes = scale.genes.p.vals,
                  scale.both = scale.both.p.vals)

saveRDS(res, "../../output/var_type/flashier_fits.rds")
