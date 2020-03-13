library(tidyverse)
library(flashier)

source("../utils.R")

pbmc <- readRDS("../../data/10x/pbmc.rds")
processed <- preprocess(pbmc,
                        min.nzcts = 10,
                        max.libsize = Inf,
                        prescale = "none")

greedy.Kmax <- 20

var.lower.bd <- pbmc[-processed$dropped.genes, ] / (exp(processed$data))^2
var.lower.bd <- t(t(var.lower.bd) / processed$size.factors^2)
var.lower.bd <- apply(var.lower.bd, 1, mean)
var.lower.bd <- min(var.lower.bd)

t.unreg <- system.time({
  fl.unreg <- flash.init(processed$data, S = sqrt(var.lower.bd), var.type = 1,
                         regularize.var = FALSE) %>%
    flash.add.greedy(Kmax = greedy.Kmax) %>%
    flash.backfit(maxiter = 100, verbose.lvl = 3)
})

cat("Calculating p-values...\n")
p.vals.unreg <- calc.p.vals(pbmc, processed, fl.unreg, var.type = 1)
cat("Done.\n")

fl.unreg$flash.fit$Y <- NULL
fl.unreg$sampler <- NULL

t.reg <- system.time({
  fl.reg <- flash.init(processed$data, S = sqrt(var.lower.bd), var.type = 1,
                       regularize.var = TRUE) %>%
    flash.add.greedy(Kmax = greedy.Kmax) %>%
    flash.backfit(maxiter = 100, verbose.lvl = 3)
})

cat("Calculating p-values...\n")
p.vals.reg <- calc.p.vals(pbmc, processed, fl.reg, var.type = 1)
cat("Done.\n")

fl.reg$flash.fit$Y <- NULL
fl.reg$sampler <- NULL

res <- list(unreg = list(elapsed.time = t.unreg[3], p.vals = p.vals.unreg, fl = fl.unreg),
            reg = list(elapsed.time = t.reg[3], p.vals = p.vals.reg, fl = fl.reg))

saveRDS(res, "../../output/var_reg/varreg_fits2.rds")
