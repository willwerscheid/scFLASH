source("../utils.R")

droplet <- readRDS("../../data/droplet.rds")
processed <- preprocess(droplet,
                        min.nzcts = 10,
                        max.libsize = 30000,
                        prescale = "both")

greedy.Kmax <- 30

do.fit <- function(prior.family) {
  sink("tmp_g.txt")
  t.g <- system.time({
    fl.g <- flashier(processed$data, var.type = 1, greedy.Kmax = greedy.Kmax,
                     prior.family = prior.family, verbose.lvl = -1)
  })
  sink()

  fl.g$flash.fit$given.tau.dim <- 1
  fl.g$flash.fit$given.tau <- 1 / min(fl.g$residuals.sd)^2

  sink("tmp_b.txt")
  t.b <- system.time({
    fl.b <- flashier(init = fl.g, fit = "backfit.only", tol = 10, maxiter = 500,
                     verbose.lvl = -1)
  })
  sink()

  output.g <- data.table::fread("tmp_g.txt")
  output.b <- data.table::fread("tmp_b.txt")

  file.remove("tmp_g.txt")
  file.remove("tmp_b.txt")

  # Only keep the final iteration for each greedy factor.
  output.g <- subset(output.g, c((output.g$Factor[1:(nrow(output.g) - 1)]
                                  != output.g$Factor[2:nrow(output.g)]), TRUE))

  fl.b$flash.fit <- NULL
  fl.b$sampler <- NULL

  return(list(t = list(greedy = t.g[3], backfit = t.b[3]),
              fl = fl.b,
              output = rbind(output.g, output.b)))
}

res <- list()

cat("Fitting point-normal priors...\n")
res$pn <- do.fit(prior.point.normal())
cat("Fitting scale-mixture-of-normal priors...\n")
res$ash <- do.fit(prior.normal.scale.mix())
cat("Fitting nonnegative priors on cells...\n")
res$snn.cell <- do.fit(c(prior.point.normal(), prior.nonnegative()))
cat("Fitting nonnegative priors on genes...\n")
res$snn.gene <- do.fit(c(prior.nonnegative(), prior.point.normal()))

saveRDS(res, "../../output/prior_type/priortype_fits.rds")
