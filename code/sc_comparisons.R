# library(ashr) -- uncomment after negative mixture proportions bug is fixed
devtools::load_all("~/Github/ashr")
# library(flashier) -- uncomment after pushing 0.1.1 updates to master
devtools::load_all("~/Github/flashier")
library(Matrix)
library(ggplot2)


# Load data ---------------------------------------------------------------

# The drop-seq dataset in Montoro et al. can be downloaded here:
#     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354
trachea <- read.table("./data/GSE103354_Trachea_droplet_UMIcounts.txt")
trachea <- as.matrix(trachea)
trachea <- Matrix(trachea)


# Set test parameters -----------------------------------------------------

# Number of times to subsample data and run fits.
ntrials <- 50

# Size of subsampled data matrices.
ncells <- 200
ngenes <- 500

# Minimum number of nonzero counts needed to include a gene.
min.cts <- 3

# Proportion of entries to impute.
prop.missing <- 0.01
nmissing <- ceiling(prop.missing * ncells * ngenes)

# Number of factors to add to each fit.
K <- 5

# Verbose flag for flashier output.
verbose <- FALSE


# Set up data frames ------------------------------------------------------

create.df <- function(names) {
  names <- paste0(names, ".")
  df.names <- outer(names, c("mse", "elbo"), FUN = paste0)
  df <- data.frame(rep(list(numeric(0)), length(df.names)))
  names(df) <- t(df.names)
  return(df)
}

# Test 1: variance structure.
var.df <- create.df(c("constant", "genewise", "fixed", "noisyA", "noisyB"))

# Test 2: data transformation.
trans.df <- create.df(c("log1p", "anscombe", "arcsin", "raw", "pearson"))

# Test 3: normalization method.
norm.df <- create.df(c("none", "fitmean", "scale"))

# Test 4: prior type.
prior.df <- create.df(c("normal", "nngeneA", "nngeneB", "nncellA", "nncellB"))

# Test 5: number of factors and backfit.
nfactors.df <- create.df(c("g1", "bf1", "g2", "bf2", "g3", "bf3"))


# Functions for populating data frames ------------------------------------

get.fitted <- function(fl) {
  # Get fitted values on the log1p scale.
  fitted.vals <- flashier:::lowrank.expand(fl$fit$EF)
  # Convert to raw counts.
  return(exp(fitted.vals) - 1)
}

log1p.mse <- function(true, preds) {
  # Threshold out negative counts.
  preds <- pmax(preds, 0)
  return(mean((log1p(true) - log1p(preds))^2))
}

get.elbo <- function(fl, adjustment) {
  return(fl$objective + adjustment)
}

all.metrics <- function(data, missing.vals,
                        fl, fitted.vals, elbo.adjustment) {
  missing.idx <- which(is.na(data))
  if (is.null(fitted.vals)) {
    fitted.vals <- get.fitted(fl)
  }
  return(list(mse = log1p.mse(missing.vals, fitted.vals[missing.idx]),
              elbo = get.elbo(fl, elbo.adjustment)))
}

get.df.row <- function(data, missing.vals,
                       fl.list,
                       mean.factors = NULL,
                       fitted.vals = NULL,
                       elbo.adjustment = NULL) {
  if (is.null(fitted.vals)) {
    fitted.vals = rep(list(NULL), length(fl.list))
  }
  if (is.null(elbo.adjustment)) {
    elbo.adjustment = rep(list(0), length(fl.list))
  }

  return(unlist(mapply(all.metrics,
                       fl = fl.list,
                       fitted.vals = fitted.vals,
                       elbo.adjustment = elbo.adjustment,
                       MoreArgs = list(data = data,
                                       missing.vals = missing.vals))))
}


# Run tests ---------------------------------------------------------------

for (i in 1:ntrials) {
  cat("TRIAL", i, "\n")

  set.seed(i)
  rand.cells <- sample(1:ncol(trachea), ncells)
  samp <- trachea[, rand.cells]

  rand.genes <- sample(which(rowSums(samp > 0) >= min.cts), ngenes)
  samp <- samp[rand.genes, ]

  missing.idx <- sample(1:length(samp), nmissing)
  missing.vals <- samp[missing.idx]
  samp[missing.idx] <- NA

  # Test 1: variance structure.
  cat("  Running variance structure tests.\n")

  fl.var0 <- flashier(log1p(samp), var.type = 0,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 1L * verbose)

  fl.var1 <- flashier(log1p(samp), var.type = 1,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 1L * verbose)

  S <- sqrt(samp) / (samp + 1)
  nz.prop <- sum(samp > 0, na.rm = TRUE) / (length(samp) - length(missing.idx))
  S[S == 0] <- sqrt(nz.prop) / (nz.prop + 1)

  fl.fixS <- flashier(log1p(samp), S = S,
                      var.type = NULL,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 1L * verbose)

  suppressMessages({
    fl.noisyA <- flashier(log1p(samp), S = S,
                          var.type = 0,
                          prior.type = "normal.mix",
                          greedy.Kmax = K + 1, backfit = "none",
                          verbose.lvl = 1L * verbose)
  })

  suppressMessages({
    fl.noisyB <- flashier(log1p(samp), S = sqrt(samp) / (samp + 1),
                          var.type = 0,
                          prior.type = "normal.mix",
                          greedy.Kmax = K + 1, backfit = "none",
                          verbose.lvl = 1L * verbose)
  })

  var.df[i, ] <- get.df.row(samp, missing.vals,
                            fl.list = list(fl.var0, fl.var1, fl.fixS,
                                           fl.noisyA, fl.noisyB))

  # Test 2: data transformation.
  cat("  Running data transformation tests.\n")

  fl.log1p <- fl.var0
  log1p.adj <- -sum(log1p(samp), na.rm = TRUE)

  fl.ans <- flashier(sqrt(samp + 0.375), var.type = 0,
                     prior.type = "normal.mix",
                     greedy.Kmax = K + 1, backfit = "none",
                     verbose.lvl = 1L * verbose)
  ans.fitted <- flashier:::lowrank.expand(fl.ans$fit$EF)^2 - 0.375
  ans.adj <- -sum(log(2) + 0.5 * log(samp + 0.375), na.rm = TRUE)

  cell.sums <- colSums(samp, na.rm = TRUE)
  props <- samp / rep(cell.sums, each = nrow(samp))
  fl.arcsin <- flashier(asin(sqrt(props)), var.type = 0,
                        prior.type = "normal.mix",
                        greedy.Kmax = K + 1, backfit = "none",
                        verbose.lvl = 1L * verbose)
  arcsin.fitted <- (sin(flashier:::lowrank.expand(fl.arcsin$fit$EF)^2)
                    * rep(cell.sums, each = nrow(samp)))
  arcsin.adj <- NA

  fl.raw <- flashier(props, var.type = 0,
                     prior.type = "normal.mix",
                     greedy.Kmax = K + 1, backfit = "none",
                     verbose.lvl = 1L * verbose)
  raw.fitted <- (flashier:::lowrank.expand(fl.raw$fit$EF)
                 * rep(cell.sums, each = nrow(samp)))
  raw.adj <- -sum(rep(log(cell.sums), each = nrow(samp))[-missing.idx])

  gene.props <- rowSums(samp, na.rm = TRUE) / sum(samp, na.rm = TRUE)
  mu <- outer(gene.props, cell.sums)
  sd.mat <- sqrt(mu - mu^2 / rep(cell.sums, each = nrow(samp)))
  resid <- (samp - mu) / sd.mat

  fl.pearson <- flashier(resid, var.type = 0,
                         prior.type = "normal.mix",
                         greedy.Kmax = K + 1, backfit = "none",
                         verbose.lvl = 1L * verbose)
  pearson.fitted <- mu + sd.mat * flashier:::lowrank.expand(fl.pearson$fit$EF)
  pearson.adj <- -sum(log(sd.mat)[-missing.idx])

  trans.df[i, ] <- get.df.row(samp, missing.vals,
                              fl.list = list(fl.log1p, fl.ans, fl.arcsin,
                                             fl.raw, fl.pearson),
                              fitted.vals = list(get.fitted(fl.log1p),
                                                 ans.fitted, arcsin.fitted,
                                                 raw.fitted, pearson.fitted),
                              elbo.adjustment = list(log1p.adj,
                                                     ans.adj, arcsin.adj,
                                                     raw.adj, pearson.adj))

  # Test 3: scaling method.
  cat("  Running scaling tests.\n")

  fl.none <- fl.var0
  none.adj <- log1p.adj

  fl.ones <- flashier(log1p(samp), var.type = 0,
                      prior.type = "normal.mix",
                      fixed.factors = c(ones.factor(2), ones.factor(1)),
                      greedy.Kmax = K,
                      backfit.after = 2, final.backfit = FALSE,
                      verbose.lvl = 1L * verbose)
  ones.adj <- log1p.adj

  scaled.samp <- samp * median(cell.sums) / rep(cell.sums, each = nrow(samp))
  fl.scale <- flashier(log1p(scaled.samp), var.type = 0,
                       prior.type = "normal.mix",
                       fixed.factors = ones.factor(2),
                       greedy.Kmax = K, backfit = "none",
                       verbose.lvl = 1L * verbose)
  scale.fitted <- (get.fitted(fl.scale)
                   * rep(cell.sums, each = nrow(samp)) / median(cell.sums))
  scale.adj <- sum(log(median(cell.sums)) - rep(log(cell.sums), each = nrow(samp))
                   - log1p(scaled.samp), na.rm = TRUE)

  norm.df[i, ] <- get.df.row(samp, missing.vals,
                             fl.list = list(fl.none, fl.ones, fl.scale),
                             fitted.vals = list(get.fitted(fl.none),
                                                get.fitted(fl.ones),
                                                scale.fitted),
                             elbo.adjustment = list(none.adj, ones.adj,
                                                    scale.adj))

  # Test 4: prior type.
  cat("  Running prior type tests.\n")

  fl.normalmix <- fl.var0

  fl.nngenes <- flashier(log1p(samp), var.type = 0,
                         prior.type = c("nonnegative", "normal.mix"),
                         greedy.Kmax = K + 1, backfit = "none",
                         verbose.lvl = 1L * verbose)

  fl.nngenes.pm <- flashier(log1p(samp), var.type = 0,
                            prior.type = c("nonnegative", "normal.mix"),
                            ash.param = list(method = "fdr"),
                            greedy.Kmax = K + 1, backfit = "none",
                            verbose.lvl = 1L * verbose)

  fl.nncells <- flashier(log1p(samp), var.type = 0,
                         prior.type = c("normal.mix", "nonnegative"),
                         greedy.Kmax = K + 1, backfit = "none",
                         verbose.lvl = 1L * verbose)

  fl.nncells.pm <- flashier(log1p(samp), var.type = 0,
                            prior.type = c("normal.mix", "nonnegative"),
                            ash.param = list(method = "fdr"),
                            greedy.Kmax = K + 1, backfit = "none",
                            verbose.lvl = 1L * verbose)

  prior.df[i, ] <- get.df.row(samp, missing.vals,
                              fl.list = list(fl.normalmix,
                                             fl.nngenes, fl.nngenes.pm,
                                             fl.nncells, fl.nncells.pm))

  # Test 5: number of factors and backfit.
  cat("  Running backfitting tests.\n")

  fl.g <- fl.var0

  fl.b <- flashier(flash.init = fl.g, backfit = "only",
                   backfit.reltol = 10,
                   verbose.lvl = 3L * verbose)

  # Add K more factors, then K more.
  fl.g2 <- flashier(flash.init = fl.g,
                    greedy.Kmax = K, backfit = "none",
                    verbose.lvl = 1L * verbose)
  fl.b2 <- flashier(flash.init = fl.b,
                    greedy.Kmax = K, backfit = "final",
                    backfit.reltol = 10,
                    verbose.lvl = 3L * verbose)
  fl.g3 <- flashier(flash.init = fl.g2,
                    greedy.Kmax = K, backfit = "none",
                    verbose.lvl = 1L * verbose)
  fl.b3 <- flashier(flash.init = fl.b2,
                    greedy.Kmax = K, backfit = "final",
                    backfit.reltol = 10,
                    verbose.lvl = 3L * verbose)

  nfactors.df[i, ] <- get.df.row(samp, missing.vals,
                                 fl.list = list(fl.g, fl.b,
                                                fl.g2, fl.b2,
                                                fl.g3, fl.b3))
}

# Test 6: incremental addition of factors

set.seed(666)
data <- log1p(trachea[rowSums(trachea > 0) >= min.cts, ])
missing.idx <- sample(1:length(data), ceiling(prop.missing * length(data)))
true.vals <- data[missing.idx]
data[missing.idx] <- NA

Kmax <- 60

fl <- flashier(data, var.type = 0,
               prior.type = "normal.mix",
               fixed.factors = c(ones.factor(1), ones.factor(2)),
               greedy.Kmax = 1, backfit = "none",
               final.nullchk = FALSE,
               verbose.lvl = 3)
mse.vec <- mse(true.vals, preds(fl, missing.idx), 1:length(true.vals))

for (k in 2:Kmax) {
  fl <- flashier(flash.init = fl,
                 greedy.Kmax = 1, backfit = "none",
                 final.nullchk = FALSE,
                 verbose.lvl = 3)
  mse.vec <- c(mse.vec,
               mse(true.vals, preds(fl, missing.idx), 1:length(true.vals)))
}

mse.diff <- c(NA, mse.vec[2:length(mse.vec)] - mse.vec[1:(length(mse.vec) - 1)])
mse.df <- data.frame(k = 1:Kmax, mse = mse.vec, mse.diff = mse.diff)


# Save results ------------------------------------------------------------

all.res <- list(var.df = var.df,
                trans.df = trans.df,
                norm.df = norm.df,
                prior.df = prior.df,
                nfactors.df = nfactors.df,
                mse.df = mse.df)
saveRDS(all.res, "./output/sc_comparisons/allres.rds")

