# library(ashr) -- uncomment after negative mixture proportions bug is fixed
devtools::load_all("~/Github/ashr")
# library(flashier) -- uncomment after pushing 0.1.1 updates to master
devtools::load_all("~/Github/flashier")
library(Matrix)
library(ggplot2)


# Load data ---------------------------------------------------------------

# Only keep basal cells from the drop-seq dataset in Montoro et al. The data
#   can be downloaded here:
#     https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354
trachea <- read.table("./data/GSE103354_Trachea_droplet_UMIcounts.txt")
trachea <- as.matrix(trachea)
trachea <- Matrix(trachea)

cell.names <- colnames(trachea)
cell.types <- as.factor(sapply(strsplit(cell.names, "_"), `[`, 3))
basal <- trachea[, cell.types == 'Basal']


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
  df.names <- outer(names, c("overall", "zero", "nonzero"), FUN = paste0)
  df <- data.frame(rep(list(numeric(0)), length(df.names)))
  names(df) <- t(df.names)
  return(df)
}

# Test 1: variance structure.
var.df <- create.df(c("constant", "genewise", "fixed", "noisyA", "noisyB"))

# Test 2: data transformation.
trans.df <- create.df(c("log1p", "anscombe", "arcsin", "raw"))

# Test 3: normalization method.
norm.df <- create.df(c("none", "fitmean", "scale"))

# Test 4: prior type.
prior.df <- create.df(c("normal", "nngeneA", "nngeneB", "nncellA", "nncellB"))

# Test 5: number of factors and backfit.
nfactors.df <- create.df(c("g1", "bf1", "g2", "bf2", "g3", "bf3"))


# Create functions for populating data frames -----------------------------

preds <- function(fl, idx) {
  return(flashier:::lowrank.expand(fl$fit$EF)[idx])
}

mse <- function(true, preds, idx) {
  return(mean((true[idx] - preds[idx])^2))
}

mse.by.lvl <- function(true, preds) {
  zero.idx <- (true.vals == 0)
  return(list(overall.mse = mse(true, preds, rep(TRUE, length(preds))),
              zero.mse = mse(true, preds, zero.idx),
              nz.mse = mse(true, preds, !zero.idx)))
}

all.mse <- function(true, preds.list) {
  return(unlist(lapply(preds.list, function(preds) mse.by.lvl(true, preds))))
}


# Run tests ---------------------------------------------------------------

for (i in 1:ntrials) {
  cat("TRIAL", i, "\n")

  set.seed(i)
  rand.cells <- sample(1:ncol(basal), ncells)
  samp <- basal[, rand.cells]

  rand.genes <- sample(which(rowSums(samp > 0) >= min.cts), ngenes)
  samp <- samp[rand.genes, ]

  missing.idx <- sample(length(samp), nmissing)
  true.vals <- log1p(samp[missing.idx])
  samp[missing.idx] <- NA

  # Test 1: variance structure.
  cat("  Running variance structure tests.\n")

  fl.var0 <- flashier(log1p(samp), var.type = 0,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 2L * verbose)

  fl.var1 <- flashier(log1p(samp), var.type = 1,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 2L * verbose)

  S <- sqrt(samp) / (samp + 1)
  nz.prop <- sum(samp > 0, na.rm = TRUE) / (length(samp) - length(missing.idx))
  S[S == 0] <- sqrt(nz.prop) / (nz.prop + 1)

  fl.fixS <- flashier(log1p(samp), S = S,
                      var.type = NULL,
                      prior.type = "normal.mix",
                      greedy.Kmax = K + 1, backfit = "none",
                      verbose.lvl = 2L * verbose)

  suppressMessages({
    fl.noisyA <- flashier(log1p(samp), S = S,
                          var.type = 0,
                          prior.type = "normal.mix",
                          greedy.Kmax = K + 1, backfit = "none",
                          verbose.lvl = 2L * verbose)
  })

  suppressMessages({
    fl.noisyB <- flashier(log1p(samp), S = sqrt(samp) / (samp + 1),
                          var.type = 0,
                          prior.type = "normal.mix",
                          greedy.Kmax = K + 1, backfit = "none",
                          verbose.lvl = 2L * verbose)
  })

  var.df[i, ] <- all.mse(true.vals,
                         list(preds(fl.var0, missing.idx),
                              preds(fl.var1, missing.idx),
                              preds(fl.fixS, missing.idx),
                              preds(fl.noisyA, missing.idx),
                              preds(fl.noisyB, missing.idx)))

  # Test 2: data transformation.
  cat("  Running data transformation tests.\n")

  fl.log1p <- fl.var0

  fl.ans <- flashier(sqrt(samp + 0.375), var.type = 0,
                     prior.type = "normal.mix",
                     greedy.Kmax = K + 1, backfit = "none",
                     verbose.lvl = 2L * verbose)
  ans.preds <- log1p(preds(fl.ans, missing.idx)^2 - 0.375)

  cell.sums <- colSums(samp, na.rm = TRUE)
  props <- samp / rep(cell.sums, each = nrow(samp))
  fl.arcsin <- flashier(asin(sqrt(props)), var.type = 0,
                        prior.type = "normal.mix",
                        greedy.Kmax = K, backfit = "none",
                        verbose.lvl = 2L * verbose)
  missing.cols <- col(samp)[missing.idx]
  arcsin.preds <- log1p(sin(preds(fl.arcsin, missing.idx))^2
                        * cell.sums[missing.cols])

  fl.raw <- flashier(props, var.type = 0,
                     prior.type = "normal.mix",
                     greedy.Kmax = K,
                     backfit.after = 2, final.backfit = FALSE,
                     verbose.lvl = 2L * verbose)
  raw.preds <- log1p(preds(fl.raw, missing.idx) * cell.sums[missing.cols])
  raw.preds[is.nan(raw.preds)] <- 0

  trans.df[i, ] <- all.mse(true.vals,
                           list(preds(fl.log1p, missing.idx),
                                ans.preds,
                                arcsin.preds,
                                raw.preds))

  # Test 3: scaling method.
  cat("  Running scaling tests.\n")

  fl.none <- fl.var0

  fl.ones <- flashier(log1p(samp), var.type = 0,
                      prior.type = "normal.mix",
                      fixed.factors = c(ones.factor(2), ones.factor(1)),
                      greedy.Kmax = K,
                      backfit.after = 2, final.backfit = FALSE,
                      verbose.lvl = 2L * verbose)

  scaled.samp <- samp * median(cell.sums) / rep(cell.sums, each = nrow(samp))
  fl.scale <- flashier(log1p(scaled.samp), var.type = 0,
                       prior.type = "normal.mix",
                       fixed.factors = ones.factor(2),
                       greedy.Kmax = K, backfit = "none",
                       verbose.lvl = 2L * verbose)
  scale.preds <- log1p((exp(preds(fl.scale, missing.idx)) - 1)
                       * cell.sums[missing.cols] / median(cell.sums))

  norm.df[i, ] <- all.mse(true.vals,
                          list(preds(fl.none, missing.idx),
                               preds(fl.ones, missing.idx),
                               scale.preds))

  # Test 4: prior type.
  cat("  Running prior type tests.\n")

  fl.normalmix <- fl.var0

  fl.nngenes <- flashier(log1p(samp), var.type = 0,
                         prior.type = c("nonnegative", "normal.mix"),
                         greedy.Kmax = K + 1, backfit = "none",
                         verbose.lvl = 2L * verbose)

  fl.nngenes.pm <- flashier(log1p(samp), var.type = 0,
                            prior.type = c("nonnegative", "normal.mix"),
                            ash.param = list(method = "fdr"),
                            greedy.Kmax = K + 1, backfit = "none",
                            verbose.lvl = 2L * verbose)

  fl.nncells <- flashier(log1p(samp), var.type = 0,
                         prior.type = c("normal.mix", "nonnegative"),
                         greedy.Kmax = K + 1, backfit = "none",
                         verbose.lvl = 2L * verbose)

  fl.nncells.pm <- flashier(log1p(samp), var.type = 0,
                            prior.type = c("normal.mix", "nonnegative"),
                            ash.param = list(method = "fdr"),
                            greedy.Kmax = K + 1, backfit = "none",
                            verbose.lvl = 2L * verbose)

  prior.df[i, ] <- all.mse(true.vals,
                           list(preds(fl.normalmix, missing.idx),
                                preds(fl.nngenes, missing.idx),
                                preds(fl.nngenes.pm, missing.idx),
                                preds(fl.nncells, missing.idx),
                                preds(fl.nncells.pm, missing.idx)))

  # Test 5: number of factors and backfit.
  cat("  Running backfitting tests.\n")

  fl.g <- fl.var0

  fl.b <- flashier(flash.init = fl.g, backfit = "only",
                   backfit.reltol = 10,
                   verbose.lvl = 3L * verbose)

  # Add K more factors, then K more.
  fl.g2 <- flashier(flash.init = fl.g,
                    greedy.Kmax = K, backfit = "none",
                    verbose.lvl = 2L * verbose)
  fl.b2 <- flashier(flash.init = fl.b,
                    greedy.Kmax = K, backfit = "final",
                    backfit.reltol = 10,
                    verbose.lvl = 3L * verbose)
  fl.g3 <- flashier(flash.init = fl.g2,
                    greedy.Kmax = K, backfit = "none",
                    verbose.lvl = 2L * verbose)
  fl.b3 <- flashier(flash.init = fl.b2,
                    greedy.Kmax = K, backfit = "final",
                    backfit.reltol = 10,
                    verbose.lvl = 3L * verbose)

  nfactors.df[i, ] <- all.mse(true.vals,
                              list(preds(fl.g, missing.idx),
                                   preds(fl.b, missing.idx),
                                   preds(fl.g2, missing.idx),
                                   preds(fl.b2, missing.idx),
                                   preds(fl.g3, missing.idx),
                                   preds(fl.b3, missing.idx)))
}

# Test 6: incremental addition of factors

set.seed(666)
data <- log1p(trachea[rowSums(trachea > 0) > 2, ])
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
