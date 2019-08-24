# All functions assume that rows correspond to genes and columns to cells.

library(Matrix)
library(ggplot2)
library(flashier)

# Plotting functions for datasets ---------------------------------------------

plot.category <- function(category, title) {
  df <- data.frame(category)
  ggplot(df, aes(x = category)) +
    geom_bar() +
    labs(x = NULL, title = title)
}

plot.libsize <- function(data, bins = 50) {
  df <- data.frame(libsize = colSums(data))
  ggplot(df, aes(x = libsize)) +
    geom_histogram(bins = bins) +
    scale_x_log10() +
    labs(x = "library size",
         y = "cell count",
         title = "Library Size of Cells")
}

plot.meanexp <- function(data, bins = 50) {
  df <- data.frame(meanexp = rowMeans(data))
  ggplot(df, aes(x = meanexp)) +
    geom_histogram(bins = bins) +
    scale_x_log10() +
    labs(x = "mean expression",
         y = "gene count",
         title = "Mean Expression of Genes")
}

plot.gene <- function(data, gene, bins = 30) {
  df <- data.frame(expression = data[gene, ])
  ggplot(df, aes(x = expression)) +
    geom_histogram(bins = bins) +
    scale_x_continuous(trans = "log1p",
                       breaks = c(0, 10^(0:floor(log10(max(data[gene, ])))))) +
    labs(x = "number of transcripts detected",
         y = "cell count",
         title = gene)
}

# Initial preprocessing -------------------------------------------------------

preprocess <- function(data, min.nzcts, max.libsize = NULL, scale.factors = NULL) {
  # Drop cells with too many transcripts.
  lib.size <- colSums(data)
  if (!is.null(max.libsize)) {
    data     <- data[,  lib.size <= max.libsize]
    if (!is.null(scale.factors)) {
      scale.factors <- scale.factors[lib.size <= max.libsize]
    }
    lib.size <- lib.size[lib.size <= max.libsize]
  }

  if (is.null(scale.factors)) {
    scale.factors <- lib.size / median(lib.size)
  }

  n.genes <- nrow(data)

  # Drop genes with transcripts detected in too few cells.
  nz.cts <- rowSums(data > 0)
  data   <- data[nz.cts >= min.nzcts, ]
  nz.cts <- nz.cts[nz.cts >= min.nzcts]

  gene.sparsity <- nz.cts / n.genes

  # Default additional pre-processing: normalize and transform.
  data <- t(t(data) / scale.factors)
  data <- log1p(data)

  # ELBO adjustment (using change-of-variables formula).
  elbo.adj <- -n.genes * sum(log(scale.factors)) - sum(data)

  return(list(data = data,
              scale.factors = scale.factors,
              elbo.adj = elbo.adj,
              gene.sparsity = gene.sparsity))
}

# Dataset-specific functions include some additional useful fields.

preprocess.droplet <- function(droplet,
                               min.nzcts = 10,
                               max.libsize = 30000,
                               scale.factors = NULL) {
  processed <- preprocess(droplet, min.nzcts, max.libsize, scale.factors)

  mouse <- sapply(strsplit(colnames(processed$data), "_"), `[`, 1)
  cell.type <- sapply(strsplit(colnames(processed$data), "_"), `[`, 3)

  processed$cell.type <- as.factor(cell.type)
  processed$mouse <- as.factor(mouse)

  return(processed)
}

# Parallel analysis. Sort of. -------------------------------------------------
est.baseline.pve <- function(data, n.trials, q = 0.9, seeds = 1:n.trials, ...) {
  pve <- rep(NA, n.trials)
  for (i in 1:n.trials) {
    set.seed(seeds[i])
    cat("Seed:", seeds[i], "\n")
    rand.data <- t(Matrix(apply(data, 1, FUN = sample)))
    fl <- flashier(rand.data, greedy.Kmax = 2, verbose.lvl = 0, ...)
    if (length(fl$pve) > 1) {
      pve[i] <- fl$pve[2]
    } else {
      pve[i] <- 0
    }
  }
  return(list(baseline.pve = quantile(pve, q), all.res = pve))
}

# 1. "Factors should be easily interpretable." -------------------------------------

plot.one.factor <- function(fl, k, notable.genes, title = NULL,
                            label.size = 8, top.n = 100, invert = FALSE,
                            gene.colors = NULL) {
  df <- data.frame(gene = rownames(fl$loadings.pm[[1]]),
                   loading = fl$loadings.pm[[1]][, k])
  if (invert) {
    df$loading <- -df$loading
  }
  df <- df[order(df$loading, decreasing = TRUE), ]
  df <- df[1:top.n, ]
  df$order <- 1:nrow(df)
  plt <- ggplot(df, aes(x = order, y = loading)) +
    geom_point() +
    scale_x_continuous(breaks = df$order[which(df$gene %in% notable.genes)],
                       labels = df$gene[which(df$gene %in% notable.genes)]) +
    labs(x = NULL, title = title) +
    lims(y = c(0, NA)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = label.size))
  for (i in 1:length(notable.genes)) {
    gene <- notable.genes[i]
    which.gene <- which(df$gene == gene)
    if (length(which.gene) == 1) {
    if (!is.null(gene.colors)) {
      the.color = gene.colors[i]
    } else {
      the.color = "red"
    }
      plt <- plt +
        geom_segment(x = which.gene, xend = which.gene,
                     y = 0, yend = df$loading[which.gene], color = the.color)
    }
  }
  plot(plt)
}

plot.factors <- function(fl, cell.types, kset = NULL, max.pt.size = 2, title = NULL) {
  if (is.null(kset)) {
    kset <- 1:fl$n.factors
  }

  # Re-normalize loadings so that factors are equally spread out.
  LL <- fl$loadings.pm[[2]][, kset, drop = FALSE]
  LL <- t(t(LL) / apply(abs(LL), 2, max))

  # Make the size of the point depend on how many of that type there are.
  sizes <- max.pt.size / sqrt(table(cell.types) / min(table(cell.types)))

  df <- reshape2::melt(LL, value.name = "loading")
  df$cell.type <- rep(as.factor(cell.types), length(kset))
  ggplot(df, aes(x = Var2, y = loading, color = cell.type)) +
    geom_jitter(position = position_jitter(0.45),
                size = rep(sizes[cell.types], length(kset))) +
    labs(title = title, x = NULL) +
    lims(y = c(-1.05, 1.05))
}

compare.factors <- function(fl1, fl2, match.n = 1, min.cor = 0,
                          incl.sparsity = TRUE, coherence.mat = NULL) {
  cormat <- crossprod(fl1$loadings.pm[[match.n]], fl2$loadings.pm[[match.n]])

  k.matches  <- rep(NA, nrow(cormat))
  match.cors <- rep(NA, nrow(cormat))

  cor.order <- order(abs(cormat), decreasing = TRUE)
  for (entry in cor.order) {
    entry.k1 <- row(cormat)[entry]
    entry.k2 <- col(cormat)[entry]
    # If neither factor has been matched yet, match them.
    if (is.na(k.matches[entry.k1]) && !(entry.k2 %in% k.matches)) {
      k.matches[entry.k1]  <- entry.k2
      match.cors[entry.k1] <- cormat[entry]
    }
  }

  res <- data.frame(fl1.k = 1:nrow(cormat), fl2.k = k.matches, cor = match.cors)
  res <- subset(res, abs(cor) > min.cor)

  if (incl.sparsity) {
    res$fl1.sparsity.genes <- calc.sparsity(fl1, 1)[res$fl1.k]
    res$fl2.sparsity.genes <- calc.sparsity(fl2, 1)[res$fl2.k]
    res$fl1.sparsity.cells <- calc.sparsity(fl1, 2)[res$fl1.k]
    res$fl2.sparsity.cells <- calc.sparsity(fl2, 2)[res$fl2.k]
  }

  if (!is.null(coherence.mat)) {
    res$fl1.coherence <- calc.coherence(fl1, coherence.mat, 1)[res$fl1.k]
    res$fl2.coherence <- calc.coherence(fl2, coherence.mat, 1)[res$fl2.k]
  }

  return(res)
}

calc.sparsity <- function(fl, n = 1, lfsr.thresh = 0.01) {
  return(colSums(fl$loadings.lfsr[[n]] > lfsr.thresh) / nrow(fl$loadings.pm[[n]]))
}

calc.coherence <- function(fl, coherence.mat, n = 1, lfsr.thresh = 0.01) {
  gene.symbols  <- rownames(fl$loadings.pm[[n]])
  avail.symbols <- rownames(coherence.mat)

  fl.lfsr <- fl$loadings.lfsr[[n]][gene.symbols %in% avail.symbols, ]

  gene.symbols <- gene.symbols[gene.symbols %in% avail.symbols]

  fl.signif <- apply(fl.lfsr, 2, function(x) which(x < lfsr.thresh))

  res <- rep(NA, length(fl.signif))
  for (k in 1:length(fl.signif)) {
    idx <- which(avail.symbols %in% gene.symbols[fl.signif[[k]]])
    submat <- coherence.mat[idx, idx]
    diag(submat) <- 0 # Don't count diagonal entries.
    res[k] <- sum(submat) / (length(submat) - length(idx))
  }

  return(res)
}

build.coherence.mat <- function(pw.mat, data) {
  gene.set <- rownames(data)
  pw.mat <- pw.mat[rownames(pw.mat) %in% gene.set, ]

  # Remove pathways with no genes (or one gene).
  pw.mat <- pw.mat[, colSums(pw.mat) > 1]

  n <- nrow(pw.mat)
  coherence.mat <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    idx <- which(pw.mat[i, ] == 1)
    if (length(idx) > 0) {
      cooccurs <- rowSums(pw.mat[, idx, drop = FALSE])
      coherence.mat[i, ] <- cooccurs / length(idx)
    } else {
      coherence.mat[i, ] <- 0
    }
  }

  # Convert to a sparse matrix and symmetrize.
  coherence.mat <- Matrix(coherence.mat)
  coherence.mat <- (coherence.mat + t(coherence.mat)) / 2

  rownames(coherence.mat) <- rownames(pw.mat)
  colnames(coherence.mat) <- rownames(pw.mat)

  return(coherence.mat)
}

plot.factor.comparison <- function(df,
                                   type = c("gene.sparsity", "cell.sparsity", "coherence")) {
  type <- match.arg(type)
  the.aes  <- switch(type,
                     gene.sparsity = aes(x = fl1.sparsity.genes, y = fl2.sparsity.genes),
                     cell.sparsity = aes(x = fl1.sparsity.cells, y = fl2.sparsity.cells),
                     coherence     = aes(x = fl1.coherence, y = fl2.coherence))

  if (type == "coherence") {
    baseline <- mean(df$fl1.coherence[which(df$fl1.k == 1)],
                     df$fl2.coherence[which(df$fl2.k == 1)])
  }
  the.lims <- switch(type,
                     gene.sparsity = c(0, 1),
                     cell.sparsity = c(0, 1),
                     coherence     = c(0.9 * baseline,
                                       max(c(df$fl1.coherence, df$fl2.coherence))))

  plt <- ggplot(subset(df, fl1.k != 1), the.aes) +
    geom_point() +
    geom_abline(slope = 1, linetype = "dashed") +
    lims(x = the.lims, y = the.lims)
  if (type == "coherence") {
    plt <- plt +
      geom_hline(yintercept = baseline, linetype = "dotted") +
      geom_vline(xintercept = baseline, linetype = "dotted")
  }
  plot(plt)
}

# 2. "Factors should be able to successfully decompose samples we haven’t seen." ---

pve.test <- function(data, flashier.call, holdout.set, Kmax) {
  n.trials <- nrow(holdout.set)

  training.pve <- matrix(NA, nrow = n.trials, ncol = Kmax)
  holdout.pve  <- matrix(NA, nrow = n.trials, ncol = Kmax)

  for (i in 1:n.trials) {
    holdout.cells <- holdout.set[i, ]
    fl <- flashier.call(data[, -holdout.cells])
    training.pve[i, ] <- fl$pve
    holdout.pve[i, ]  <- calc.holdout.pve(fl, data[, holdout.cells])
  }
  return(data.frame(training.pve = training.pve, holdout.pve = holdout.pve))
}

# Not PVE strictly speaking since factors aren't orthogonal... but close.
calc.holdout.pve <- function(fl, holdout.data) {
  LL <- fl$loadings.pm[[1]]
  FF <- solve(crossprod(LL), crossprod(LL, holdout.data))
  return(rowSums(FF^2) / sum(holdout.data^2))
}

plot.pve.test.res <- function(res, baseline) {
  df     <- reshape2::melt(res, value.name = "pve")
  df$k   <- as.numeric(sapply(strsplit(as.character(df$variable), "[.]"), `[`, 3))
  df$set <- as.factor(sapply(strsplit(as.character(df$variable), "[.]"), `[`, 1))
  ggplot(df, aes(x = k, y = pve, color = set)) +
    geom_jitter(position = position_jitter(0.2)) +
    scale_y_log10() +
    geom_hline(yintercept = baseline, linetype = "dashed")
}

# 3. "The distribution of residuals should appear reasonable." ---------------------
zscore.resids <- function(data, fl) {
  if (!is.null(fl$residuals.sd) && fl$flash.fit$est.tau.dim < 2) {
    z <- as.vector(data - fitted(fl)) / fl$residuals.sd
  } else if (!is.null(fl$residuals.sd)) {
    z <- as.vector(t(data - fitted(fl))) / fl$residuals.sd
  } else {
    z <- (data - fitted(fl)) * sqrt(fl$flash.fit$tau[[1]])
    z <- as.vector(t(z) * sqrt(fl$flash.fit$tau[[2]]))
  }
  return(z)
}

# First and second moments of a truncated lognormal distribution with support
#   on [1, \infty].
etrunclnorm <- function(mu, sigma) {
  return(exp(mu + sigma^2 / 2) * pnorm((mu + sigma^2) / sigma)
         / pnorm(mu / sigma))
}

e2trunclnorm <- function(mu, sigma) {
  return(exp(2 * mu + 2 * sigma^2) * pnorm((mu + 2 * sigma^2) / sigma)
         / pnorm(mu / sigma))
}

# First and second moments of Y / lambda, where Y is a lognormal distribution
#   that has been shifted left by 1 and then constrained to have nonnegative
#   support (by moving all mass on negative values to a point mass at zero).
calc.EY.over.lambda <- function(mu, sigma) {
  return(exp(mu + sigma^2 / 2) * pnorm((mu + sigma^2) / sigma) - pnorm(mu / sigma))
}

calc.EY2.over.lambda2 <- function(mu, sigma, EY.over.lambda) {
  return(exp(2 * mu + 2 * sigma^2) * pnorm((mu + 2 * sigma^2) / sigma)
         - pnorm(mu / sigma) - 2 * EY.over.lambda)
}

# Currently only works for constant and gene-wise variance structures.
calc.p.vals <- function(data, sf, mu, sigma) {
  # Get the original (untransformed) data.
  data <- t(t(exp(data) - 1) / sf)

  # Calculate the moments of the fitted model, adjusting the model so that values
  #   less than zero are collapsed to a point mass at zero.
  EY.over.l <- calc.EY.over.lambda(mu, sigma)
  EY2  <- t(sf^2 * t(calc.EY2.over.lambda2(mu, sigma, EY.over.l)))
  EY   <- t(sf * t(EY.over.l))
  VarY <- EY2 - EY^2

  pois.data  <- data[EY > VarY]
  pois.EY    <- EY[EY > VarY]
  pois.c     <- runif(length(pois.data))
  pois.pvals <- (pois.c * ppois(pois.data, pois.EY)
                 + (1 - pois.c) * ppois(pois.data - 1, pois.EY))

  NB.data  <- data[EY < VarY]
  NB.p     <- EY[EY < VarY] / VarY[EY < VarY]
  NB.r     <- EY[EY < VarY] * NB.p / (1 - NB.p)
  NB.c     <- runif(length(NB.data))
  NB.pvals <- (NB.c * pnbinom(NB.data, NB.r, NB.p)
               + (1 - NB.c) * pnbinom(NB.data - 1, NB.r, NB.p))

  p.vals  <- c(pois.pvals, NB.pvals)
  p.vals  <- ceiling(100 * p.vals)
  p.table <- table(p.vals)

  probs <- p.table / sum(p.table)
  KL.div  <- sum(probs * log(100 * probs))

  return(list(table = p.table, KL.divergence = KL.div))
}

plot.p.vals <- function(p.vals) {
  df <- data.frame(p.vals$table)
  df$p.vals <- as.numeric(df$p.vals) / 100 - .005
  df$Freq <- df$Freq / sum(df$Freq)
  title <- paste("KL divergence relative to uniform:", round(p.vals$KL.divergence, 3))
  ggplot(df, aes(x = p.vals, y = Freq)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.01, linetype = "dashed") +
    labs(x = "p-value", y = NULL, title = title) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
