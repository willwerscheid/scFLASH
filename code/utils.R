# All functions assume that rows correspond to genes and columns to cells.

library(Matrix)
library(ggplot2)
library(flashier)


# Plotting functions for datasets -----------------------------------------------------

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
         title = "Library size of cells")
}

plot.meanexp <- function(data, bins = 50) {
  df <- data.frame(meanexp = rowMeans(data))
  ggplot(df, aes(x = meanexp)) +
    geom_histogram(bins = bins) +
    scale_x_log10() +
    labs(x = "mean expression",
         y = "gene count",
         title = "Mean expression of genes")
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


# Initial preprocessing ---------------------------------------------------------------

preprocess <- function(data,
                       min.nzcts,
                       max.libsize = NULL,
                       size.factors = NULL,
                       prescale = c("none", "genes", "cells", "both")) {
  prescale <- match.arg(prescale)

  retlist <- list()

  # Drop cells with too many transcripts.
  lib.size <- colSums(data)
  if (!is.null(max.libsize)) {
    dropped.cells <- which(lib.size > max.libsize)
    if (length(dropped.cells) > 0) {
      retlist$dropped.cells <- dropped.cells
      data <- data[, -dropped.cells]
      lib.size <- lib.size[-dropped.cells]
      if (!is.null(size.factors)) {
        size.factors <- size.factors[-dropped.cells]
      }
    }
  }
  n.cells <- ncol(data)

  if (is.null(size.factors)) {
    # Default size factors.
    size.factors <- lib.size / median(lib.size)
  } else if (length(size.factors) == 1) {
    size.factors <- rep(size.factors, n.cells)
  }

  # Drop genes with transcripts detected in too few cells.
  nz.cts <- rowSums(data > 0)
  dropped.genes <- which(nz.cts < min.nzcts)
  if (length(dropped.genes) > 0) {
    retlist$dropped.genes <- dropped.genes
    data <- data[-dropped.genes, ]
    nz.cts <- nz.cts[-dropped.genes]
  }
  n.genes <- nrow(data)

  retlist <- c(retlist, list(lib.size = lib.size,
                             size.factors = size.factors,
                             gene.sparsity = nz.cts / n.cells))

  # Scale and log-transform.
  data <- t(t(data) / size.factors)
  data <- log1p(data)

  # ELBO adjustment (using change-of-variables formula).
  elbo.adj <- -n.genes * sum(log(size.factors)) - sum(data)

  # Additional scaling to approximate a Kronecker variance structure.
  if (prescale != "none") {
    fl.kron <- flashier(data, var.type = c(1, 2), greedy.Kmax = 1,
                        final.nullchk = FALSE, verbose.lvl = 0)
    gene.sds <- fl.kron$residuals.sd[[1]]
    cell.sds <- fl.kron$residuals.sd[[2]]
    if (prescale %in% c("genes", "both")) {
      retlist$gene.prescaling.factors <- gene.sds
      data <- data / gene.sds
      elbo.adj <- elbo.adj - n.cells * sum(log(gene.sds))
    }
    if (prescale %in% c("cells", "both")) {
      retlist$cell.prescaling.factors <- cell.sds
      data <- t(t(data) / cell.sds)
      elbo.adj <- elbo.adj - n.genes * sum(log(cell.sds))
    }
  }

  retlist$elbo.adj <- elbo.adj
  retlist$data <- data

  return(retlist)
}

# Dataset-specific functions include defaults and return additional useful fields.

preprocess.droplet <- function(droplet,
                               min.nzcts = 10,
                               max.libsize = 30000,
                               size.factors = NULL) {
  processed <- preprocess(droplet, min.nzcts, max.libsize, size.factors)

  mouse <- sapply(strsplit(colnames(processed$data), "_"), `[`, 1)
  cell.type <- sapply(strsplit(colnames(processed$data), "_"), `[`, 3)

  processed$cell.type <- as.factor(cell.type)
  processed$mouse <- as.factor(mouse)

  processed$ionocyte.genes <- c("Ascl3", "Atp6v1c2", "Atp6v0d2", "Cftr", "Foxi1")

  processed$goblet.genes <- c("Gp2", "Tff1", "Tff2", "Muc5b", "Lman1l", "P2rx4", "Muc5ac",
                              "Dcpp1", "Dcpp2", "Dcpp3", "Lipf")
  processed$goblet1.genes <- c("Tff1", "Tff2", "Muc5b", "Lman1l", "P2rx4", "Muc5ac")
  processed$goblet2.genes <- c("Dcpp1", "Dcpp2", "Dcpp3", "Lipf")

  processed$all.goblet.genes <- c(processed$goblet1.genes, processed$goblet2.genes)
  processed$all.goblet.colors <- c(rep("blue", length(processed$goblet1.genes)),
                                   rep("red", length(processed$goblet2.genes)))

  processed$tuft.genes <- c("Il25", "Tslp", "Pou2f3", "Gnb3", "Gng13",
                            "Alox5ap", "Ptprc", "Spib", "Sox9")
  processed$tuft1.genes <- c("Gnb3", "Gng13", "Itpr3", "Plcb2", "Gnat3", "Ovol3",
                             "Commd1", "Cited2", "Atp1b1", "Fxyd6", "Pou2f3")
  processed$tuft2.genes <- c("Alox5ap", "Ptprc", "Spib", "Sox9", "Mgst3", "Gpx2",
                             "Ly6e", "S100a11", "Cst3", "Cd24a", "B2m", "Dclk1",
                             "Sdc4", "Il13ra1")

  processed$all.tuft.genes <- with(processed, union(c(tuft1.genes, tuft2.genes),
                                                    tuft.genes))
  processed$all.tuft.colors <- with(processed, {
    c(rep("blue", length(tuft1.genes)),
      rep("red", length(tuft2.genes)),
      rep("black", length(setdiff(tuft.genes, c(tuft1.genes, tuft2.genes)))))
  })

  processed$hillock.genes <- c("Krt13", "Krt4", "Ecm1", "S100a11", "Cldn3", "Lgals3",
                               "Anxa1", "S100a6", "Upk3bl", "Aqp5", "Anxa2", "Crip1",
                               "Gsto1", "Tppp3")

  return(processed)
}


# Plotting functions for flashier fits ------------------------------------------------

plot.factors <- function(res, cell.types,
                         kset = NULL, max.pt.size = 2, title = NULL) {
  # Sort loadings according to proportion of variance explained.
  if (is.null(kset)) {
    kset <- setdiff(order(res$fl$pve, decreasing = TRUE), which(res$fl$pve == 0))
  }

  if (is.null(res$cell.prescaling.factors)) {
    res$cell.prescaling.factors <- 1
  }

  # Re-normalize loadings so that factors are equally spread out.
  LL <- res$fl$loadings.pm[[2]][, kset, drop = FALSE]
  LL <- LL * res$cell.prescaling.factors
  LL <- t(t(LL) / apply(abs(LL), 2, max))

  # To make it easier to compare factors, flip them to make the largest
  #   loadings positive.
  flip <- 2 * (colSums(LL > 0.75) > colSums(LL < -0.75)) - 1
  LL <- t(t(LL) * flip)

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

get.orig.k <- function(fl, k) {
  return(order(fl$pve, decreasing = TRUE)[k])
}

plot.one.factor <- function(fl, k, notable.genes, gene.prescaling.factors = 1,
                            title = NULL, label.size = 8, top.n = 100,
                            invert = FALSE, gene.colors = NULL) {
  df <- data.frame(gene = rownames(fl$loadings.pm[[1]]),
                   loading = fl$loadings.pm[[1]][, k] * gene.prescaling.factors)
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

plot.factor.subpops <- function(fl, k, subpop1.genes, subpop2.genes,
                                subpop1.label, subpop2.label, title = NULL) {
  all.genes <- c(subpop1.genes, subpop2.genes)
  df <- data.frame(gene = all.genes,
                   type = c(rep(subpop1.label, length(subpop1.genes)),
                            rep(subpop2.label, length(subpop2.genes))),
                   loading = fl$loadings.pm[[1]][all.genes, k])
  ggplot(df, aes(x = reorder(gene, -loading), y = loading, fill = type)) +
    geom_bar(stat = "identity") +
    labs(x = NULL, title = title) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}


# Calculating and plotting p-values ---------------------------------------------------

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

calc.p.vals <- function(orig.data, processed, fl, var.type) {
  set.seed(666)

  # Get the original (untransformed) data.
  if (!is.null(processed$dropped.genes)) {
    orig.data <- orig.data[-processed$dropped.genes, ]
  }
  if (!is.null(processed$dropped.cells)) {
    orig.data <- orig.data[, -processed$dropped.cells]
  }

  # Calculate mu and sigma.
  mu <- fitted(fl)
  if (!is.null(processed$gene.prescaling.factors)) {
    mu <- mu * processed$gene.prescaling.factors
  }
  if (!is.null(processed$cell.prescaling.factors)) {
    mu <- t(t(mu) * processed$cell.prescaling.factors)
  }

  n.genes <- nrow(processed$data)
  n.cells <- ncol(processed$data)
  if (length(var.type) > 1) {
    gene.sds <- fl$residuals.sd[[1]]
    cell.sds <- fl$residuals.sd[[2]]
  } else if (var.type == 0) {
    gene.sds <- rep(fl$residuals.sd, n.genes)
    cell.sds <- rep(1, n.cells)
  } else if (var.type == 1) {
    gene.sds <- fl$residuals.sd
    cell.sds <- rep(1, n.cells)
  } else if (var.type == 2) {
    gene.sds <- rep(1, n.genes)
    cell.sds <- fl$residuals.sd
  }

  if (!is.null(processed$gene.prescaling.factors))
    gene.sds <- gene.sds * processed$gene.prescaling.factors
  if (!is.null(processed$cell.prescaling.factors))
    cell.sds <- cell.sds * processed$cell.prescaling.factors

  sigma <- outer(gene.sds, cell.sds)

  # Calculate the moments of the fitted model, adjusting the model so that values
  #   less than zero are collapsed to a point mass at zero.
  lambda <- processed$size.factors
  EY.over.l <- calc.EY.over.lambda(mu, sigma)
  EY2 <- t(t(calc.EY2.over.lambda2(mu, sigma, EY.over.l)) * lambda^2)
  EY <- t(t(EY.over.l) * lambda)
  VarY <- EY2 - EY^2

  pois.data <- orig.data[EY >= VarY]
  pois.EY <- EY[EY >= VarY]
  pois.c <- runif(length(pois.data))
  pois.pvals <- (pois.c * ppois(pois.data, pois.EY)
                 + (1 - pois.c) * ppois(pois.data - 1, pois.EY))
  pois.llik <- sum(dpois(pois.data, pois.EY, log = TRUE))

  NB.data <- orig.data[EY < VarY]
  NB.p <- EY[EY < VarY] / VarY[EY < VarY]
  NB.r <- EY[EY < VarY] * NB.p / (1 - NB.p)
  NB.c <- runif(length(NB.data))
  NB.pvals <- (NB.c * pnbinom(NB.data, NB.r, NB.p)
               + (1 - NB.c) * pnbinom(NB.data - 1, NB.r, NB.p))
  NB.llik <- sum(dnbinom(NB.data, NB.r, NB.p, log = TRUE))

  p.vals  <- c(pois.pvals, NB.pvals)
  p.vals  <- ceiling(100 * p.vals)
  p.table <- table(p.vals)

  probs <- p.table / sum(p.table)
  KL.div  <- sum(probs * log(100 * probs))
  llik <- pois.llik + NB.llik

  return(list(table = p.table, KL.divergence = KL.div, llik = llik))
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
    scale_y_continuous(breaks = seq(0, max(df$Freq), by = 0.002)) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}


# Aligning and comparing factors from different fits ----------------------------------

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


# Not used ----------------------------------------------------------------------------

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
