## ============================================================================
##
##   fit             short summary: mean and sd of each block
##   summary(fit)    credible intervals, covariance matrix of each block
##   coef(fit)       posterior means
##   vcov(fit)       posterior covariances, one matrix per block
##
## q is mean-field across blocks, q(theta) = prod_k q_k(theta_k), so
## correlations only exist inside a multivariate block, between blocks 
## they are zero by construction.
## ============================================================================


# Names of the coordinates of block k: "theta" or "beta[1]", ..., "beta[9]" for
# models like the logistic example
coord_names <- function(name, d) {
  if (d == 1) name else paste0(name, "[", 1:d, "]")
}


# Short description of a distribution in the family, e.g. "beta(alpha = 147.5, beta = 91.9)".
# Blocks with many parameters (a multivariate normal) only show their dimension.
family_label <- function(family, eta) {
  v <- unlist(family$par(eta))
  if (length(v) <= 4) {
    sprintf("%s(%s)", family$name, paste(names(v), "=", sapply(v, format, digits = 4), collapse = ", "))
  } else {
    sprintf("%s, dimension %d", family$name, family$theta_dim)
  }
}


# Mean and sd of each coordinate of block k
block_table <- function(fit, k) {
  mc <- fit$mean_cov[[k]]
  cbind(mean = mc$mean, sd = sqrt(diag(as.matrix(mc$cov))))
}


## print(fit)
print.qvmp_fit <- function(x, digits = 3, ...) {
  blocks <- names(x$families)

  cat(sprintf("Q-VMP fit: %d block%s, %s (rho = %g)\n", length(blocks),
              if (length(blocks) > 1) "s" else "",
              if (x$converged) sprintf("converged in %d iterations", x$iters)
              else sprintf("NOT converged after %d iterations", x$iters),
              x$rho))
  cat(sprintf("ELBO: %.2f\n", elbo(x)))

  for (k in seq_along(blocks)) {
    fam <- x$families[[k]]
    cat(sprintf("\nBlock %s: %s\n", blocks[k], family_label(fam, x$eta[[k]])))
    tab <- block_table(x, k)
    rownames(tab) <- coord_names(blocks[k], fam$theta_dim)
    print(round(tab, digits))
  }
  invisible(x)
}



## summary(fit): computes the tables, print.summary.qvmp_fit() shows them
summary.qvmp_fit <- function(object, level = 0.95, ...) {
  blocks <- names(object$families)
  probs  <- c((1 - level) / 2, (1 + level) / 2)

  per_block <- lapply(seq_along(blocks), function(k) {
    fam <- object$families[[k]]
    eta <- object$eta[[k]]

    # mean, sd and the credible interval from the quantiles of q_k
    tab <- cbind(block_table(object, k),
                 fam$quantile(eta, probs[1]),
                 fam$quantile(eta, probs[2]))
    colnames(tab)[3:4] <- paste0(format(100 * probs), "%")
    rownames(tab) <- coord_names(blocks[k], fam$theta_dim)

    # covariance matrix inside the block (only for multivariate blocks)
    S <- if (fam$theta_dim > 1) as.matrix(object$mean_cov[[k]]$cov)
    if (!is.null(S)) dimnames(S) <- list(rownames(tab), rownames(tab))

    list(posterior = family_label(fam, eta),
         prior     = family_label(fam, object$priors[[k]]),
         table     = tab,
         cov       = S)
  })
  names(per_block) <- blocks

  structure(list(blocks = per_block, level = level, elbo = elbo(object),
                 iters = object$iters, converged = object$converged, rho = object$rho,
                 Ftrace = object$Ftrace),
            class = "summary.qvmp_fit")
}


print.summary.qvmp_fit <- function(x, digits = 3, ...) {
  K <- length(x$blocks)

  cat(sprintf("Q-VMP fit: %d block%s, %s (rho = %g)\n", K, if (K > 1) "s" else "",
              if (x$converged) sprintf("converged in %d iterations", x$iters)
              else sprintf("NOT converged after %d iterations", x$iters),
              x$rho))
  cat(sprintf("ELBO: %.2f\n", x$elbo))
  if (!is.null(x$Ftrace)) {
    cat(sprintf("Free energy: %.2f -> %.2f (%s)\n", x$Ftrace[1], x$Ftrace[length(x$Ftrace)],
                if (all(diff(x$Ftrace) <= 1e-8)) "decreasing" else "not monotone"))
  }

  for (b in names(x$blocks)) {
    blk <- x$blocks[[b]]
    cat(sprintf("\nBlock %s\n", b))
    cat(sprintf("  posterior: %s\n", blk$posterior))
    cat(sprintf("  prior:     %s\n", blk$prior))
    print(round(blk$table, digits))

    # lower triangle of the covariance matrix, variances on the diagonal
    if (!is.null(blk$cov)) {
      # fixed decimals, enough to show `digits` significant digits of the largest entry
      decimals <- max(0, digits - 1 - floor(log10(max(abs(blk$cov)))))
      S <- formatC(blk$cov, format = "f", digits = decimals)
      S[upper.tri(S)] <- ""
      cat("\nCovariance:\n")
      print(S, quote = FALSE, right = TRUE)
    }
  }

  if (K > 1) cat("\nBlocks are independent under q (mean-field): no correlation between blocks.\n")
  invisible(x)
}



## coef(fit), vcov(fit)
# Posterior means, one named vector
coef.qvmp_fit <- function(object, ...) {
  blocks <- names(object$families)
  unlist(lapply(seq_along(blocks), function(k) {
    setNames(object$mean_cov[[k]]$mean, coord_names(blocks[k], object$families[[k]]$theta_dim))
  }))
}


# Posterior covariances, one matrix per block (no cross-block covariances under mean-field)
vcov.qvmp_fit <- function(object, ...) {
  blocks <- names(object$families)
  out <- lapply(seq_along(blocks), function(k) {
    S  <- as.matrix(object$mean_cov[[k]]$cov)
    nm <- coord_names(blocks[k], object$families[[k]]$theta_dim)
    dimnames(S) <- list(nm, nm)
    S
  })
  setNames(out, blocks)
}


