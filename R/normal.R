## ============================================================================
## Multivariate normal exponential family
##
## q(x) = N_p(mu, Sigma)
##
## Natural parameters:
##   eta1 = Sigma^{-1} mu
##   eta2 = -1/2 Sigma^{-1}
##
## Log-partition function, with the precision Lambda = Sigma^{-1} = -2 eta2:
##   A(eta) = 1/2 eta1' Lambda^{-1} eta1 - 1/2 log|Lambda|
##
## eta2 is symmetric, so only its lower triangle (vech) is stored:
##   eta = c(eta1, vech(eta2)),   dim = p + p(p + 1)/2
##
## p = 1 gives the univariate normal N(mu, s2)
## ============================================================================

#' Construct a multivariate Normal exponential family
#'
#' Defines a multivariate Normal distribution as an exponential family
#' for use as a variational factor.
#'
#' @param p Dimension of the multivariate Normal distribution.
#'
#' @return A list containing the exponential-family functions and
#' parameter information.
#'
#' @keywords internal

fam_mvnormal <- function(p) {
 
  # Positions of the lower triangle (diagonal included) in the p x p matrix
  lower <- which(lower.tri(matrix(0, p, p), diag = TRUE))

  # For every entry of the p x p matrix, its position in vech(eta2).
  index <- matrix(0, p, p)
  index[lower] <- seq_along(lower)
  index <- pmax(index, t(index))

  # vech -> symmetric p x p matrix 
  unvech <- function(v) {matrix(v[index], p, p)}

  # Natural parameter -> precision Lambda = -2 eta2
  precision <- function(eta) {unvech(-2 * eta[-(1:p)])}

  list(
    name = "mvnormal",

    # Number of natural parameters
    dim = p + length(lower),

    # Number of model parameters represented by this family
    theta_dim = p,

    # Log-partition function A(eta)
    A = function(eta) {
      eta1   <- eta[1:p]
      Lambda <- precision(eta)

      0.5 * sum(eta1 * solve(Lambda, eta1)) - 0.5 * determinant(Lambda, logarithm = TRUE)$modulus
    },

    # Mean and covariance
    mean_cov = function(eta) {
      eta1   <- eta[1:p]
      Lambda <- precision(eta)

      covariance <- solve(Lambda)
      mean       <- as.vector(covariance %*% eta1)

      list(mean = mean, cov = covariance)
    },

    # Usual parameters -> natural parameters
    eta = function(mu, Sigma) {
      Lambda <- solve(Sigma)

      eta1 <- Lambda %*% mu
      eta2 <- -0.5 * Lambda

      c(eta1, eta2[lower])
    },

    # Natural parameters -> usual parameters
    par = function(eta) {
      eta1   <- eta[1:p]
      Lambda <- precision(eta)

      Sigma <- solve(Lambda)
      mu    <- as.vector(Sigma %*% eta1)

      list(mu = mu, Sigma = Sigma)
    },

    # Marginal quantile of each coordinate, used for the credible intervals in summary()
    quantile = function(eta, prob) {
      Sigma <- solve(precision(eta))
      mu    <- as.vector(Sigma %*% eta[1:p])

      qnorm(prob, mu, sqrt(diag(Sigma)))
    },

    # Check whether eta is in the natural parameter space:
    # the precision must be positive definite
    valid = function(eta) {
      Lambda <- precision(eta)

      !inherits(try(chol(Lambda), silent = TRUE), "try-error")
    },

    # Valid starting value used to build the AD tapes
    # Corresponds to N(0, I) using usual parameters, here it is using natural parameters 
    template = c(rep(0, p), (-0.5 * diag(p))[lower])
  )
}
