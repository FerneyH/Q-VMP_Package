## ============================================================================
## Beta exponential family
##
## q(x) = Beta(alpha, beta)
##
## Natural parameters:
##   eta = (alpha - 1, beta - 1)
##
## Mean and variance:
##   E[X]   = alpha / (alpha + beta)
##   Var[X] = alpha * beta /((alpha + beta)^2 * (alpha + beta + 1))
## ============================================================================
#' Defines the Beta distribution as an exponential family for use
#' as a variational factor.
#' @return A list containing the exponential-family functions and
#' parameter information.
#'
#' @keywords internal

fam_beta <- list(
  
  name = "beta",
  
  # Number of natural parameters
  dim = 2,
  
  # Number of model parameters represented by this family
  theta_dim = 1,
  
  # Log-partition function A(eta)
  A = function(eta) {
    alpha <- eta[1] + 1
    beta  <- eta[2] + 1
    
    lgamma(alpha) + lgamma(beta) - lgamma(alpha + beta)
  },
  
  
  # Mean and covariance 
  mean_cov = function(eta) {
    alpha <- eta[1] + 1
    beta  <- eta[2] + 1
    
    mean <- alpha / (alpha + beta)
    
    variance <-alpha * beta /((alpha + beta)^2 * (alpha + beta + 1))
    
    list(mean = mean, cov  = variance)
  },
  
  # Usual parameters -> natural parameters
  eta = function(alpha, beta) {
    c(alpha - 1,beta - 1)
  },
  
  
  # Natural parameters -> usual parameters
  par = function(eta) {
    c(alpha = eta[1] + 1,beta  = eta[2] + 1)
  },
  
  
  # Quantile of q, used for the credible intervals in summary()
  quantile = function(eta, prob) {
    qbeta(prob, eta[1] + 1, eta[2] + 1)
  },


  # Check whether eta is in the natural parameter space
  valid = function(eta) {
    all(eta > -1)
  },
  
  
  # Valid starting value used to build the AD tapes
  # Corresponds to Beta(1, 1)
  template = c(0, 0)
)