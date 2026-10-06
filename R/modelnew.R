## ============================================================================
## Q-VMP model construction
##
## A model contains:
##
##   Families: Exponential-family variational distributions. 
##
##   Tapes: AD functions for the log-partition function and Fisher information.
##
##   Expected Log-likelihood: AD functions for the expected log-likelihood.
##
## The model is constructed once and can then be fitted with different
## priors, starting values, or learning rates.
##
## Each family (beta.R, normal.R) is a list with
##   dim               length of the natural parameter eta (as stored by the engine)
##   theta_dim         length of the block theta_k it describes
##   A(eta)            log-partition (RTMB-compatible, so AD gives grad A and hess A)
##   mean_cov(eta)     list(mean = E theta_k, cov = Cov theta_k), used by the Taylor
##                     expansion (RTMB-compatible, it is taped inside E_q[loglik])
##   eta(...)          usual parameters -> natural parameter
##   par(eta)          natural parameter -> usual parameters
##   valid(eta)        is eta inside the natural parameter space?
##   template          any valid eta; only its structure is used, to build the AD tapes
## ============================================================================

# Need to create the import for this
library(RTMB)


#' Construct AD functions for an exponential family
#'
#' @param family Exponential-family specification.
#' @param eta Natural parameter used to construct the AD tape.
#'
#' @return A list containing the log-partition function, its gradient,
#' and the Fisher information AD functions.
#'
#' @keywords internal
make_family_tapes <- function(family, eta) {
  d <- family$dim
  
  # A(eta)
  A_tape <- MakeTape(family$A,eta)
  
  # grad A(eta)
  gradA_tape <- A_tape$jacfun()
  
  # Hessian A(eta) = Fisher information
  fisher_tape <- gradA_tape$jacfun()
  
  list(A = function(eta) {A_tape(eta)},
       gradA = function(eta) {as.vector(gradA_tape(eta))},
       fisher = function(eta) {matrix(fisher_tape(eta),d,d)}
  )
}


#' Construct a Q-VMP model
#'
#' @param loglik Log-likelihood function.
#' @param families List of exponential-family specifications.
#'
#' @return A Q-VMP model object.
#'
#' @keywords internal

qvmp_model <- function(loglik, families) {
  
  # prototype vector showing the tape what size and shape the input will be
  templates <- lapply(families,`[[`,"template")
  
  # class "qvmp_model": printing shows one line (print.R) instead of the tapes
  structure(
    list(families = families,
      # AD functions for the variational families
      tapes = lapply(seq_along(families),function(k) {make_family_tapes(families[[k]],templates[[k]])}),
      # Expected log-likelihood
      expected = make_expected_loglik(loglik,families,templates)
    ),
    class = "qvmp_model")
}


