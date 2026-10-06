## ============================================================================
## One Bayesian Learning Rule step for an exponential-family block
##
## eta <- (1 - rho) eta + rho * natgrad_eta E_q[log p]
##
## With E_q[log p] = E_q[loglik] + E_q[log pi] and pi in the same family as q,
## equivalent to:
##
##  eta <- (1 - rho) eta + rho (eta0 + (Fisher)^{-1}*grad_eta E_q[loglik])
##
## ============================================================================
#' One Bayesian Learning Rule step
#'
#' Performs one natural-gradient Bayesian Learning Rule update for
#' an exponential-family variational block.
#'
#' @param tapes AD functions for the variational family.
#' @param eta Current natural parameter.
#' @param eta0 Prior natural parameter.
#' @param grad_eta Gradient of the expected log-likelihood with
#' respect to the natural parameter.
#' @param rho Learning rate.
#'
#' @return Updated natural parameter.
#'
#' @keywords internal


blr_step <- function(tapes, eta, eta0, grad_eta, rho) {
  ## Hessian/Fisher by AD
  # Symmetric positive definite: fisher = R'R (Cholesky parameterization)
  fisher<-tapes$fisher(eta)
  R <- chol(fisher)
  
  ## Inverse Fisher * Gradient: solve R'z = grad_eta, then R x = z
  natural_gradient <- backsolve(R, backsolve(R, grad_eta, transpose = TRUE))
  
  ## Updates
  (1 - rho) * eta + rho * (eta0 + natural_gradient)
}