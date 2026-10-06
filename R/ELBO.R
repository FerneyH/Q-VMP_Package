## ============================================================================
## free energy ommiting constants 
## ELBO
## ============================================================================
#' 
#' @param tape AD functions for the variational family.
#' @param eta Variational natural parameter.
#' @param eta0 Prior natural parameter.
#'
#' @return Numeric value.
#'
#' @keywords internal
prior_variational_term <- function(tape,eta,eta0){
  sum((eta - eta0)*tape$gradA(eta))-tape$A(eta)
}


#' Computes the Taylor-approximated variational free energy,
#' up to the constant A(eta0).
#'
#' @param eta_list List of variational natural parameters.
#' @param tapes List of family AD functions.
#' @param priors List of prior natural parameters.
#' @param expected_loglik Expected log-likelihood function.
#'
#' @return Numeric free-energy value.
#'
#' @keywords internal
free_energy_objective <- function(eta_list,tapes,priors,expected_loglik) {
  # Prior contribution
  prior_variational_term_total <- sum(sapply(seq_along(tapes),function(k){prior_variational_term(tapes[[k]],eta_list[[k]],priors[[k]])}))
  # Negative expected log-likelihood
  prior_variational_term_total - expected_loglik(unlist(eta_list))
}


#' Compute the ELBO (negative free energy)
#'
#' The ELBO is the negative of the variational free energy.
#'
#' @param fit A fitted Q-VMP model.
#' 
#' @return Numeric ELBO value.
#'
#' @keywords internal
elbo <- function(fit) {-(fit$F + fit$A0)}