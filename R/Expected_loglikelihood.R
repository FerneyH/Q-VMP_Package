## ============================================================================
## The expected log-likelihood approximation is:
##
##   E_q[log p(y | theta)] ~ log L(y | m) + 1/2 tr{ H(m) S },
##
## where:
##   m = E_q(theta),
##   S = Var_q(theta),
##   H(m) = Hessian of the log-likelihood evaluated at the variational mean.
## ============================================================================


#' Construct the expected log-likelihood
#'
#' @param loglik Log-likelihood function.
#' @param families List of exponential-family specifications.
#' @param eta Initial natural parameters used to build the AD tapes.
#'
#' @return A list containing the expected log-likelihood, its gradient,
#' and natural-parameter block indices.
#'
#' @keywords internal
make_expected_loglik <- function(loglik,
                                 families,
                                 eta){
  # number of q_k
  K <- length(families)
  
  # Indices for natural parameters
  eta_idx <- block_index(sapply(families,`[[`,"dim"))
  
  # Indices for model parameters, some q_k are multivariate
  theta_idx <- block_index(sapply(families,`[[`,"theta_dim"))
  
  # Total dimension
  n_theta <- sum(sapply(families,`[[`,"theta_dim"))
  
  # prototype vector showing the tape what size and shape the input will be, Using initial value for theta using the mean 
  theta0 <- unlist(lapply(seq_len(K),function(k){families[[k]]$mean_cov(eta[[k]])$mean}))
  
  # AD tape for the log-likelihood
  loglik_tape <- MakeTape(loglik,theta0)
  
  # Hessian of the log-likelihood
  hessian_tape <- loglik_tape$jacfun()$jacfun()
  
  # Taylor approximation
  expected_loglik <- function(p) {
    # Mean and covariance for each variational block
    mean_cov <- lapply(seq_len(K),function(k) {families[[k]]$mean_cov(p[eta_idx[[k]]])})  
    
    # Overall mean vector
    mean_theta <- do.call(c,lapply(mean_cov,`[[`,"mean"))
    
    # Hessian evaluated at the current mean
    H <- matrix(hessian_tape(mean_theta),n_theta,n_theta)
    
    # 1/2 tr(H S)
    correction <- 0
    for (k in seq_len(K)) {
      idx <- theta_idx[[k]]
      # Hessian
      H_k <- H[idx, idx]
      # Covariance under the variational distribution
      S_k <- mean_cov[[k]]$cov
      
      correction <- correction +sum(H_k * S_k)
    }
    
    # Expected loglikelihood approximation
    loglik_tape(mean_theta) + 0.5 * correction
  }
  
  # AD tape for the expected log-likelihood
  expected_tape <- MakeTape(expected_loglik,unlist(eta))
  
  # Return
  list(# Expected log-likelihood
       E = function(p) {expected_tape(p)},
       # Gradient of the expected loglikelihood 
       E_grad = function(p) {as.vector(expected_tape$jacobian(p))},
       # index
       eta_idx = eta_idx
      )
}
 
# helper that tells the code where each block of parameters 
# is located inside one long vector.
# For example:  if sizes <- c(3, 2), This means:
# Block 1 (variational family 1) has 3 parameters
# Block 2 (variational family 2) has 2 parameters

block_index <- function(sizes) {
  ends <- cumsum(sizes)
  lapply(seq_along(sizes), function(k) (ends[k] - sizes[k] + 1):ends[k])
}



