## ============================================================================
## qvmp_fit(): run the BLR until the variational parameters stop moving
##
## The model (AD tapes) is built once with qvmp_model() (modelnew.R).
##
## Inputs
##   model     a Q-VMP model from qvmp_model()
##   priors    list, the prior natural parameter eta0 for each block (same family as q)
##   init      list, starting natural parameter for each block (default: the prior)
##   rho       fixed learning rate
##   trace     keep F at every iteration (fit$Ftrace)
## ============================================================================

qvmp_fit <- function(model, priors, init = priors, rho = 0.8, max_iters = 1000, tol = 1e-6,
                     trace = FALSE) {
  families <- model$families
  tapes    <- model$tapes
  ex       <- model$expected
  K        <- length(families)
  eta      <- init

  Ftrace <- if (trace) free_energy_objective(eta, tapes, priors, ex$E)
  for (iter in 1:max_iters) {
    old  <- unlist(eta)
    # grad_eta E_q[loglik], all blocks
    grad <- ex$E_grad(old)
    for (k in 1:K) {
      eta[[k]] <- blr_step(tapes[[k]], eta[[k]], priors[[k]], grad[ex$eta_idx[[k]]], rho)
      if (!families[[k]]$valid(eta[[k]])) stop(sprintf("block %d left its parameter space", k))
    }
    if (trace) Ftrace <- c(Ftrace, free_energy_objective(eta, tapes, priors, ex$E))
    if (max(abs(unlist(eta) - old)) < tol) break
  }
  # class "qvmp_fit": printing shows a short summary (print.R), all pieces stay in the list
  structure(
    list(eta = eta,
         par = setNames(lapply(1:K, function(k) families[[k]]$par(eta[[k]])), names(families)),
         mean_cov = setNames(lapply(1:K, function(k) families[[k]]$mean_cov(eta[[k]])), names(families)),
         # free energy up to the constants A(eta0)
         F = free_energy_objective(eta, tapes, priors, ex$E), Ftrace = Ftrace,
         A0 = sum(sapply(1:K, function(k) tapes[[k]]$A(priors[[k]]))),
         iters = iter, converged = iter < max_iters, rho = rho,
         families = families, priors = priors),
    class = "qvmp_fit")
}
