# Pseudo Inverse fitting method
fitting_method_pseudo_inverse <- function(design, outcome) {
  return(solve(t(design) %*% design, t(design) %*% outcome))
}

# BFGS fitting method
fitting_method_bfgs <- function(design, outcome, option) {
  beta_init <- solve(t(design) %*% design, t(design) %*% outcome)
  
  optim_res <- optim(
    beta_init, 
    fn = function(b) neg_log_likelihood(b, design, outcome),
    gr = function(b) neg_grad(b, design, outcome),
    method = "BFGS",
    control = if (!is.null(option$control)) option$control else list(fnscale = 1)
  )
  
  if (optim_res$convergence != 0) {
    warning("BFGS optimization did not converge. Consider checking your data, initial values, or optimization settings.")
  }
  
  return(optim_res$par)
}
