# Pseudo Inverse fitting method
fitting_method_pseudo_inverse <- function(design, outcome) {
  return(solve(t(design) %*% design, t(design) %*% outcome))
}

fitting_method_bfgs <- function(design, outcome, option) {
  model_type <- if (!is.null(option$model)) option$model else "linear"
  
  # import the nll cal and neg grad based on model type
  likelihood_funcs <- get_likelihood_functions(model_type)
  
  if (model_type == "logistic") {
    beta_init <- rep(0, ncol(design))  
  } else if (model_type == "linear") {
    XtX <- t(design) %*% design
    Xty <- t(design) %*% outcome
    beta_init <- tryCatch(
      solve(XtX, Xty),  
      error = function(e) {
        warning("Matrix inversion failed. Using zeros as initial beta.")
        return(rep(0, ncol(design)))  
      }
    )
  } else {
    stop("Unsupported model type: ", model_type)
  }
  
  # bfgs
  optim_res <- optim(
    beta_init, 
    fn = function(b) likelihood_funcs$log_likelihood(b, design, outcome),
    gr = function(b) likelihood_funcs$gradient(b, design, outcome),
    method = "BFGS",
    control = if (!is.null(option$control)) option$control else list(fnscale = 1, maxit = 1000)
  )
  
  if (optim_res$convergence != 0) {
    warning("BFGS optimization did not converge. Consider checking your data, initial values, or optimization settings.")
  }
  
  return(optim_res$par)
}