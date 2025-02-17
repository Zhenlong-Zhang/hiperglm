#-------------------------------------------------------------
# -------------------Pseudo Inverse fitting method
#-------------------------------------------------------------
fitting_method_pseudo_inverse <- function(design, outcome) {
  return(solve(t(design) %*% design, t(design) %*% outcome))
}
#-------------------------------------------------------------
# --------------------bfgs
#-------------------------------------------------------------
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
#-------------------------------------------------------------
#-------------- -----newton
#-------------------------------------------------------------
fitting_method_newton <- function(design, outcome, option) {
  max_iter <- if (!is.null(option$max_iter)) option$max_iter else 20  # max
  tol <- if (!is.null(option$tol)) option$tol else 1e-6
  
  model_type <- if (!is.null(option$model)) option$model else "logistic"
  likelihood_funcs <- get_likelihood_functions(model_type)
  
  beta <- rep(0, ncol(design))
  last_update_norm <- Inf
  
  for (i in 1:max_iter) {
    grad <- likelihood_funcs$gradient(beta, design, outcome)
    H <- likelihood_funcs$hessian(beta, design, outcome)
    # in case it fails
    beta_update <- tryCatch(
      solve(H, grad),  
      error = function(e) {
        warning("Hessian is singular, using small step gradient update instead.")
        return(grad * 0.01) 
      }
    )
    beta <- beta - beta_update
    update_norm <- sum(abs(beta_update))  
    # warning
    if (update_norm < tol) {
      if (i < 3) {
        warning("Newton's method converged in very few iterations (", i, "). Consider checking data conditioning.")
      } else {
        message("Newton's method converged at iteration ", i)
      }
      break
    }
    if (i == max_iter) {
      warning("Newton's method reached max iterations (", max_iter, ") without full convergence.")
    }
    last_update_norm <- update_norm
  }
  return(beta)
}
