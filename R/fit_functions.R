source("R/likelihood_calculation.R")
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
  max_iter <- if (!is.null(option$max_iter)) option$max_iter else 50  
  epsilon_abs <- if (!is.null(option$epsilon_abs)) option$epsilon_abs else 1e-6  
  epsilon_rel <- if (!is.null(option$epsilon_rel)) option$epsilon_rel else 1e-6  
  epsilon_small <- 1e-8  

  model_type <- if (!is.null(option$model)) option$model else "logistic"
  likelihood_funcs <- get_likelihood_functions(model_type)

  beta <- rep(0, ncol(design))
  log_likelihood_old <- likelihood_funcs$log_likelihood(beta, design, outcome)

  for (i in 1:max_iter) {
    grad <- likelihood_funcs$gradient(beta, design, outcome)
    H <- likelihood_funcs$hessian(beta, design, outcome)

    beta_update <- tryCatch(
      solve(H, grad),  
      error = function(e) {
        warning("Hessian is singular, using small step gradient update instead.")
        return(grad * 0.01) 
      }
    )

    beta <- beta - beta_update 

    log_likelihood_new <- likelihood_funcs$log_likelihood(beta, design, outcome)
    change <- abs(log_likelihood_new - log_likelihood_old)  
    rel_change <- change / (abs(log_likelihood_old) + epsilon_small)  

    # converge check
    if (change < epsilon_abs || rel_change < epsilon_rel) {
      message("Newton's method converged at iteration ", i)
      return(beta)
    }

    log_likelihood_old <- log_likelihood_new  # update log-likelihood
  }

  # warning if not
  warning("Newton's method reached the maximum iteration limit of ", max_iter, " without full convergence.")
  return(beta)
}
