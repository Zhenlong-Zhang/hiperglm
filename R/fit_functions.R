#' Internal Fit Function for hiperglm
#'
#' @section Pseudo Inverse Method:
#'
#' These functions implement the pseudo inverse fitting method.
#'
#' @rdname fitting_methods
#' @keywords internal
fitting_method_pseudo_inverse <- function(design, outcome) {
  return(solve(t(design) %*% design, t(design) %*% outcome))
}


#' @section BFGS Fitting Method:
#'
#' Implements the BFGS optimization method for fitting.
#'
#' @rdname fitting_methods
#' @keywords internal
fitting_method_bfgs <- function(design, outcome, option) {
  model_type <- if (!is.null(option$model)) option$model else "linear"
  neg_likelihood_funcs <- get_neg_likelihood_functions(model_type)
  
  if (model_type == "logistic") {
    p <- mean(outcome)
    if (p > 0 && p < 1) {
      init_intercept <- log(p / (1 - p))
      beta_init <- rep(0, ncol(design))
      beta_init[1] <- init_intercept
    } else {
      beta_init <- rep(0, ncol(design))
    }
  } else if (model_type == "linear") {
    beta_init <- rep(0, ncol(design))
  }
  
  optim_res <- optim(
    beta_init, 
    fn = function(b) neg_likelihood_funcs$neg_log_likelihood(b, design, outcome),
    gr = function(b) neg_likelihood_funcs$neg_gradient(b, design, outcome),
    method = "BFGS",
    control = if (!is.null(option$control)) option$control else list(fnscale = 1, maxit = 1000)
  )
  
  if (optim_res$convergence != 0) {
    warning("BFGS optimization did not converge. Consider checking your data, initial values, or optimization settings.")
  }
  return(optim_res$par)
}

#' @section Newton Fitting Method:
#'
#' Implements Newton's method for logistic regression.
#'
#' @rdname fitting_methods
#' @keywords internal
fitting_method_newton <- function(design, outcome, option) {
  max_iter <- if (!is.null(option$max_iter)) option$max_iter else 50  
  epsilon_abs <- if (!is.null(option$epsilon_abs)) option$epsilon_abs else 1e-6  
  epsilon_rel <- if (!is.null(option$epsilon_rel)) option$epsilon_rel else 1e-6  
  epsilon_small <- 1e-8  

  model_type <- if (!is.null(option$model)) option$model else "logistic"
  neg_likelihood_funcs <- get_neg_likelihood_functions(model_type)
  
  if (model_type == "logistic") {
    p <- mean(outcome)
    if (p > 0 && p < 1) {
      init_intercept <- log(p / (1 - p))
      beta <- rep(0, ncol(design))
      beta[1] <- init_intercept
    } else {
      beta <- rep(0, ncol(design))
    }
  } else if (model_type == "linear") {
    beta <- rep(0, ncol(design))
  }
  
  neg_log_likelihood_old <- neg_likelihood_funcs$neg_log_likelihood(beta, design, outcome)
  iter <- 0
  converged <- FALSE

  while (!converged && iter < max_iter) {
    iter <- iter + 1L
    
    neg_grad <- neg_likelihood_funcs$neg_gradient(beta, design, outcome)
    H <- neg_likelihood_funcs$hessian(beta, design, outcome)

    beta_update <- tryCatch(
      solve(H, neg_grad),  
      error = function(e) {
        warning("Hessian is singular, using small step gradient update instead.")
        return(neg_grad * 0.01) 
      }
    )

    beta <- beta - beta_update 
    neg_log_likelihood_new <- neg_likelihood_funcs$neg_log_likelihood(beta, design, outcome)
    change <- abs(neg_log_likelihood_new - neg_log_likelihood_old)  
    rel_change <- change / (abs(neg_log_likelihood_old) + epsilon_small)  

    converged <- (change < epsilon_abs || rel_change < epsilon_rel)
    neg_log_likelihood_old <- neg_log_likelihood_new  
  }

  if (converged) {
    message("Newton's method converged at iteration ", iter)
  } else {
    warning("Newton's method reached the maximum iteration limit (", max_iter, ") without full convergence.")
  }
  
  return(beta)
}
