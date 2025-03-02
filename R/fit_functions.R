#' Internal Fit Function for hiperglm
#'
#' @section QR Decomposition Fitting Method:
#'
#' Implements the least squares solution using QR decomposition for linear.
#'
#' @rdname fitting_methods
#' @keywords internal
fitting_method_qr <- function(design, outcome) {
  qr_decomp <- qr(design)
  beta <- qr.solve(qr_decomp, outcome)
  return(beta)
}

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
  
  if (model_type == "logistic") {
    p <- mean(outcome)
    if (p > 0 && p < 1) {
      init_intercept <- log(p / (1 - p))
      beta_init <- rep(0, ncol(design))
      beta_init[1] <- init_intercept
    } else {
      beta_init <- rep(0, ncol(design))
    }
    neg_log_likelihood <- neg_log_likelihood_logistic
    neg_gradient <- neg_gradient_logistic
  } else if (model_type == "linear") {
    beta_init <- rep(0, ncol(design))
    neg_log_likelihood <- neg_log_likelihood_linear
    neg_gradient <- neg_gradient_linear
  }
  
  optim_res <- optim(
    beta_init, 
    fn = function(b) neg_log_likelihood(b, design, outcome),
    gr = function(b) neg_gradient(b, design, outcome),
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
  
  if (model_type == "logistic") {
    p <- mean(outcome)
    if (p > 0 && p < 1) {
      init_intercept <- log(p / (1 - p))
      beta <- rep(0, ncol(design))
      beta[1] <- init_intercept
    } else {
      beta <- rep(0, ncol(design))
    }
    neg_log_likelihood <- neg_log_likelihood_logistic
    neg_gradient <- neg_gradient_logistic
    hessian <- hessian_logistic
  } else if (model_type == "linear") {
    beta <- rep(0, ncol(design))
    neg_log_likelihood <- neg_log_likelihood_linear
    neg_gradient <- neg_gradient_linear
    hessian <- NULL
  }
  
  neg_log_likelihood_old <- neg_log_likelihood(beta, design, outcome)
  converged <- FALSE
  
  for (iter in 1:max_iter) {
    step_result <- take_one_newton_step(
      beta, 
      design, 
      outcome, 
      neg_gradient, 
      hessian, 
      neg_log_likelihood, 
      solver = "qr",
      epsilon_small = epsilon_small
    )
    new_beta <- step_result$beta
    new_neg_log_likelihood <- step_result$neg_log_likelihood

    change <- abs(new_neg_log_likelihood - neg_log_likelihood_old)
    rel_change <- change / (abs(neg_log_likelihood_old) + epsilon_small)
    
    if (change < epsilon_abs || rel_change < epsilon_rel) {
      beta <- new_beta
      converged <- TRUE
      
      break
    }
    
    beta <- new_beta
    neg_log_likelihood_old <- new_neg_log_likelihood  
  }
  
  if (!converged) {
    warning("Newton's method reached the maximum iteration limit (", max_iter, ") without full convergence.")
  }
  
  return(beta)
}
