#' Check if two numeric vectors are close to each other
#' @param v First numeric vector
#' @param w Second numeric vector
#' @param abs_tol Absolute tolerance
#' @param rel_tol Relative tolerance
#' @return Boolean indicating whether all elements are within tolerance
#' @keywords internal
are_all_close <- function(v, w, abs_tol = 1e-6, rel_tol = 1e-6) {
  abs_diff <- abs(v - w)
  are_all_within_atol <- all(abs_diff < abs_tol)
  are_all_within_rtol <- all(abs_diff < rel_tol * pmax(abs(v), abs(w)))
  return(are_all_within_atol && are_all_within_rtol)
}
#' Generate simulated data for testing
#' @param n_obs Number of observations
#' @param n_pred Number of predictors
#' @param model Model type (default: "linear")
#' @param intercept Optional intercept value
#' @param coef_true True coefficients (optional)
#' @param design Design matrix (optional)
#' @param seed Random seed for reproducibility
#' @param signal_to_noise Signal-to-noise ratio
#' @return A list with `design`, `outcome`, and `coef_true`
#' @keywords internal
simulate_data <- function(
    n_obs, n_pred, model = "linear", intercept = NULL, 
    coef_true = NULL, design = NULL, seed = NULL, option = list()
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if ((model != "linear")  && !is.null(option$signal_to_noise)) {
    warning(paste(
      "The `signal_to_noise` option is currently unsupported for",
      "non-linear models and will be ignored."
    ))
  }
  if (is.null(coef_true)) {
    coef_true <- rnorm(n_pred, sd = 1 / sqrt(n_pred))
  } 
  if (is.null(design)) {
    design <- matrix(rnorm(n_obs * n_pred), nrow = n_obs, ncol = n_pred)
  }
  if (!is.null(intercept)) {
    if (!is.numeric(intercept)) {
      stop("The intercept argument must be numeric.")
    }
    coef_true <- c(intercept, coef_true)
    design <- cbind(rep(1, n_obs), design)
  }
  expected_mean <- as.vector(design %*% coef_true)
  if (model == 'linear') {
    signal_to_noise <- option$signal_to_noise
    if (is.null(signal_to_noise)) {
      signal_to_noise <- 0.1
    }
    noise_magnitude <- sqrt(var(expected_mean) / signal_to_noise^2)
    noise <- noise_magnitude * rnorm(n_obs)
    outcome <- expected_mean + noise 
  } else {
    n_trial <- option$n_trial
    prob <- 1 / (1 + exp(-expected_mean))
    if (is.null(n_trial)) {
      outcome <- rbinom(n_obs, 1, prob)
    } else {
      n_success <- rbinom(n_obs, n_trial, prob)
      outcome <- list(n_success = n_success, n_trial = n_trial)
    }
  }
  return(list(design = design, outcome = outcome, coef_true = coef_true))
}

#' Compute numerical gradient using finite difference method (centered difference)
#' @param func Function to compute gradient for
#' @param x Point at which gradient is evaluated
#' @param dx Small perturbation for finite difference
#' @return Numerical gradient (vector)
#' @keywords internal
approx_grad <- function(func, x, dx = .Machine$double.eps^(1/3)) {
  numerical_grad <- rep(0, length(x))
  for (i in seq_along(x)) {
    e_i <- rep(0, length(x))
    e_i[i] <- 1
    
    f_plus  <- func(x + dx * e_i)
    f_minus <- func(x - dx * e_i)
    
    numerical_grad[i] <- (f_plus - f_minus) / (2 * dx)
  }
  return(numerical_grad)
}


#' Compute numerical Hessian using finite difference method (centered difference)
#' @param f Function to compute Hessian for
#' @param beta Coefficient vector at which Hessian is evaluated
#' @param epsilon Small perturbation for finite difference (default: 1e-5)
#' @return Approximate Hessian matrix
#' @keywords internal
numerical_hessian <- function(f, beta, epsilon = 1e-5) {  
  n <- length(beta)
  H_approx <- matrix(0, n, n)
  
  perturbation_scale <- epsilon  
  
  for (i in 1:n) {
    for (j in 1:n) {
      beta_ij1 <- beta
      beta_ij2 <- beta
      beta_ij3 <- beta
      beta_ij4 <- beta

      beta_ij1[i] <- beta_ij1[i] + perturbation_scale
      beta_ij1[j] <- beta_ij1[j] + perturbation_scale

      beta_ij2[i] <- beta_ij2[i] + perturbation_scale
      beta_ij2[j] <- beta_ij2[j] - perturbation_scale

      beta_ij3[i] <- beta_ij3[i] - perturbation_scale
      beta_ij3[j] <- beta_ij3[j] + perturbation_scale

      beta_ij4[i] <- beta_ij4[i] - perturbation_scale
      beta_ij4[j] <- beta_ij4[j] - perturbation_scale

      H_approx[i, j] <- (f(beta_ij1) - f(beta_ij2) - f(beta_ij3) + f(beta_ij4)) / (4 * perturbation_scale^2)
    }
  }

  return(H_approx)
}

#' Test convergence via confidence interval coverage
#'
#' This internal function simulates data to check that the computed confidence intervals
#' achieve the nominal coverage level. It indirectly tests the convergence of the estimation
#' algorithm (e.g., QR for linear models or Newton's method for logistic models).
#'
#' @param model_type "linear" or "logistic".
#' @param option Options to be passed to hiper_glm.
#' @param m Number of simulations (default: 500).
#' @param n_obs Number of observations per simulation (default: 100).
#' @param n_pred Number of predictors (default: 3).
#' @param nominal Nominal coverage probability (default: 0.95).
#' @param threshold Tolerance threshold (default: 3).
#'
#' @return No return value. Uses expect_true to assert that the observed coverage is within tolerance.
#'
#' @keywords internal
test_confint_coverage <- function(model_type, option, m = 500, n_obs = 100, n_pred = 3, 
                                  nominal = 0.95, threshold = 3) {
  count <- 0
  
  sim_model <- if(model_type == "linear") "linear" else "logit"
  fit_model <- if(model_type == "linear") "linear" else "logistic"
  
  for(i in 1:m) {
    data <- simulate_data(n_obs = n_obs, n_pred = n_pred, model = sim_model)
    fit <- hiper_glm(data$design, data$outcome, model = fit_model, option = option)
    beta_hat <- coef(fit)
    Sigma_hat <- vcov(fit)
    
    diff <- data$coef_true - beta_hat
    stat <- sum(backsolve(chol(Sigma_hat), diff)^2)
    
    q <- qchisq(nominal, df = length(beta_hat))
    if(stat < q) count <- count + 1
  }
  
  coverage <- count / m
  sigma <- sqrt(nominal * (1 - nominal) / m)
  expect_true(abs(coverage - nominal) < threshold * sigma)
}