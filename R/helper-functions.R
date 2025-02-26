#' Sigmoid (expit) function
#' @description Computes the expit (logistic sigmoid) function.
#' @param x Numeric vector or matrix.
#' @return Transformed values using the expit function.
#' @keywords internal
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}

#' Perform One Newton Update Step (Internal)
#'
#' This function performs a single Newton update step given the current parameter
#' estimate, design matrix, and outcome. It computes the negative gradient and Hessian,
#' then updates the parameter vector accordingly. If the Hessian is singular,
#' a small step gradient update is used.
#'
#' @param beta A numeric vector representing the current coefficients.
#' @param design A numeric matrix representing the design matrix.
#' @param outcome A numeric vector representing the outcome.
#' @param neg_gradient A function to compute the negative gradient given beta, design, and outcome.
#' @param hessian A function to compute the Hessian matrix given beta, design, and outcome.
#' @param neg_log_likelihood A function to compute the negative log-likelihood given beta, design, and outcome.
#' @param epsilon_small A small numeric value to avoid division by zero (default: 1e-8).
#'
#' @return A list containing:
#' \describe{
#'   \item{beta}{The updated coefficient vector.}
#'   \item{neg_log_likelihood}{The negative log-likelihood evaluated at the updated beta.}
#' }
#'
#' @keywords internal
take_one_newton_step <- function(beta, design, outcome, neg_gradient, hessian, neg_log_likelihood, epsilon_small = 1e-8) {
  neg_grad <- neg_gradient(beta, design, outcome)
  H <- hessian(beta, design, outcome)
  beta_update <- tryCatch(
    solve(H, neg_grad),
    error = function(e) {
      warning("Hessian is singular, using small step gradient update instead.")
      return(neg_grad * 0.01)
    }
  )  
  new_beta <- beta - beta_update
  new_neg_log_likelihood <- neg_log_likelihood(new_beta, design, outcome) 
  return(list(beta = new_beta, neg_log_likelihood = new_neg_log_likelihood))
}
