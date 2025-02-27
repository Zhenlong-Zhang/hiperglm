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
#' This function performs a single Newton update step using either a QR-based weighted
#' least squares solver or a pseudo-inverse solver. This facilitates unit-testing one step
#' of Newton's method.
#'
#' @param beta Numeric vector of current coefficients.
#' @param design Numeric design matrix.
#' @param outcome Numeric vector of outcomes.
#' @param neg_gradient Function to compute the negative gradient.
#' @param hessian Function to compute the Hessian matrix.
#' @param neg_log_likelihood Function to compute the negative log-likelihood.
#' @param solver Character string: "qr" (default) or "pseudo".
#' @param epsilon_small Small numeric value to avoid division by zero (default: 1e-8).
#'
#' @return A list with elements:
#' \describe{
#'   \item{beta}{The updated coefficient vector.}
#'   \item{neg_log_likelihood}{The negative log-likelihood at the updated beta.}
#' }
#'
#' @keywords internal
take_one_newton_step <- function(beta, design, outcome, neg_gradient, hessian, neg_log_likelihood,
                                  solver = "qr", epsilon_small = 1e-8) {
  solver <- match.arg(solver, choices = c("qr", "pseudo"))
  
  p <- expit(as.vector(design %*% beta))
  
  if (solver == "pseudo") {
    neg_grad <- neg_gradient(beta, design, outcome) 
    delta <- tryCatch(
      solve(hessian(beta, design, outcome), neg_grad),
      error = function(e) {
        warning("Hessian is singular, using small step gradient update instead.")
        neg_grad * 0.01
      }
    )
    new_beta <- beta - delta
  } else if (solver == "qr") {
    W_vec <- p * (1 - p)
    W_vec[W_vec < 1e-8] <- 1e-8
    sqrt_W <- sqrt(W_vec)
    tilde_X <- design * sqrt_W
    tilde_y <- (outcome - p) / sqrt_W
    qr_decomp <- qr(tilde_X)
    R <- qr.R(qr_decomp)
    R <- R[, order(qr_decomp$pivot)]
    delta <- qr.coef(qr_decomp, tilde_y)
    if (any(is.na(delta))) {
      warning("QR decomposition failed, falling back to pseudo-inverse solver.")
      neg_grad <- neg_gradient(beta, design, outcome)
      delta <- tryCatch(
        solve(hessian(beta, design, outcome), neg_grad),
        error = function(e) {
          warning("Hessian is singular, using small step gradient update instead.")
          neg_grad * 0.01
        }
      )
    }
    new_beta <- beta + delta
  }
  
  new_neg_log_likelihood <- neg_log_likelihood(new_beta, design, outcome)
  
  return(list(beta = new_beta, neg_log_likelihood = new_neg_log_likelihood))
}

#' Extract the R Factor from a QR Decomposition
#'
#' This function performs a QR decomposition on a given matrix \eqn{X}, extracts the upper-triangular
#' factor \eqn{R}, and reorders its columns according to the pivoting information. This is useful for
#' ensuring that \eqn{X \approx Q \times R} holds, even when column pivoting is used.
#'
#' @param X A numeric matrix.
#'
#' @return The reordered upper-triangular matrix \eqn{R}.
#' @keywords internal
extract_R <- function(X) {
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R <- R[, order(qr_decomp$pivot)]
  return(R)
}