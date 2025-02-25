#' --------linear model
#'-------------------------------------------------------------
# Compute negative log-likelihood for linear regression
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Response variable
# @return Negative log-likelihood value
neg_log_likelihood_linear <- function(beta, X, y) {
  residuals <- y - X %*% beta
  return(sum(residuals^2) / 2)
}
# Compute gradient of the negative log-likelihood
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Response variable
# @return Gradient vector
neg_gradient_linear <- function(beta, X, y) {
  residuals <- y - X %*% beta
  return(-t(X) %*% residuals)
}
#' ------------------------logic model
#'-------------------------------------------------------------
# Compute negative log-likelihood for logistic regression
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Binary response variable (0 or 1)
# @return Negative log-likelihood value
neg_log_likelihood_logistic <- function(beta, X, y) {
  p <- 1 / (1 + exp(-X %*% beta))  
  log_lik <- sum(y * log(p) + (1 - y) * log(1 - p))  
  return(-log_lik)  
}
# Compute gradient of the negative log-likelihood for logistic regression
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Binary response variable (0 or 1)
# @return Gradient vector
neg_gradient_logistic <- function(beta, X, y) {
  p <- 1 / (1 + exp(-X %*% beta))  
  grad <- t(X) %*% (y - p)
  return(-grad)  
}
#' --------------------------- hessian
#'-------------------------------------------------------------
# Compute Hessian matrix for logistic regression
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Binary response variable (0 or 1)
# @return Hessian matrix
hessian_logistic <- function(beta, X, y) {
  p <- 1 / (1 + exp(-X %*% beta))

  p <- pmax(pmin(p, 1 - 1e-8), 1e-8)
  W_diag <- as.vector(p * (1 - p)) 

  H <- crossprod(X, W_diag * X)  
  return(H)
}



#' --------------choose based on model
#'-------------------------------------------------------------
get_neg_likelihood_functions <- function(model) {
  neg_likelihoods <- list(
    linear = list(
      neg_log_likelihood = neg_log_likelihood_linear,
      neg_gradient = neg_gradient_linear,
      hessian = NULL  # linear deos not need Hessian
    ),
    logistic = list(
      neg_log_likelihood = neg_log_likelihood_logistic,
      neg_gradient = neg_gradient_logistic,
      hessian = hessian_logistic  # Hessian applu Newton’s 
    )
  )
  
  if (!model %in% names(neg_likelihoods)) {
    stop("Unsupported model type: ", model)
  }
  
  return(neg_likelihoods[[model]])
}
