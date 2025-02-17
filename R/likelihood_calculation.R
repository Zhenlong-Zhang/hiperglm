# Compute negative log-likelihood for linear regression
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Response variable
# @return Negative log-likelihood value
neg_log_likelihood <- function(beta, X, y) {
  residuals <- y - X %*% beta
  return(sum(residuals^2) / 2)
}

# Compute gradient of the negative log-likelihood
# @param beta Coefficients (vector)
# @param X Design matrix
# @param y Response variable
# @return Gradient vector
neg_grad <- function(beta, X, y) {
  residuals <- y - X %*% beta
  return(-t(X) %*% residuals)
}
