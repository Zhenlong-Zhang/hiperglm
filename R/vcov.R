#' Variance-Covariance Matrix for hiper_glm Objects
#'
#' @param object An object of class "hiperglm".
#' @param ... Further arguments passed to or from other methods (not used here).
#'
#' @return A matrix representing the variance-covariance of the estimated coefficients.
#' @export
vcov.hiperglm <- function(object, ...) {
  # should only work for linear. i think this warning should be here because
  if (object$model != "linear") {
    stop("vcov() is only implemented for linear models in this step.")
  }

  X <- object$design    
  y <- object$outcome   
  beta_hat <- object$coefficients
  n <- nrow(X)
  p <- ncol(X)

  residuals <- y - X %*% beta_hat
  sigma2 <- (n / (n - p)) * sum(residuals^2)

  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R <- R[, order(qr_decomp$pivot)]

  xtx_inv <- chol2inv(R)
  vcov_mat <- sigma2 * xtx_inv
  return(vcov_mat)
}
