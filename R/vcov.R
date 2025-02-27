#' Variance-Covariance Matrix for hiper_glm Objects
#'
#' Computes the variance-covariance matrix of the estimated coefficients.
#' For linear models, the matrix is given by \(\hat{\sigma}^2 (X^\top X)^{-1}\)
#' with \(\hat{\sigma}^2 = \frac{n}{n-p} \| y - X\hat{\beta} \|^2\).
#' For logistic models, the observed Fisher information is \(X^\top W X\),
#' where \(W\) is a diagonal matrix with entries \(p_i(1-p_i)\).
#' The inverse is computed from the "R" factor of the QR decomposition of
#' \(X\) (linear) or \(\sqrt{W}X\) (logistic).
#'
#' @param object An object of class "hiperglm".
#' @param ... Further arguments passed to or from other methods (not used here).
#'
#' @return A matrix representing the variance-covariance of the estimated coefficients.
#' @export
vcov.hiperglm <- function(object, ...) {
  if (object$model == "linear") {
    X <- object$design
    y <- object$outcome
    beta_hat <- object$coefficients
    n <- nrow(X)
    p <- ncol(X)
    residuals <- y - X %*% beta_hat
    sigma2 <- (n / (n - p)) * sum(residuals^2)
    R <- extract_R(X)
    xtx_inv <- chol2inv(R)
    vcov_mat <- sigma2 * xtx_inv
    return(vcov_mat)
  } else if (object$model == "logistic") {
    X <- object$design
    beta_hat <- object$coefficients
    p <- expit(as.vector(X %*% beta_hat))
    sqrt_W <- sqrt(p * (1 - p))
    tilde_X <- X * sqrt_W
    R <- extract_R(tilde_X)
    fisher_inv <- chol2inv(R)
    return(fisher_inv)
  } else {
    stop("vcov() is only implemented for linear and logistic models.")
  }
}
