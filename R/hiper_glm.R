#' Fit a linear model using either pseudo-inverse or BFGS optimization
#' @param X Design matrix (n x p)
#' @param y Response variable (n x 1)
#' @param method Optimization method: "pseudo_inverse" or "BFGS"
#' @return A `hiperglm` object containing estimated coefficients
#' @export
hiper_glm <- function(X, y, method = c("pseudo_inverse", "BFGS")) {
  method <- match.arg(method)
  
  if (method == "pseudo_inverse") {
    beta_hat <- solve(t(X) %*% X, t(X) %*% y)  
  } else if (method == "BFGS") {
    log_likelihood <- function(beta, X, y) {
      residuals <- y - X %*% beta
      return(sum(residuals^2) / 2)  
    }

    grad <- function(beta, X, y) {
      residuals <- y - X %*% beta
      return(-t(X) %*% residuals)  
    }

    beta_init <- solve(t(X) %*% X, t(X) %*% y)  
    optim_res <- optim(beta_init, log_likelihood, grad, X = X, y = y, method = "BFGS", control = list(fnscale = 1))
    beta_hat <- optim_res$par
  }

  return(structure(list(coefficients = beta_hat, method = method), class = "hiperglm"))
}
