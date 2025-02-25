#' Sigmoid (expit) function
#' @description Computes the expit (logistic sigmoid) function.
#' @param x Numeric vector or matrix.
#' @return Transformed values using the expit function.
#' @keywords internal
expit <- function(x) {
  return(1 / (1 + exp(-x)))
}