#' Extract model coefficients from a hiper_glm object
#' @param object A `hiperglm` object
#' @param ... Additional arguments
#' @return A numeric vector of coefficients
#' @export
coef.hiperglm <- function(object, ...) {
  if (!is.list(object) || is.null(object$coefficients)) {
    stop("Invalid hiperglm object: coefficients not found.")
  }
  return(object$coefficients)
}
