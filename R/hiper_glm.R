#' Generalized Linear Model Estimation
#'
#' This function fits a generalized linear model using different optimization methods.
#'
#' @param design A matrix representing the design matrix (independent variables).
#' @param outcome A vector representing the outcome variable (dependent variable).
#' @param model A character string specifying the model type: "logistic" or "linear".
#' If NULL, the function will infer the model type.
#' @param option A list of options including the optimization method ("QR", "BFGS", or "Newton"),
#' maximum iterations, and convergence tolerances.
#'
#' @return A list containing the estimated coefficients, method used, and model type.
#' @export

hiper_glm <- function(design, outcome, model = NULL, option = list()) {
  if (is.null(model)) {
    if (all(outcome %in% c(0, 1))) {
      model <- "logistic"
    } else {
      model <- "linear"
    }
  }

  if (nrow(design) != length(outcome)) {
    stop("Error: The number of rows in 'design' must match the length of 'outcome'.")
  }
  if (!is.matrix(design)) {
    stop("Error: 'design' must be a matrix.")
  }
  if (!is.vector(outcome)) {
    stop("Error: 'outcome' must be a vector.")
  }

  allowed_methods <- c("QR", "BFGS", "Newton")

  if (is.null(option)) {
    option <- list()
  }
  
  option <- modifyList(
    list(method = "QR", model = model, max_iter = 50, epsilon_abs = 1e-6, epsilon_rel = 1e-6),
    option
  )
  
  method <- match.arg(option$method, allowed_methods)

  if (!is.null(option$control) && !is.list(option$control)) {
    stop("Error: 'option$control' must be a list if provided.")
  }

  if (method == "QR") {
    if (model == "logistic") {
      stop("Error: QR method is not supported for logistic regression.")
    }
    beta_hat <- tryCatch(
      fitting_method_qr(design, outcome),
      error = function(e) {
        stop("Error in QR fitting: ", e$message)
      }
    )
  } else if (method == "BFGS") {
    beta_hat <- fitting_method_bfgs(design, outcome, option)
  } else if (method == "Newton") {
    if (model != "logistic") {
      stop("Error: Newton's method is only implemented for logistic regression.")
    }
    beta_hat <- fitting_method_newton(design, outcome, option)
  }

  return(structure(
  list(
    coefficients = beta_hat,
    method = method,
    model = model,
    design = design,
    outcome = outcome
  ),
  class = "hiperglm"
  ))
}