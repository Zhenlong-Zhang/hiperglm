# Load the fitting methods
source("R/fit_functions.R")

hiper_glm <- function(design, outcome, model = NULL, option = list()) {
  #  Logistic if outcome 0/1
  if (is.null(model)) {
    if (all(outcome %in% c(0, 1))) {
      model <- "logistic"
    } else {
      model <- "linear"
    }
  }
  # inputs check
  if (nrow(design) != length(outcome)) {
    stop("Error: The number of rows in 'design' must match the length of 'outcome'.")
  }
  if (!is.matrix(design)) {
    stop("Error: 'design' must be a matrix.")
  }
  if (!is.vector(outcome)) {
    stop("Error: 'outcome' must be a vector.")
  }
  allowed_methods <- c("pseudo_inverse", "BFGS", "Newton")
  # null is defaly to be pse
  if (is.null(option)) {
    option <- list()
  }
  option <- modifyList(list(method = "pseudo_inverse", model = model), option)
  method <- match.arg(option$method, allowed_methods)
  if (!is.null(option$control) && !is.list(option$control)) {
    stop("Error: 'option$control' must be a list if provided.")
  }
  if (method == "pseudo_inverse") {
    if (model == "logistic") {
      stop("Error: Pseudo-inverse method is not supported for logistic regression.")
    }
    beta_hat <- tryCatch(
      fitting_method_pseudo_inverse(design, outcome),
      error = function(e) {
        stop("Error in pseudo-inverse fitting: ", e$message)
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

  return(structure(list(coefficients = beta_hat, method = method, model = model), class = "hiperglm"))
}
