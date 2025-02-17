# Load the fitting methods
source("R/fit_functions.R")


hiper_glm <- function(design, outcome, model = "linear", option = list()) {
  
  # Check input validity
  if (nrow(design) != length(outcome)) {
    stop("Error: The number of rows in 'design' must match the length of 'outcome'.")
  }
  if (!is.matrix(design)) {
    stop("Error: 'design' must be a matrix.")
  }
  if (!is.vector(outcome)) {
    stop("Error: 'outcome' must be a vector.")
  }
  
  # Allowed methods
  allowed_methods <- c("pseudo_inverse", "BFGS")
  option <- modifyList(list(method = "pseudo_inverse"), option)
  method <- match.arg(option$method, allowed_methods)
  
  # Control checks
  if (!is.null(option$control) && !is.list(option$control)) {
    stop("Error: 'option$control' must be a list if provided or if none a default method of pseudo inverse will be performed.")
  }

  # Fit the model
  if (method == "pseudo_inverse") {
    beta_hat <- fitting_method_pseudo_inverse(design, outcome)  
  } else if (method == "BFGS") {
    beta_hat <- fitting_method_bfgs(design, outcome, option)
  }

  return(structure(list(coefficients = beta_hat, method = method, model = model), class = "hiperglm"))
}
