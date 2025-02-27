test_that("vcov.hiperglm produces correct output for linear models", {
  set.seed(123)
  X_linear <- cbind(1, rnorm(20))
  beta_linear <- c(2, 3)
  y_linear <- as.vector(X_linear %*% beta_linear + rnorm(20))
  
  fit_linear <- list(design = X_linear, outcome = y_linear, coefficients = beta_linear, model = "linear")
  class(fit_linear) <- "hiperglm"
  
  vcov_linear <- vcov.hiperglm(fit_linear)
  
  expect_equal(dim(vcov_linear), c(ncol(X_linear), ncol(X_linear)))
  expect_true(are_all_close(vcov_linear, t(vcov_linear), abs_tol = 1e-6, rel_tol = 1e-6))
})

test_that("vcov.hiperglm produces correct output for logistic models", {
  set.seed(123)
  X_logistic <- cbind(1, rnorm(20))
  beta_logistic <- c(-1, 2)
  p_logistic <- expit(as.vector(X_logistic %*% beta_logistic))
  y_logistic <- rbinom(20, 1, p_logistic)
  
  fit_logistic <- list(design = X_logistic, outcome = y_logistic, coefficients = beta_logistic, model = "logistic")
  class(fit_logistic) <- "hiperglm"
  
  vcov_logistic <- vcov.hiperglm(fit_logistic)
  
  expect_equal(dim(vcov_logistic), c(ncol(X_logistic), ncol(X_logistic)))
  expect_true(are_all_close(vcov_logistic, t(vcov_logistic), abs_tol = 1e-6, rel_tol = 1e-6))
})
