test_that("Analytical gradient matches numerical gradient for BFGS (Linear Regression)", {
  data <- simulate_data(n_obs = 100, n_pred = 5, seed = 42)

  beta_test <- rep(0, ncol(data$design))  # initialize parameters
  grad_analytical <- neg_gradient_linear(beta_test, data$design, data$outcome)
  grad_numerical <- approx_grad(function(b) neg_log_likelihood_linear(b, data$design, data$outcome), beta_test)

  expect_true(are_all_close(grad_analytical, grad_numerical, abs_tol = 1e-4, rel_tol = 1e-4))
})

test_that("Analytical gradient matches numerical gradient for BFGS (Logistic Regression)", {
  set.seed(42)
  X_logistic <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  y_logistic <- rbinom(100, 1, 0.5)

  beta_test <- rep(0, ncol(X_logistic))  # initialize parameters
  grad_analytical <- neg_gradient_logistic(beta_test, X_logistic, y_logistic)
  grad_numerical <- approx_grad(function(b) neg_log_likelihood_logistic(b, X_logistic, y_logistic), beta_test)

  expect_true(are_all_close(grad_analytical, grad_numerical, abs_tol = 1e-4, rel_tol = 1e-4))
})
