test_that("Analytical gradient matches numerical gradient for BFGS", {
  data <- simulate_data(n_obs = 100, n_pred = 5, seed = 42)

  beta_test <- rep(0, ncol(data$design))  # initialize parameters
  grad_analytical <- neg_grad(beta_test, data$design, data$outcome)
  grad_numerical <- approx_grad(function(b) neg_log_likelihood(b, data$design, data$outcome), beta_test)

  expect_true(are_all_close(grad_analytical, grad_numerical, abs_tol = 1e-4, rel_tol = 1e-4))
})
