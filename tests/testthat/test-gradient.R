test_that("Analytical gradient matches numerical gradient for BFGS", {
  data <- simulate_data(n_obs = 100, n_pred = 5, seed = 42)

  log_likelihood <- function(beta) {
    residuals <- data$outcome - data$design %*% beta
    return(sum(residuals^2) / 2)  # **负对数似然**
  }

  grad <- function(beta) {
    residuals <- data$outcome - data$design %*% beta
    return(-t(data$design) %*% residuals)  # **梯度**
  }

  beta_test <- rep(0, ncol(data$design))  # **初始化参数**
  grad_analytical <- grad(beta_test)
  grad_numerical <- approx_grad(log_likelihood, beta_test)

  expect_true(are_all_close(grad_analytical, grad_numerical, abs_tol = 1e-4, rel_tol = 1e-4))
})
