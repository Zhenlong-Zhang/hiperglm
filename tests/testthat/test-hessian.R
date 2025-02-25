test_that("Hessian calculation for logistic regression is correct", {
  set.seed(42)  

  n_obs <- 100
  n_pred <- 5
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 42)
  X_test <- data$design
  y_test <- data$outcome
  beta_test <- rep(0, ncol(X_test)) 

  H_analytic <- hessian_logistic(beta_test, X_test, y_test)

  H_numerical <- numerical_hessian(
    function(b) neg_log_likelihood_logistic(b, X_test, y_test), 
    beta_test
  )

  hessian_error <- abs(H_analytic - H_numerical)

  expect_true(sum(hessian_error) < 5e-2)
})

