test_that("MLE estimated via pseudo-inverse and BFGS are close", {
  data <- simulate_data(n_obs = 100, n_pred = 5, seed = 42)

  #  pseudo-inverse for MLE
  model_pinv <- hiper_glm(data$design, data$outcome, option = list(method = "pseudo_inverse"))
  beta_pinv <- coef(model_pinv)

  # BFGS for MLE
  model_bfgs <- hiper_glm(data$design, data$outcome, option = list(method = "BFGS"))
  beta_bfgs <- coef(model_bfgs)

  # compare MLE 
  expect_true(are_all_close(beta_pinv, beta_bfgs, abs_tol = 1e-4, rel_tol = 1e-4))
})
