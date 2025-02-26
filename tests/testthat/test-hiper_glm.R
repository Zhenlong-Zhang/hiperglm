test_that("MLE estimated via pseudo-inverse and BFGS are close(Linear Regression)", {
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

test_that("BFGS correctly estimates MLE for Logistic Regression", {
  set.seed(42)
  X_logistic <- cbind(1, matrix(rnorm(100 * 5), nrow = 100, ncol = 5))
  y_logistic <- rbinom(100, 1, 0.5)

  model_bfgs <- hiper_glm(X_logistic, y_logistic, model = "logistic", option = list(method = "BFGS"))
  beta_bfgs <- coef(model_bfgs)

  # makesure beta_bfgs reasonable
  expect_type(beta_bfgs, "double")
  expect_equal(length(beta_bfgs), ncol(X_logistic))

  # compare
  glm_model <- glm(y_logistic ~ X_logistic - 1, family = binomial)
  glm_coef <- coef(glm_model)

  expect_true(are_all_close(beta_bfgs, glm_coef, abs_tol = 1e-4, rel_tol = 1e-4))
})
