test_that("Newton's method converges properly for logistic regression", {
  set.seed(42)
  n_obs <- 100
  n_pred <- 5
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 42)

  # run Newton’s 
  model_newton <- hiper_glm(data$design, data$outcome, option = list(method = "Newton", epsilon_abs = 1e-6, epsilon_rel = 1e-6, max_iter = 100))
  beta_newton <- coef(model_newton)

  expect_type(beta_newton, "double")
  expect_equal(length(beta_newton), ncol(data$design))

  # compare
  glm_model <- glm(data$outcome ~ data$design - 1, family = binomial)
  glm_coef <- coef(glm_model)

  expect_true(are_all_close(beta_newton, glm_coef, abs_tol = 1e-3, rel_tol = 1e-3))
})

test_that("Newton's method warns when non-convergent", {
  set.seed(42)
  n_obs <- 100
  n_pred <- 5
  data <- simulate_data(n_obs, n_pred, model = "logit", seed = 42)

  expect_warning(
    hiper_glm(data$design, data$outcome, option = list(method = "Newton", max_iter = 2)),
    regexp = "Newton's method reached the maximum iteration limit"
  )
})

test_that("QR and pseudo solvers yield similar Newton step", {

  data <- simulate_data(n_obs = 100, n_pred = 1, model = "logistic", seed = 123)
  X <- data$design
  outcome <- data$outcome
  beta_init <- rep(0, ncol(X))
  
  step_qr <- take_one_newton_step(
    beta = beta_init,
    design = X,
    outcome = outcome,
    neg_gradient = neg_gradient_logistic,
    hessian = hessian_logistic,
    neg_log_likelihood = neg_log_likelihood_logistic,
    solver = "qr",
    epsilon_small = 1e-8
  )
  
  step_pseudo <- take_one_newton_step(
    beta = beta_init,
    design = X,
    outcome = outcome,
    neg_gradient = neg_gradient_logistic,
    hessian = hessian_logistic,
    neg_log_likelihood = neg_log_likelihood_logistic,
    solver = "pseudo",
    epsilon_small = 1e-8
  )
  
  expect_true(are_all_close(step_qr$beta, step_pseudo$beta, abs_tol = 1e-6, rel_tol = 1e-6))
})
