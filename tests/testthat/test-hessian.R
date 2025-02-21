numerical_hessian <- function(f, beta, epsilon = 1e-5) {  
  n <- length(beta)
  H_approx <- matrix(0, n, n)

  for (i in 1:n) {
    for (j in 1:n) {
      perturbation_scale <- max(abs(beta), 1) * epsilon  
      
      beta_ij1 <- beta
      beta_ij2 <- beta
      beta_ij3 <- beta
      beta_ij4 <- beta

      beta_ij1[i] <- beta_ij1[i] + perturbation_scale
      beta_ij1[j] <- beta_ij1[j] + perturbation_scale

      beta_ij2[i] <- beta_ij2[i] + perturbation_scale
      beta_ij2[j] <- beta_ij2[j] - perturbation_scale

      beta_ij3[i] <- beta_ij3[i] - perturbation_scale
      beta_ij3[j] <- beta_ij3[j] + perturbation_scale

      beta_ij4[i] <- beta_ij4[i] - perturbation_scale
      beta_ij4[j] <- beta_ij4[j] - perturbation_scale

      H_approx[i, j] <- (f(beta_ij1) - f(beta_ij2) - f(beta_ij3) + f(beta_ij4)) / (4 * perturbation_scale^2)
    }
  }

  return(H_approx)
}

# test

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

