test_that("Confidence intervals for linear model have correct coverage", {
  skip_if_not(Sys.getenv("MANUAL_TEST") != "")

  set.seed(123)
  m <- 500
  nominal <- 0.95
  p <- 3     
  count <- 0
  
  for(i in 1:m) {
    data <- simulate_data(n_obs = 100, n_pred = p, model = "linear")
    fit <- hiper_glm(data$design, data$outcome, model = "linear", option = list(method = "QR"))
    beta_hat <- coef(fit)
    Sigma_hat <- vcov(fit)

    q <- qchisq(nominal, df = length(beta_hat))
    
    diff <- data$coef_true - beta_hat

    Sigma_inv_half <- solve(chol(Sigma_hat))
    stat <- sum((Sigma_inv_half %*% diff)^2)
    
    if(stat < q) count <- count + 1
  }
  
  coverage <- count / m

  sigma <- sqrt(nominal * (1 - nominal) / m)
  

  expect_true(abs(coverage - nominal) < 3 * sigma)
})

test_that("Confidence intervals for logistic model have correct coverage", {
  skip_if_not(Sys.getenv("MANUAL_TEST") != "")
  
  set.seed(123)
  m <- 500 
  nominal <- 0.95
  p <- 3  
  count <- 0
  
  for(i in 1:m) {
    data <- simulate_data(n_obs = 100, n_pred = p, model = "logit")
    fit <- hiper_glm(data$design, data$outcome, model = "logistic", 
                     option = list(method = "Newton", max_iter = 100))
    beta_hat <- coef(fit)
    Sigma_hat <- vcov(fit)
    
    q <- qchisq(nominal, df = length(beta_hat))
    diff <- data$coef_true - beta_hat
    Sigma_inv_half <- solve(chol(Sigma_hat))
    stat <- sum((Sigma_inv_half %*% diff)^2)
    
    if(stat < q) count <- count + 1
  }
  
  coverage <- count / m
  sigma <- sqrt(nominal * (1 - nominal) / m)
  expect_true(abs(coverage - nominal) < 3 * sigma)
})
