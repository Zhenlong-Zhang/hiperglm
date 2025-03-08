test_that("eigenLeastSquaresWeighted returns correct weighted LS results", {
  set.seed(42)
  n <- 100
  p <- 3
  X <- matrix(rnorm(n * p), n, p)
  beta_true <- c(1, 2, -1)
  y <- as.vector(X %*% beta_true + rnorm(n, sd = 0.5))
  w <- runif(n, 0.5, 1.5)
  
  baseline <- get_weighted_ls_baseline(X, y, w)
  res_cpp <- eigenLeastSquaresWeighted(X, y, w)
  
  expect_true(are_all_close(as.numeric(res_cpp$coef), as.numeric(baseline$coef)))
  expect_true(are_all_close(as.numeric(res_cpp$fitted), as.numeric(baseline$fitted)))
  expect_true(are_all_close(as.numeric(res_cpp$resid), as.numeric(baseline$resid)))
})
