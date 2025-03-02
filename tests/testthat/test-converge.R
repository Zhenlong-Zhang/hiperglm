test_that("Confidence intervals for linear model have correct coverage", {
  skip_if_not(Sys.getenv("MANUAL_TEST") != "")
  set.seed(123)
  test_confint_coverage("linear", option = list(method = "QR"))
})

test_that("Confidence intervals for logistic model have correct coverage", {
  skip_if_not(Sys.getenv("MANUAL_TEST") != "")
  set.seed(123)
  test_confint_coverage("logistic", option = list(method = "Newton", max_iter = 100))
})
