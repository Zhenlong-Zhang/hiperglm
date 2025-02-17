# Load necessary functions
source("utils.R")

# Test 1: Check if values are close (should return TRUE)
test_that("are_all_close returns TRUE for values that are close", {
  expect_true(are_all_close(1.000001, 1.000002, rel_tol = 1e-5, abs_tol = 1e-5))
})

# Test 2: Check if large difference returns FALSE (should return FALSE)
test_that("are_all_close returns FALSE for values with a big difference", {
  expect_false(are_all_close(1, 2, rel_tol = 1e-5))
})

# Test 3: Check if absolute error is above tolerance (should return FALSE)
test_that("are_all_close returns FALSE when absolute error exceeds tolerance", {
  expect_false(are_all_close(1, 1.1, abs_tol = 0.05))
})
