# where to get the function needed for testing
source("utils.R")

test_that("are_all_close works correctly", {
  # test 1: value close，should be TRUE
  expect_true(are_all_close(1.000001, 1.000002, rel_tol = 1e-5, abs_tol = 1e-5))

  # test 2: big difference，should be  FALSE
  expect_false(are_all_close(1, 2, rel_tol = 1e-5))

  # test 3: abs is wrong should be FALSE
  expect_false(are_all_close(1, 1.1, abs_tol = 0.05))
})
