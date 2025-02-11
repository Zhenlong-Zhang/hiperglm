test_that("are_all_close works correctly", {
  are_all_close <- function(a, b, rel_tol = 1e-5, abs_tol = 1e-8) {
    diff <- abs(a - b)
    return(all(diff < abs_tol | diff / (abs(b) + abs_tol) < rel_tol))
  }

  # test 1: value close，should be TRUE
  expect_true(are_all_close(1.000001, 1.000002, rel_tol = 1e-5))

  # test 2: big difference，should be  FALSE
  expect_false(are_all_close(1, 2, rel_tol = 1e-5))

  # test 3: abs is wrong should be FALSE
  expect_false(are_all_close(1, 1.1, abs_tol = 0.05))
})
