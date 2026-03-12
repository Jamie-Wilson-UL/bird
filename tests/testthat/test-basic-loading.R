test_that("package loads", {
  expect_true(requireNamespace("bird", quietly = TRUE))
  
  # Check that basic functions are available
  expect_true(exists("impute"))
  expect_true(exists("bayesian_impute"))
  expect_true(exists("bayes_np_impute"))
  expect_true(exists("complete"))
  expect_true(exists("export"))
})
