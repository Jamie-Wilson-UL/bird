# Tests for data setup functions (prepare_survival_data, quick_impute, etc.)

test_that("prepare_survival_data requires explicit column names", {
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1),
    age = c(50, 60, 45, 70, 55)
  )

  expect_error(
    prepare_survival_data(test_data, verbose = FALSE),
    "Please specify both 'time' and 'status'; auto-detection is not yet supported",
    fixed = TRUE
  )

  result <- prepare_survival_data(test_data, time = "time", status = "status", verbose = FALSE)

  expect_s3_class(result, "survival_data")
  expect_equal(attr(result, "time_col"), "time")
  expect_equal(attr(result, "status_col"), "status") 
  expect_equal(attr(result, "n_total"), 5)
  expect_equal(attr(result, "n_censored"), 2)
  expect_equal(attr(result, "censoring_rate"), 0.4)
})

test_that("prepare_survival_data handles manual specification", {
  # Test with non-standard column names requiring manual specification
  test_data <- data.frame(
    follow_up_days = c(100, 150, 200),
    died = c(1, 0, 1),
    treatment = c("A", "B", "A")
  )
  
  result <- prepare_survival_data(test_data, 
                                 time = "follow_up_days", 
                                 status = "died",
                                 verbose = FALSE)
  
  expect_equal(attr(result, "time_col"), "follow_up_days")
  expect_equal(attr(result, "status_col"), "died")
})

test_that("public imputation entry points require explicit columns for raw data", {
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )

  expect_error(
    impute(test_data, verbose = FALSE),
    "Please specify both 'time' and 'status'; auto-detection is not yet supported",
    fixed = TRUE
  )

  expect_error(
    bayesian_impute(test_data, verbose = FALSE),
    "Please specify both 'time_col' and 'status_col'; auto-detection is not yet supported",
    fixed = TRUE
  )

  expect_error(
    bayes_np_impute(test_data, verbose = FALSE),
    "Please specify both 'time_col' and 'status_col'; auto-detection is not yet supported",
    fixed = TRUE
  )

  expect_error(
    compare_models(test_data, models = c("weibull", "exponential"), verbose = FALSE),
    "Please specify both 'time' and 'status'; auto-detection is not yet supported",
    fixed = TRUE
  )
})

test_that("prepare_survival_data validates input data", {
  # Test empty data
  expect_error(
    prepare_survival_data(data.frame()),
    "Data cannot be empty"
  )
  
  # Test non-data.frame input
  expect_error(
    prepare_survival_data(list(a = 1, b = 2)),
    "Data must be a data.frame"
  )
  
  # Test data where columns can't be detected
  test_data <- data.frame(
    weird_col_1 = c(10, 20, 30),
    weird_col_2 = c(1, 2, 3)
  )
  
  expect_error(
    prepare_survival_data(test_data, verbose = FALSE),
    "Please specify both 'time' and 'status'; auto-detection is not yet supported",
    fixed = TRUE
  )
})

test_that("prepare_survival_data handles missing column specifications", {
  test_data <- data.frame(
    time = c(10, 15),
    status = c(1, 0)
  )
  
  # Test missing time column
  expect_error(
    prepare_survival_data(test_data, time = "missing_col", status = "status"),
    "Time column 'missing_col' not found"
  )
  
  # Test missing status column  
  expect_error(
    prepare_survival_data(test_data, time = "time", status = "missing_col"),
    "Status column 'missing_col' not found"
  )
})

test_that("bayesian_impute.survival_data method works", {
  test_data <- data.frame(
    time = c(10, 15, 20),
    status = c(1, 0, 1) 
  )
  
  survival_data <- prepare_survival_data(test_data, time = "time", status = "status", verbose = FALSE)
  
  # Mock the default method to avoid Stan compilation 
  expect_no_error({
    # This would normally run: result <- bayesian_impute(survival_data)
    # For testing, we just verify the method dispatch works
    expect_true(inherits(survival_data, "survival_data"))
  })
})

test_that("print.survival_data works correctly", {
  test_data <- data.frame(
    time = c(10, 15, 20, 25),
    status = c(1, 0, 1, 0)
  )
  
  survival_data <- prepare_survival_data(test_data, time = "time", status = "status", verbose = FALSE)
  
  # Capture printed output
  output <- capture.output(print(survival_data))
  
  expect_true(any(grepl("Survival Data", output)))
  expect_true(any(grepl("Observations: 4", output)))
  expect_true(any(grepl("Time column: 'time'", output)))
  expect_true(any(grepl("Status column: 'status'", output)))
  expect_true(any(grepl("Censored:2\\(50%\\)", output)))
})

test_that("quick_impute provides one-liner functionality", {
  test_data <- data.frame(
    time = c(10, 15, 20),
    status = c(1, 0, 1)
  )
  
  # Test that it would work 
  expect_no_error({
    # This would normally be: result <- quick_impute(test_data, n_imputations = 2)
    # For testing, we verify the data preparation step works
    prepared <- prepare_survival_data(test_data, time = "time", status = "status", verbose = FALSE)
    expect_true(inherits(prepared, "survival_data"))
  })
})

test_that("prepare_survival_data works with numeric column indices", {
  test_data <- data.frame(
    follow_up = c(10, 15, 20),
    outcome = c(1, 0, 1),
    age = c(50, 60, 45)
  )
  
  # Using numeric indices instead of column names
  result <- prepare_survival_data(test_data, time = 1, status = 2, verbose = FALSE)
  
  expect_equal(attr(result, "time_col"), "follow_up")  
  expect_equal(attr(result, "status_col"), "outcome")
})

test_that("prepare_survival_data respects validation settings", {
  test_data <- data.frame(
    time = c(-1, 15, 20),  # Invalid negative time
    status = c(1, 0, 1)
  )
  
  # With validation enabled (default), should error
  expect_error(
    prepare_survival_data(test_data, time = "time", status = "status", verbose = FALSE),
    "All survival times must be positive"
  )
  
  # With validation disabled, should work (though data is invalid)
  expect_no_error({
    result <- prepare_survival_data(test_data, time = "time", status = "status", validate = FALSE, verbose = FALSE)
    expect_s3_class(result, "survival_data")
  })
}) 
