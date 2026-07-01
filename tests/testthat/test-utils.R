test_that("format_time works correctly", {
  expect_equal(format_time(30), "30.0 seconds")
  expect_equal(format_time(90), "1.5 minutes")
  expect_equal(format_time(3900), "1.1 hours")
})

test_that("safe_eval handles errors appropriately", {
  # Should work normally
  expect_equal(safe_eval(2 + 2), 4)
  
  # Should catch errors
  expect_error(
    safe_eval(stop("test error"), "Custom error"),
    "Custom error: test error"
  )
})

test_that("bind_rows_fill handles missing columns", {
  a <- data.frame(x = 1:2, y = c("a", "b"))
  b <- data.frame(y = "c", z = TRUE)

  out <- bind_rows_fill(list(a, b))

  expect_equal(names(out), c("x", "y", "z"))
  expect_equal(nrow(out), 3)
  expect_true(is.na(out$x[3]))
  expect_true(is.na(out$z[1]))
  expect_true(out$z[3])
})

test_that("package environment exists", {
  expect_true(exists(".bird_env", envir = asNamespace("bird")))
  expect_true(is.environment(.bird_env))
}) 
