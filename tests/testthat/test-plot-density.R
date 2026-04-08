test_that("plot_density_comparison supports dataset_id selection", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("logspline")
  
  base_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    status = c(1, 1, 1, 1, 1, 0, 0, 1, 1, 1)
  )
  
  d1 <- data.frame(time = c(1, 2, 3, 4, 5, 6.4, 7.1, 8, 9, 10), status = rep(1L, 10))
  d2 <- data.frame(time = c(1, 2, 3, 4, 5, 6.8, 7.4, 8, 9, 10), status = rep(1L, 10))
  d3 <- data.frame(time = c(1, 2, 3, 4, 5, 7.2, 7.8, 8, 9, 10), status = rep(1L, 10))
  
  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1, d2, d3),
      model_info = list(time_unit = "days")
    ),
    class = "bayesian_imputation"
  )
  
  p <- plot.bayesian_imputation(mock_object, type = "density", dataset_id = 3)
  expect_s3_class(p, "ggplot")
  expect_true(all(c("x", "y") %in% names(p$data)))
  expect_false("type" %in% names(p$data))

  p_multi <- plot.bayesian_imputation(mock_object, type = "density", dataset_id = c(1, 3))
  expect_s3_class(p_multi, "ggplot")
  expect_true(".dataset" %in% names(p_multi$data))
})

test_that("plot_density_comparison validates dataset_id", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("logspline")
  
  base_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 1, 1, 0, 1, 1)
  )
  d1 <- data.frame(time = c(1, 2, 3, 4.4, 5, 6), status = rep(1L, 6))
  d2 <- data.frame(time = c(1, 2, 3, 4.8, 5, 6), status = rep(1L, 6))
  
  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1, d2)
    ),
    class = "bayesian_imputation"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "density", dataset_id = 0),
    "dataset_id must be between 1 and 2"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "density", dataset_id = 3),
    "dataset_id must be between 1 and 2"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "density", dataset_id = c(1, 1.2)),
    "dataset_id values must be integers"
  )
})
