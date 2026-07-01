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

  p_default <- plot.bayesian_imputation(mock_object, type = "density")
  expect_equal(max(p_default$data$x), max(d1$time))
})

test_that("density helper matches default logspline without boundary plateau", {
  skip_if_not_installed("logspline")

  control_times <- c(1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23)
  grid <- density_plot_grid(control_times, n = 400)

  actual <- logspline_density_on_grid(control_times, grid)
  expected <- logspline::dlogspline(grid, logspline::logspline(control_times))

  expect_equal(grid[1], 0)
  expect_equal(tail(grid, 1), max(control_times))
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_lt(actual[1], max(actual) / 2)
})

test_that("density helper uses log scale automatically for larger positive survival data", {
  skip_if_not_installed("logspline")

  values <- rep(c(30, 60, 90, 180, 270, 365, 540, 730, 1000), times = c(8, 14, 20, 35, 42, 38, 25, 12, 6))
  grid <- density_plot_grid(values, n = 400)

  auto_density <- logspline_density_on_grid(values, grid)
  log_density <- logspline_density_on_grid(values, grid, density_scale = "log")
  time_density <- logspline_density_on_grid(values, grid, density_scale = "time")

  expect_equal(auto_density, log_density, tolerance = 1e-10)
  expect_equal(auto_density[1], 0)
  expect_gt(max(time_density), 0)
})

test_that("selected dataset labels distinguish ids from counts", {
  shared_ids <- c("6-MP" = 97L, control = 97L)
  different_ids <- c("6-MP" = 12L, control = 97L)
  default_ids <- resolve_comparison_dataset_ids(
    dataset_id = NULL,
    keys = c("6-MP", "control"),
    n_by_key = c(100, 100)
  )

  expect_equal(selected_dataset_label(shared_ids), "selected dataset 97")
  expect_equal(
    selected_dataset_label(different_ids),
    "selected datasets: 6-MP=12, control=97"
  )
  expect_equal(default_ids, c("6-MP" = 1L, control = 1L))
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
