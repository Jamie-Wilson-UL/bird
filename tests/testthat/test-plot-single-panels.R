test_that("single-panel histogram/survival/boxplot plots work", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  base_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8),
    status = c(1, 1, 1, 0, 1, 0, 1, 1)
  )

  d1 <- data.frame(time = c(1, 2, 3, 4.8, 5, 6.7, 7, 8), status = rep(1L, 8))
  d2 <- data.frame(time = c(1, 2, 3, 5.2, 5, 7.1, 7, 8), status = rep(1L, 8))

  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1, d2),
      model_info = list(time_unit = "days")
    ),
    class = "bayesian_imputation"
  )

  p_hist <- plot.bayesian_imputation(mock_object, type = "histogram", dataset_id = 2)
  expect_s3_class(p_hist, "ggplot")
  expect_lt(p_hist$layers[[1]]$stat_params$bins, 30)

  p_hist_alias <- plot.bayesian_imputation(mock_object, type = "hist", dataset_id = 2)
  expect_s3_class(p_hist_alias, "ggplot")

  p_surv_main <- plot.bayesian_imputation(mock_object, type = "survival", dataset_id = 2)
  expect_s3_class(p_surv_main, "ggplot")
  expect_false(any(vapply(p_surv_main$layers, function(layer) identical(layer$aes_params$linetype, "dotted"), logical(1))))
  expect_equal(p_surv_main$layers[[1]]$data$time[1], 0)
  expect_equal(p_surv_main$layers[[1]]$data$survival[1], 1)

  p_surv_with_median <- plot.bayesian_imputation(
    mock_object,
    type = "survival",
    dataset_id = 2,
    show_median = TRUE
  )
  expect_s3_class(p_surv_with_median, "ggplot")
  expect_true(any(vapply(p_surv_with_median$layers, function(layer) identical(layer$aes_params$linetype, "dotted"), logical(1))))

  p_hist_multi <- plot.bayesian_imputation(mock_object, type = "histogram", dataset_id = c(1, 2))
  expect_s3_class(p_hist_multi, "ggplot")
  expect_true(".dataset" %in% names(p_hist_multi$data))
  expect_lt(p_hist_multi$layers[[1]]$stat_params$bins, 30)

  p_surv_multi <- plot.bayesian_imputation(mock_object, type = "survival", dataset_id = c(1, 2))
  expect_s3_class(p_surv_multi, "ggplot")
  expect_gte(length(p_surv_multi$layers), 1)

  p_box <- plot.bayesian_imputation(mock_object, type = "boxplot", dataset_id = 2)
  expect_s3_class(p_box, "ggplot")

  p_box_alias <- plot.bayesian_imputation(mock_object, type = "box", dataset_id = 2)
  expect_s3_class(p_box_alias, "ggplot")

  p_box_multi <- plot.bayesian_imputation(mock_object, type = "boxplot", dataset_id = c(1, 2))
  expect_s3_class(p_box_multi, "ggplot")
  expect_true(".dataset" %in% names(p_box_multi$data))
})

test_that("plotting helpers build baseline survival data and adaptive bins", {
  skip_if_not_installed("survival")

  km <- survival::survfit(
    survival::Surv(time, status) ~ 1,
    data = data.frame(time = c(4, 7, 10), status = c(1, 1, 1))
  )

  curve_data <- km_survival_plot_data(km, group = "A")
  expect_equal(curve_data$time[1], 0)
  expect_equal(curve_data$survival[1], 1)
  expect_equal(curve_data$group[1], "A")

  expect_equal(adaptive_histogram_bins(rep(5, 10)), 1)
  expect_equal(adaptive_histogram_bins(c(1, 2)), 2)
  expect_lt(adaptive_histogram_bins(c(1:21, 80)), 30)
})

test_that("plot labels use fitted object time units", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  base_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 1, 0, 1, 0, 1)
  )
  d1 <- data.frame(time = c(1, 2, 3.5, 4, 5.5, 6), status = rep(1L, 6))

  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1),
      model_info = list(time_unit = "weeks")
    ),
    class = "bayesian_imputation"
  )

  p_hist <- plot.bayesian_imputation(mock_object, type = "histogram", dataset_id = 1)
  expect_equal(p_hist$labels$x, "Time (weeks)")

  p_surv <- plot.bayesian_imputation(mock_object, type = "survival", dataset_id = 1)
  expect_equal(p_surv$labels$x, "Time (weeks)")

  p_box <- plot.bayesian_imputation(mock_object, type = "boxplot", dataset_id = 1)
  expect_equal(p_box$labels$y, "Time (weeks)")

  group_object <- structure(
    list(
      group_names = "A",
      group_results = list(A = mock_object)
    ),
    class = "bayesian_imputation_groups"
  )
  expect_equal(time_axis_label(group_object), "Time (weeks)")

  model_object <- structure(
    list(
      model_names = "weibull",
      model_results = list(weibull = mock_object)
    ),
    class = "bird_model_comparison"
  )
  expect_equal(time_axis_label(model_object, "Survival time"), "Survival time (weeks)")
})

test_that("single-panel plot types validate dataset_id", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  base_data <- data.frame(
    time = c(1, 2, 3, 4),
    status = c(1, 0, 1, 1)
  )
  d1 <- data.frame(time = c(1, 2.5, 3, 4), status = rep(1L, 4))

  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1),
      model_info = list(time_unit = "days")
    ),
    class = "bayesian_imputation"
  )

  expect_error(
    plot.bayesian_imputation(mock_object, type = "histogram", dataset_id = 0),
    "dataset_id must be between 1 and 1"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "survival", dataset_id = 2),
    "dataset_id must be between 1 and 1"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "survival", dataset_id = 1.5),
    "dataset_id values must be integers"
  )
  expect_error(
    plot.bayesian_imputation(mock_object, type = "boxplot", dataset_id = 3),
    "dataset_id must be between 1 and 1"
  )
})

test_that("single-model palette overrides are accepted", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")

  base_data <- data.frame(
    time = c(1, 2, 3, 4, 5, 6),
    status = c(1, 1, 0, 1, 0, 1)
  )
  d1 <- data.frame(time = c(1, 2, 3.8, 4, 5.6, 6), status = rep(1L, 6))
  d2 <- data.frame(time = c(1, 2, 4.1, 4, 5.9, 6), status = rep(1L, 6))

  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(d1, d2),
      model_info = list(time_unit = "days")
    ),
    class = "bayesian_imputation"
  )

  p1 <- plot.bayesian_imputation(
    mock_object,
    type = "histogram",
    dataset_id = 1,
    palette = c(hist_fill = "#4daf4a")
  )
  expect_s3_class(p1, "ggplot")

  p2 <- plot.bayesian_imputation(
    mock_object,
    type = "survival",
    dataset_id = 1,
    palette = c(survival_curve = "#e41a1c")
  )
  expect_s3_class(p2, "ggplot")

  p3 <- plot.bayesian_imputation(
    mock_object,
    type = "density",
    dataset_id = 1,
    color = "#377eb8"
  )
  expect_s3_class(p3, "ggplot")

  p4 <- plot.bayesian_imputation(
    mock_object,
    type = "boxplot",
    dataset_id = 1,
    color = "#984ea3"
  )
  expect_s3_class(p4, "ggplot")
})
