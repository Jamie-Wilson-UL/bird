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

  p_hist_alias <- plot.bayesian_imputation(mock_object, type = "hist", dataset_id = 2)
  expect_s3_class(p_hist_alias, "ggplot")

  p_surv_main <- plot.bayesian_imputation(mock_object, type = "survival", dataset_id = 2)
  expect_s3_class(p_surv_main, "ggplot")

  p_hist_multi <- plot.bayesian_imputation(mock_object, type = "histogram", dataset_id = c(1, 2))
  expect_s3_class(p_hist_multi, "ggplot")
  expect_true(".dataset" %in% names(p_hist_multi$data))

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
