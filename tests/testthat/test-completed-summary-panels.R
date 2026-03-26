make_mock_imputation_object <- function() {
  original_data <- data.frame(
    time = c(2, 4, 6, 8, 10, 12),
    status = c(1, 1, 0, 1, 0, 1)
  )

  d1 <- data.frame(time = c(2, 4, 7, 8, 11, 12), status = rep(1L, 6))
  d2 <- data.frame(time = c(2, 4, 7.5, 8, 10.5, 12), status = rep(1L, 6))

  structure(
    list(
      original_data = original_data,
      imputed_datasets = list(d1, d2),
      time_col = "time",
      status_col = "status"
    ),
    class = "bayesian_imputation"
  )
}

make_mock_group_object <- function() {
  g1 <- make_mock_imputation_object()
  g2 <- make_mock_imputation_object()
  g3 <- make_mock_imputation_object()

  # make groups slightly different
  g2$imputed_datasets[[1]]$time <- g2$imputed_datasets[[1]]$time + 1
  g3$imputed_datasets[[1]]$time <- g3$imputed_datasets[[1]]$time + 2

  structure(
    list(
      group_names = c("A", "B", "C"),
      group_results = list(A = g1, B = g2, C = g3),
      original_data = g1$original_data
    ),
    class = "bayesian_imputation_groups"
  )
}

make_mock_model_comparison_object <- function() {
  m1 <- make_mock_imputation_object()
  m2 <- make_mock_imputation_object()
  m3 <- make_mock_imputation_object()

  m2$imputed_datasets[[1]]$time <- m2$imputed_datasets[[1]]$time + 0.75
  m3$imputed_datasets[[1]]$time <- m3$imputed_datasets[[1]]$time + 1.5

  structure(
    list(
      model_names = c("weibull", "lognormal", "nonparametric"),
      model_results = list(weibull = m1, lognormal = m2, nonparametric = m3),
      original_data = m1$original_data,
      time_col = "time",
      status_col = "status",
      model_info = list(time_unit = "days")
    ),
    class = "bird_model_comparison"
  )
}

test_that("single completed summary accepts custom panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("gridExtra")

  fit <- make_mock_imputation_object()
  g <- plot.bayesian_imputation(
    fit,
    type = "completed_dataset_summary",
    dataset_id = 1,
    panels = c("hist", "boxplot")
  )

  expect_s3_class(g, "gtable")
})

test_that("single completed summary validates panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("gridExtra")

  fit <- make_mock_imputation_object()
  expect_error(
    plot.bayesian_imputation(
      fit,
      type = "completed_dataset_summary",
      dataset_id = 1,
      panels = c("hist", "banana")
    ),
    "Invalid `panels`"
  )
})

test_that("group completed summary supports custom panels with >2 groups", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("gridExtra")

  fit_grp <- make_mock_group_object()
  expect_warning(
    g <- plot.bayesian_imputation_groups(
      fit_grp,
      type = "completed_dataset_summary",
      dataset_id = 1,
      panels = c("hist", "survival", "boxplot")
    ),
    "visually crowded"
  )
  expect_s3_class(g, "gtable")
})

test_that("model comparison completed summary supports custom panels", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("gridExtra")

  cmp <- make_mock_model_comparison_object()
  expect_warning(
    g <- plot.bird_model_comparison(
      cmp,
      type = "completed_dataset_summary",
      dataset_id = 1,
      panels = c("hist", "survival", "boxplot")
    ),
    "visually crowded"
  )
  expect_s3_class(g, "gtable")
})

test_that("model comparison supports boxplots_comparison type", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("survival")
  skip_if_not_installed("gridExtra")

  cmp <- make_mock_model_comparison_object()
  g <- plot.bird_model_comparison(
    cmp,
    type = "boxplots_comparison",
    n_max = 2
  )
  expect_s3_class(g, "gtable")
})

test_that("model comparison supports density type", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("logspline")

  cmp <- make_mock_model_comparison_object()
  p <- plot.bird_model_comparison(
    cmp,
    type = "density",
    dataset_id = 1
  )
  expect_s3_class(p, "ggplot")
})
