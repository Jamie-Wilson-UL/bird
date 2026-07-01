test_that("generate_complete_datasets produces monotone imputations", {
  set.seed(123)
  data <- data.frame(
    time = c(5, 10, 8, 6),
    status = c(1, 0, 0, 1)
  )
  posterior_imputations <- matrix(
    c(9, 11,
      10, 12,
      11, 13,
      12, 14),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
  )
  n_imputations <- 3
  completed <- generate_complete_datasets(data, "time", "status", posterior_imputations, n_imputations)
  expect_equal(length(completed), n_imputations)
  censored_idx <- which(data$status == 0)
  for (dataset in completed) {
    expect_true(all(dataset$status == 1))
    expect_true(all(dataset$time[censored_idx] >= data$time[censored_idx]))
    expect_true(all(dataset$original_time == data$time))
    expect_true(all(dataset$original_status == data$status))
  }
})

test_that("generate_complete_datasets_np standardizes no-censoring datasets", {
  data <- data.frame(
    time = c(5, 8, 12),
    status = c(1, 1, 1)
  )

  completed <- generate_complete_datasets_np(
    data,
    time_col = "time",
    status_col = "status",
    posterior_imputations = matrix(numeric(0), nrow = 4, ncol = 0),
    censored_indices = integer(0),
    n_imputations = 2
  )

  expect_equal(length(completed), 2)
  expect_true(all(c(
    "dataset_id", "original_status", "status",
    "was_censored", "time", "original_time"
  ) %in% names(completed[[1]])))
  expect_equal(completed[[1]]$dataset_id, rep(1L, nrow(data)))
  expect_false(any(completed[[1]]$was_censored))
  expect_equal(completed[[1]]$original_time, data$time)
  expect_equal(completed[[1]]$original_status, data$status)
})

test_that("calculate_imputation_summary handles censored observations", {
  base_data <- data.frame(
    time = c(5, 7, 10, 12),
    status = c(1, 0, 0, 1)
  )
  dataset1 <- base_data
  dataset1$time[c(2,3)] <- c(9, 11)
  dataset1$status[] <- 1L
  dataset1$original_time <- base_data$time
  dataset1$original_status <- base_data$status
  dataset1$was_censored <- base_data$status == 0
  dataset1$dataset_id <- 1
  dataset2 <- dataset1
  dataset2$time[c(2,3)] <- c(8, 13)
  dataset2$dataset_id <- 2
  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(dataset1, dataset2)
    ),
    class = "bayesian_imputation"
  )
  summary <- calculate_imputation_summary(mock_object)
  expect_true(summary$has_imputations)
  expect_equal(summary$n_censored, 2)
  expect_length(summary$original_censored_summary, 4)
  expect_length(summary$imputed_times_summary, 4)
  expect_length(summary$imputation_gains_summary, 4)
  expect_true(summary$imputation_gains_summary["min"] >= 0)
})

test_that("calculate_clinical_metrics returns coherent metrics", {
  skip_if_not_installed("survival")
  skip_if_not_installed("posterior")
  base_data <- data.frame(
    time = c(5, 7, 10, 12),
    status = c(1, 0, 0, 1)
  )
  dataset1 <- base_data
  dataset1$time <- c(5, 9, 11, 12)
  dataset1$status[] <- 1L
  dataset1$original_time <- base_data$time
  dataset1$original_status <- base_data$status
  dataset1$was_censored <- base_data$status == 0
  dataset1$dataset_id <- 1
  dataset2 <- dataset1
  dataset2$time <- c(5, 8, 12, 12)
  dataset2$dataset_id <- 2
  mock_object <- structure(
    list(
      original_data = base_data,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(dataset1, dataset2),
      model_info = list(
        distribution = "weibull",
        time_unit = "days"
      ),
      posterior_samples = posterior::as_draws_df(
        data.frame(shape = rnorm(100, 1.5, 0.1), scale = rnorm(100, 10, 0.5))
      )
    ),
    class = "bayesian_imputation"
  )
  metrics <- calculate_clinical_metrics(mock_object, time_points = c(6, 12))
  expect_true(is.list(metrics$original))
  expect_true(is.list(metrics$imputed))
  expect_equal(length(metrics$imputed$median_summary), 6)
  expect_equal(nrow(metrics$imputed$survival_summary), 4)
  expect_equal(metrics$time_points, c(6, 12))
  # interpretation field was removed in recent updates
})
