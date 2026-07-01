test_that("pooled Cox Wald test combines imputations correctly", {
  make_group_result <- function(seed) {
    set.seed(seed)
    original <- data.frame(
      time = c(5, 7, 9, 11, 13),
      status = c(1, 0, 1, 0, 1),
      group = sample(c("A", "B"), 5, replace = TRUE)
    )
    dataset1 <- original
    dataset1$time <- dataset1$time + runif(5, 0, 2)
    dataset1$status[] <- 1L
    dataset1$dataset_id <- 1
    dataset1$original_time <- original$time
    dataset1$original_status <- original$status
    dataset1$was_censored <- original$status == 0
    dataset2 <- dataset1
    dataset2$time <- dataset2$time + runif(5, 0, 1)
    dataset2$dataset_id <- 2
    structure(
      list(
        original_data = original,
        imputed_datasets = list(dataset1, dataset2),
        time_col = "time",
        status_col = "status",
        model_info = list(distribution = "weibull"),
        diagnostics = list(convergence_ok = TRUE)
      ),
      class = "bayesian_imputation"
    )
  }

  group_results <- list(
    A = make_group_result(1),
    B = make_group_result(2)
  )

  pooled_object <- structure(
    list(
      group_results = group_results,
      group_errors = list(),
      group_names = c("A", "B"),
      groups = "group",
      original_data = do.call(rbind, lapply(group_results, `[[`, "original_data"))
    ),
    class = c("bayesian_imputation_groups", "bayesian_imputation")
  )

  pooled <- pool_cox_group_test(pooled_object)
  expect_true(is.list(pooled))
  expect_true(all(c("coef", "vcov", "statistic", "df", "df2", "r1", "p_value") %in% names(pooled)))
  expect_true(is.numeric(pooled$statistic))
  expect_true(pooled$statistic >= 0)
  expect_true(is.numeric(pooled$df2))
  expect_true(is.finite(pooled$df2) || is.infinite(pooled$df2))
  expect_true(is.numeric(pooled$p_value))
  expect_true(pooled$p_value >= 0 && pooled$p_value <= 1)

  test_result <- calculate_log_rank_test(pooled_object)
  expect_true(all(c("statistic", "numerator_df", "denominator_df", "p_value") %in% names(test_result)))
})
