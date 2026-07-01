## export ---------------------------------------------------------------

#' Tests for export functionality without context()

quiet_short_chain_stan_warnings <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      msg <- conditionMessage(w)
      is_expected_short_chain_diag <- grepl(
        "largest R-hat|Bulk Effective Samples Size|Tail Effective Samples Size",
        msg,
        ignore.case = TRUE
      )
      if (is_expected_short_chain_diag) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

make_mock_imputation_for_export <- function() {
  original <- data.frame(
    time = c(10, 15, 20, 25),
    status = c(1L, 0L, 1L, 0L),
    arm = c("A", "A", "B", "B")
  )

  dataset1 <- data.frame(
    time = c(10, 18, 20, 30),
    original_time = original$time,
    original_status = original$status,
    status = c(1L, 1L, 1L, 1L),
    was_censored = original$status == 0L,
    dataset_id = 1L,
    arm = original$arm
  )

  dataset2 <- data.frame(
    time = c(10, 19, 20, 33),
    original_time = original$time,
    original_status = original$status,
    status = c(1L, 1L, 1L, 1L),
    was_censored = original$status == 0L,
    dataset_id = 2L,
    arm = original$arm
  )

  structure(
    list(
      original_data = original,
      time_col = "time",
      status_col = "status",
      imputed_datasets = list(dataset1, dataset2)
    ),
    class = "bayesian_imputation"
  )
}

test_that("export saves mock single-model datasets without fitting", {
  td <- withr::local_tempdir()
  res <- make_mock_imputation_for_export()

  csv_paths <- NULL
  utils::capture.output({
    csv_paths <- export(res, file.path(td, "mock_csv"), format = "csv")
  })
  expect_length(csv_paths, 2)
  expect_true(all(file.exists(csv_paths)))

  csv1 <- utils::read.csv(csv_paths[[1]], stringsAsFactors = FALSE)
  expect_true(all(c(".imp", "dataset_id", "time", "imputed_time", "status", "original_status", "was_censored") %in% names(csv1)))
  expect_equal(nrow(csv1), nrow(res$original_data))
  expect_equal(unique(csv1$.imp), 1)

  rds_path <- NULL
  utils::capture.output({
    rds_path <- export(res, file.path(td, "mock_rds"), format = "rds", datasets = "first")
  })
  expect_length(rds_path, 1)
  expect_true(file.exists(rds_path))
  rds1 <- readRDS(rds_path)
  expect_true(is.data.frame(rds1))
  expect_equal(unique(rds1$.imp), 1)
})

test_that("export saves mock grouped datasets without fitting", {
  td <- withr::local_tempdir()
  group_a <- make_mock_imputation_for_export()
  group_b <- make_mock_imputation_for_export()
  group_a$original_data$grp <- "A"
  group_b$original_data$grp <- "B"
  for (i in seq_along(group_a$imputed_datasets)) {
    group_a$imputed_datasets[[i]]$grp <- "A"
    group_b$imputed_datasets[[i]]$grp <- "B"
  }

  res <- structure(
    list(
      group_names = c("A", "B"),
      group_results = list(A = group_a, B = group_b),
      group_errors = list(),
      groups = "grp"
    ),
    class = c("bayesian_imputation_groups", "bayesian_imputation")
  )

  path <- NULL
  utils::capture.output({
    path <- export(res, file.path(td, "mock_group"), format = "csv", datasets = "first")
  })
  expect_length(path, 1)
  expect_true(file.exists(path))

  dat <- utils::read.csv(path, stringsAsFactors = FALSE)
  expect_equal(nrow(dat), nrow(group_a$original_data) + nrow(group_b$original_data))
  expect_true("grp" %in% names(dat))
  expect_equal(sort(unique(dat$grp)), c("A", "B"))
  expect_equal(unique(dat$.imp), 1)
})

test_that("parametric export saves expected CSVs (all and first)", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  if (!tolower(Sys.getenv("BIRD_FITTING_TESTS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_FITTING_TESTS=true to run fitting-based export test")
  }
  # Skip entirely on CI to avoid C++ toolchain variability and compilation time
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (requires Stan compilation)")

  td <- withr::local_tempdir()
  file_base_all <- file.path(td, "weibull_results_single")
  file_base_first <- file.path(td, "weibull_results_first")

  set.seed(100)
  res <- quiet_short_chain_stan_warnings(
    bayesian_impute(
      survival::lung,
      n_imputations = 2,
      distribution = "weibull",
      mcmc_options = list(iter_warmup = 50, iter_sampling = 50, chains = 2, adapt_delta = 0.9, max_treedepth = 10),
      verbose = FALSE
    )
  )

  # Export all datasets to CSV (expects two files, one per dataset)
  paths_all <- export(res, file_base_all, format = "csv")
  expect_length(paths_all, 2)
  for (p in paths_all) {
    expect_true(file.exists(p))
    dat <- utils::read.csv(p, stringsAsFactors = FALSE)
    expect_true(all(c("time", "imputed_time", "status", "original_status", "was_censored", ".imp") %in% names(dat)))
    expect_equal(nrow(dat), nrow(res$original_data))
  }

  # Export first dataset to CSV (single file)
  path_first <- export(res, file_base_first, format = "csv", datasets = "first")
  expect_length(path_first, 1)
  expect_true(file.exists(path_first))
  dat1 <- utils::read.csv(path_first, stringsAsFactors = FALSE)
  expect_true(all(c("time", "imputed_time", "status", "original_status", "was_censored") %in% names(dat1)))
  expect_equal(nrow(dat1), nrow(res$original_data))
})

test_that("group export (combined) saves expected CSV", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  if (!tolower(Sys.getenv("BIRD_FITTING_TESTS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_FITTING_TESTS=true to run fitting-based export test")
  }
  # Skip entirely on CI to avoid C++ toolchain variability and compilation time
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (requires Stan compilation)")

  td <- withr::local_tempdir()
  file_base <- file.path(td, "weibull_results_groups")

  set.seed(101)
  resg <- quiet_short_chain_stan_warnings(
    bayesian_impute(
      survival::lung,
      groups = "sex",
      n_imputations = 2,
      distribution = "weibull",
      mcmc_options = list(iter_warmup = 50, iter_sampling = 50, chains = 2, adapt_delta = 0.9, max_treedepth = 10),
      verbose = FALSE
    )
  )

  # Combined export
  paths <- export(resg, file_base, format = "csv", groups = "combined")
  expect_length(paths, 1)
  expect_true(file.exists(paths))
  dat <- utils::read.csv(paths, stringsAsFactors = FALSE)
  # Expect group column preserved
  expect_true("sex" %in% names(dat))
  # Long format has stacked imputations when exporting all; allow any >= nrow
  expect_true(nrow(dat) >= nrow(resg$original_data))
})

test_that("nonparametric export saves expected CSV (first)", {
  skip_on_cran()
  skip_if_not_installed("survival")
  if (!tolower(Sys.getenv("BIRD_FITTING_TESTS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_FITTING_TESTS=true to run fitting-based export test")
  }

  td <- withr::local_tempdir()
  file_base <- file.path(td, "np_results_first")

  set.seed(102)
  res_np <- bayes_np_impute(
    survival::lung,
    n_imputations = 2,
    mcmc = list(nburn = 20, nsave = 20, nskip = 1, ndisplay = 20),
    verbose = FALSE
  )

  path_first <- export(res_np, file_base, format = "csv", datasets = "first")
  expect_length(path_first, 1)
  expect_true(file.exists(path_first))
  dat <- utils::read.csv(path_first, stringsAsFactors = FALSE)
  expect_true(all(c("time", "imputed_time", "status", "original_status", "was_censored") %in% names(dat)))
  expect_equal(nrow(dat), nrow(res_np$original_data))
})
