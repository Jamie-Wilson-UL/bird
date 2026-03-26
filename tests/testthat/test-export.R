## export ---------------------------------------------------------------

#' Tests for export functionality without context()
test_that("parametric export saves expected CSVs (all and first)", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  # Skip entirely on CI to avoid C++ toolchain variability and compilation time
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (requires Stan compilation)")

  td <- withr::local_tempdir()
  file_base_all <- file.path(td, "weibull_results_single")
  file_base_first <- file.path(td, "weibull_results_first")

  set.seed(100)
  res <- bayesian_impute(
    survival::lung,
    n_imputations = 2,
    distribution = "weibull",
    mcmc_options = list(iter_warmup = 50, iter_sampling = 50, chains = 2, adapt_delta = 0.9, max_treedepth = 10),
    verbose = FALSE
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
  # Skip entirely on CI to avoid C++ toolchain variability and compilation time
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (requires Stan compilation)")

  td <- withr::local_tempdir()
  file_base <- file.path(td, "weibull_results_groups")

  set.seed(101)
  resg <- bayesian_impute(
    survival::lung,
    groups = "sex",
    n_imputations = 2,
    distribution = "weibull",
    mcmc_options = list(iter_warmup = 50, iter_sampling = 50, chains = 2, adapt_delta = 0.9, max_treedepth = 10),
    verbose = FALSE
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
