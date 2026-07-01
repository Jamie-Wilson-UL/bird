test_that("bayes_np_impute runs and returns expected structure", {
  skip_on_cran()
  skip_if_not_installed("survival")
  if (!tolower(Sys.getenv("BIRD_FITTING_TESTS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_FITTING_TESTS=true to run fitting-based nonparametric test")
  }

  set.seed(2025)
  result <- bayes_np_impute(
    survival::lung,
    n_imputations = 2,
    mcmc = list(nburn = 10, nsave = 10, nskip = 1, ndisplay = 10),
    verbose = FALSE
  )

  expect_s3_class(result, "bayesian_imputation")
  expect_equal(result$model_info$distribution, "nonparametric_lddp")
  expect_equal(length(result$imputed_datasets), 2)
  expect_true(is.list(result$fit))
  expect_true(is.matrix(result$fit$save.state$ysave))
  expect_true(result$diagnostics$convergence_ok)
  
  # Clean up any dppackage output files created during the test
  dppackage_files <- list.files(pattern = "dppackage.*\\.out$", full.names = TRUE)
  if (length(dppackage_files) > 0) {
    unlink(dppackage_files)
  }
})
