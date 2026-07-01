## real-dataset stress checks ----------------------------------------------

quiet_short_chain_stan_warnings_local <- function(expr) {
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

test_that("real survival datasets work across fit/complete/export", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  if (!tolower(Sys.getenv("BIRD_REAL_DATASET_STRESS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_REAL_DATASET_STRESS=true to run real-dataset stress tests")
  }
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (stress-style real dataset test)")

  withr::local_options(list(bird.stan_backend = "precompiled"))

  datasets <- list(
    ovarian = list(data = survival::ovarian, time = "futime", status = "fustat"),
    veteran = list(data = survival::veteran, time = "time", status = "status"),
    aml = list(data = survival::aml, time = "time", status = "status")
  )

  td <- withr::local_tempdir()

  for (nm in names(datasets)) {
    cfg <- datasets[[nm]]
    dat <- cfg$data

    fit <- quiet_short_chain_stan_warnings_local(
      impute(
        dat,
        time = cfg$time,
        status = cfg$status,
        distribution = "weibull",
        n_imputations = 2,
        mcmc_options = list(
          iter_warmup = 50,
          iter_sampling = 50,
          chains = 2,
          adapt_delta = 0.9,
          max_treedepth = 10
        ),
        verbose = FALSE
      )
    )

    d1 <- complete(fit, dataset = 1)
    expect_equal(nrow(d1), nrow(dat))
    expect_true(all(c("time", "imputed_time", cfg$status, "original_status", "was_censored") %in% names(d1)))

    cens_idx <- which(d1$was_censored)
    if (length(cens_idx) > 0) {
      expect_true(all(d1$imputed_time[cens_idx] >= d1$time[cens_idx], na.rm = TRUE))
    }

    out_path <- export(
      fit,
      file.path(td, paste0("stress_", nm)),
      format = "csv",
      datasets = "first"
    )
    expect_length(out_path, 1)
    expect_true(file.exists(out_path))
    exported <- utils::read.csv(out_path, stringsAsFactors = FALSE)
    expect_equal(nrow(exported), nrow(dat))
    expect_true(all(c("time", "imputed_time", cfg$status, "original_status", "was_censored") %in% names(exported)))
  }
})

test_that("real datasets support grouped and nonparametric paths with export", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  if (!tolower(Sys.getenv("BIRD_REAL_DATASET_STRESS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_REAL_DATASET_STRESS=true to run real-dataset stress tests")
  }
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (stress-style real dataset test)")

  withr::local_options(list(bird.stan_backend = "precompiled"))
  td <- withr::local_tempdir()

  # Grouped parametric path on veteran
  fit_grp <- quiet_short_chain_stan_warnings_local(
    impute(
      survival::veteran,
      time = "time",
      status = "status",
      groups = "trt",
      distribution = "weibull",
      n_imputations = 2,
      mcmc_options = list(
        iter_warmup = 50,
        iter_sampling = 50,
        chains = 2,
        adapt_delta = 0.9,
        max_treedepth = 10
      ),
      verbose = FALSE
    )
  )

  grp_d1 <- complete(fit_grp, dataset = 1, groups = "combined")
  expect_equal(nrow(grp_d1), nrow(survival::veteran))
  expect_true("trt" %in% names(grp_d1))

  grp_out <- export(
    fit_grp,
    file.path(td, "stress_veteran_group"),
    format = "csv",
    groups = "combined",
    datasets = "first"
  )
  expect_length(grp_out, 1)
  expect_true(file.exists(grp_out))

  # Nonparametric path on ovarian (small and fast enough for stress check)
  fit_np <- impute(
    survival::ovarian,
    model = "nonparametric",
    time = "futime",
    status = "fustat",
    n_imputations = 2,
    mcmc = list(nburn = 20, nsave = 20, nskip = 1, ndisplay = 20),
    verbose = FALSE
  )

  np_d1 <- complete(fit_np, dataset = 1)
  expect_equal(nrow(np_d1), nrow(survival::ovarian))
  expect_true(all(c("time", "imputed_time", "fustat", "original_status", "was_censored") %in% names(np_d1)))

  np_out <- export(
    fit_np,
    file.path(td, "stress_ovarian_np"),
    format = "csv",
    datasets = "first"
  )
  expect_length(np_out, 1)
  expect_true(file.exists(np_out))
})

test_that("optional heavy real-dataset stress run", {
  skip_on_cran()
  skip_if_not_installed("survival")
  skip_if_not_installed("rstan")
  if (Sys.getenv("CI") == "true") skip("Skipping on CI (heavy stress test)")
  if (!tolower(Sys.getenv("BIRD_HEAVY_STRESS", "false")) %in% c("1", "true", "yes")) {
    skip("Set BIRD_HEAVY_STRESS=true to run heavy stress test")
  }

  withr::local_options(list(bird.stan_backend = "precompiled"))
  td <- withr::local_tempdir()

  # Larger real dataset from survival package
  dat <- survival::colon
  dat <- dat[is.finite(dat$time) & is.finite(dat$status), , drop = FALSE]

  fit <- impute(
    dat,
    time = "time",
    status = "status",
    distribution = "weibull",
    n_imputations = 10,
    mcmc_options = list(
      chains = 4,
      iter_warmup = 500,
      iter_sampling = 500,
      adapt_delta = 0.9,
      max_treedepth = 12
    ),
    verbose = TRUE
  )

  d1 <- complete(fit, dataset = 1)
  expect_equal(nrow(d1), nrow(dat))
  expect_true(all(c("time", "imputed_time", "status", "original_status", "was_censored") %in% names(d1)))

  out_path <- export(
    fit,
    file.path(td, "stress_colon_heavy"),
    format = "csv",
    datasets = "first"
  )
  expect_length(out_path, 1)
  expect_true(file.exists(out_path))
})
