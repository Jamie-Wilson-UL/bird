test_that("validate_survival_data works with valid data", {
  # Create test data
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1),
    x = 1:5
  )
  
  result <- validate_survival_data(test_data, "time", "status")
  
  expect_equal(result$n_total, 5)
  expect_equal(result$n_observed, 3)
  expect_equal(result$n_censored, 2)
  expect_equal(result$censoring_rate, 0.4)
  expect_equal(result$observed_times, c(10, 20, 30))
  expect_equal(result$censored_times, c(15, 25))
})

test_that("validate_survival_data handles logical status", {
  test_data <- data.frame(
    time = c(10, 15, 20),
    status = c(TRUE, FALSE, TRUE)
  )
  
  result <- validate_survival_data(test_data, "time", "status")
  
  expect_equal(result$n_observed, 2)
  expect_equal(result$n_censored, 1)
})

test_that("validate_survival_data catches missing columns", {
  test_data <- data.frame(time = c(10, 15, 20))
  
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "Status variable 'status' not found in data"
  )
  
  expect_error(
    validate_survival_data(test_data, "duration", "status"),
    "Time variable 'duration' not found in data"
  )
})

test_that("validate_survival_data catches invalid time values", {
  # Negative times
  test_data <- data.frame(
    time = c(-1, 15, 20),
    status = c(1, 0, 1)
  )
  
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "All survival times must be positive"
  )
  
  # Zero times
  test_data$time[1] <- 0
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "All survival times must be positive"
  )
  
  # Missing times
  test_data$time[1] <- NA
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "Missing values in time variable not allowed"
  )
})

test_that("validate_survival_data catches invalid status values", {
  # Invalid status values
  test_data <- data.frame(
    time = c(10, 15, 20, 25),
    status = c(1, 2, 1, 99)
  )
  
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "Status variable must contain only"
  )
  
  # Missing status
  test_data$status[2] <- NA
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "Missing values in status variable not allowed"
  )
})

test_that("validate_survival_data catches insufficient events", {
  # Only one event
  test_data <- data.frame(
    time = c(10, 15, 20),
    status = c(1, 0, 0)
  )
  
  expect_error(
    validate_survival_data(test_data, "time", "status"),
    "At least 2 observed events required for estimation"
  )
})

test_that("validate_survival_data warns about no censoring", {
  # No censored observations
  test_data <- data.frame(
    time = c(10, 15, 20),
    status = c(1, 1, 1)
  )
  
  expect_warning(
    validate_survival_data(test_data, "time", "status"),
    "No censored observations found"
  )
})

test_that("validate_imputation_params works correctly", {
  result <- validate_imputation_params(50, "weibull", NULL)
  
  expect_equal(result$n_imputations, 50)
  expect_equal(result$distribution, "weibull")
  expect_true(is.list(result$mcmc_options))
  expect_equal(result$mcmc_options$chains, 4)
})

test_that("validate_imputation_params catches invalid inputs", {
  # Invalid n_imputations
  expect_error(
    validate_imputation_params(0, "weibull", NULL),
    "Assertion on 'n_imputations' failed"
  )
  
  expect_error(
    validate_imputation_params(1001, "weibull", NULL),
    "Assertion on 'n_imputations' failed"
  )
  
  # Invalid distribution
  expect_error(
    validate_imputation_params(50, "normal", NULL),
    "Distribution 'normal' not supported"
  )
})

test_that("validate_imputation_params warns about low imputations", {
  expect_warning(
    validate_imputation_params(5, "weibull", NULL),
    "Using fewer than 10 imputations may give unreliable results"
  )
})

test_that("prepare_stan_data formats correctly", {
  # Create test data
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )
  
  priors <- get_default_priors("weibull")
  
  stan_data <- prepare_stan_data(test_data, "time", "status", "weibull", priors)
  
  expect_equal(stan_data$n_obs, 3)
  expect_equal(stan_data$n_cens, 2)
  expect_equal(length(stan_data$t_obs), 3)
  expect_equal(length(stan_data$t_cens), 2)
  expect_equal(stan_data$prior_family, 1L)
  expect_true(all(c(
    "scale_prior_shape", "scale_prior_rate",
    "mu_log_shape", "sd_log_shape",
    "mu_log_scale", "sd_log_scale",
    "shape_upper", "scale_upper"
  ) %in% names(stan_data)))
})

test_that("prepare_stan_data supports Weibull lognormal prior family", {
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )

  priors <- get_default_priors("weibull")
  priors$prior_family <- "lognormal"

  stan_data <- prepare_stan_data(test_data, "time", "status", "weibull", priors)

  expect_equal(stan_data$prior_family, 2L)
  expect_true(all(c(
    "scale_prior_shape", "scale_prior_rate",
    "mu_log_shape", "sd_log_shape",
    "mu_log_scale", "sd_log_scale"
  ) %in% names(stan_data)))
})

test_that("prepare_stan_data fills unused Weibull prior family defaults", {
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )

  priors <- list(
    prior_family = "gamma_scale",
    scale_prior_shape = 3,
    scale_prior_rate = 0.02
  )

  stan_data <- prepare_stan_data(test_data, "time", "status", "weibull", priors)

  expect_equal(stan_data$prior_family, 1L)
  expect_equal(stan_data$scale_prior_rate, 0.02)
  expect_true(all(c(
    "mu_log_shape", "sd_log_shape",
    "mu_log_scale", "sd_log_scale"
  ) %in% names(stan_data)))
})

test_that("prepare_stan_data works with exponential distribution", {
  # Create test data
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )
  
  priors <- get_default_priors("exponential")
  
  stan_data <- prepare_stan_data(test_data, "time", "status", "exponential", priors)
  
  expect_equal(stan_data$n_obs, 3)
  expect_equal(stan_data$n_cens, 2)
  expect_true(all(c("rate_prior_shape", "rate_prior_rate") %in% names(stan_data)))
})

test_that("prepare_stan_data works with lognormal distribution", {
  # Create test data
  test_data <- data.frame(
    time = c(10, 15, 20, 25, 30),
    status = c(1, 0, 1, 0, 1)
  )
  
  priors <- get_default_priors("lognormal")
  
  stan_data <- prepare_stan_data(test_data, "time", "status", "lognormal", priors)
  
  expect_equal(stan_data$n_obs, 3)
  expect_equal(stan_data$n_cens, 2)
  expect_true(all(c("mu_prior_mean", "mu_prior_sd", "sigma_prior_sd") %in% names(stan_data)))
}) 
