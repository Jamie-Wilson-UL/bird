test_that("get_default_priors returns correct structure", {
  priors <- get_default_priors("weibull")
  
  expect_true(is.list(priors))
  expect_true(all(c(
    "prior_family",
    "scale_prior_shape", "scale_prior_rate",
    "mu_log_shape", "sd_log_shape",
    "mu_log_scale", "sd_log_scale"
  ) %in% names(priors)))
  
  # Check reasonable default values
  expect_equal(priors$prior_family, "gamma_scale")
  expect_true(priors$scale_prior_shape > 0)
  expect_true(priors$scale_prior_rate > 0)
  expect_equal(priors$mu_log_shape, 0)
  expect_true(priors$sd_log_shape > 0)
  expect_equal(priors$mu_log_scale, 0)
  expect_true(priors$sd_log_scale > 0)
})

test_that("get_default_priors handles unknown distribution", {
  expect_error(
    get_default_priors("unknown"),
    "Unknown distribution: unknown"
  )
})

test_that("get_adaptive_priors returns matched Weibull Gamma scale prior", {
  priors <- get_adaptive_priors("weibull", c(30, 45, 60, 90, 120, 180))

  lognormal_mean <- function(mu, sigma) exp(mu + 0.5 * sigma^2)
  lognormal_var <- function(mu, sigma) {
    (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
  }

  scale_mean <- lognormal_mean(priors$mu_log_scale, priors$sd_log_scale)
  scale_var <- lognormal_var(priors$mu_log_scale, priors$sd_log_scale)

  expect_equal(priors$prior_family, "gamma_scale")
  expect_equal(priors$scale_prior_shape, scale_mean^2 / scale_var)
  expect_equal(priors$scale_prior_rate, scale_mean / scale_var)
})

test_that("normalize_weibull_priors preserves custom values and fills defaults", {
  priors <- normalize_weibull_priors(list(
    prior_family = "gamma_scale",
    scale_prior_shape = 3,
    scale_prior_rate = 0.02
  ))

  expect_equal(priors$prior_family, "gamma_scale")
  expect_equal(priors$scale_prior_shape, 3)
  expect_equal(priors$scale_prior_rate, 0.02)
  expect_true(all(c(
    "mu_log_shape", "sd_log_shape",
    "mu_log_scale", "sd_log_scale"
  ) %in% names(priors)))
})

test_that("stan_models_ready function works", {
  result <- stan_models_ready()
  expect_true(is.logical(result))
})
