test_that("Stan backend auto mode matches model availability", {
  withr::local_options(list(bird.stan_backend = "auto"))
  expected <- if (is.null(get_precompiled_stan_model("weibull"))) "runtime" else "precompiled"
  expect_identical(select_stan_backend("weibull"), expected)
})

test_that("Stan backend option can force runtime", {
  withr::local_options(list(bird.stan_backend = "runtime"))
  expect_identical(select_stan_backend("weibull"), "runtime")
})

test_that("Invalid Stan backend option falls back safely", {
  withr::local_options(list(bird.stan_backend = "not-a-backend"))
  expected <- if (is.null(get_precompiled_stan_model("weibull"))) "runtime" else "precompiled"
  expect_warning(
    backend <- select_stan_backend("weibull"),
    "Invalid option bird.stan_backend"
  )
  expect_identical(backend, expected)
})

test_that("Legacy backend option key is still honored", {
  withr::local_options(list(bird.stan_backend = NULL, bayessurvival.stan_backend = "runtime"))
  expect_identical(select_stan_backend("weibull"), "runtime")
})
