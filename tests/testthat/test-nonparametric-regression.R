test_that("nonparametric fingerprint is stable under fixed seed", {
  skip_if_not_installed("survival")

  set.seed(20260211)
  res <- bayes_np_impute(
    survival::lung,
    n_imputations = 2,
    mcmc = list(nburn = 30, nsave = 40, nskip = 1, ndisplay = 40),
    verbose = FALSE
  )

  cens <- res$original_data$status == 0
  d1 <- res$imputed_datasets[[1]]
  d2 <- res$imputed_datasets[[2]]
  ysave <- res$fit$save.state$ysave
  thetasave <- res$fit$save.state$thetasave

  fingerprint <- c(
    mean1 = mean(d1$time[cens]),
    mean2 = mean(d2$time[cens]),
    q95_1 = as.numeric(stats::quantile(d1$time[cens], 0.95)),
    q95_2 = as.numeric(stats::quantile(d2$time[cens], 0.95)),
    median_ysave_cens = stats::median(as.numeric(ysave[, cens])),
    alpha_last = tail(thetasave[, "alpha"], 1)
  )

  expected <- c(
    mean1 = 692.666822,
    mean2 = 794.101559,
    q95_1 = 1827.134692,
    q95_2 = 1726.301755,
    median_ysave_cens = 6.367891,
    alpha_last = 1.159034
  )

  expect_equal(fingerprint, expected, tolerance = 5e-3)
})
