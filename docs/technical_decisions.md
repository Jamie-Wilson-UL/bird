# Technical Decisions/Design Rationale

This document describes key design choices, default settings and assumptions in `bird`.

---

## 1. Prior Specification

When you run a parametric model without supplying priors, the package first tries to learn reasonable hyperparameters from the observed event times. If there are too few usable observed events, it falls back to generic weakly informative settings; if the observed times are very homogeneous, it uses conservative spread safeguards. Nonparametric models use fixed defaults inherited from DPpackage.

### 1.1 Parametric Models: Default Data-Adaptive Priors 

The parametric models (`weibull`, `exponential`, `lognormal`) aim to let the data speak wherever possible. If you don't supply a `priors` argument:

1. The code pulls out the observed event times (the uncensored observations).
2. It calls `get_adaptive_priors(distribution, t_obs)`, which inspects the median, spread and scale of those times.
3. The function sets hyperparameters that center the prior on plausible parameter ranges for the dataset while keeping enough spread for the likelihood to dominate.
4. If there are fewer than 2 valid observed event times, or the observed median is invalid, it falls back to `get_default_priors(distribution)` which gives generic weakly informative settings.

The adaptive priors are empirical Bayes style defaults. They are not intended to encode strong prior knowledge, and they are not fitted as a separate model. They use simple robust summaries of the observed event times to put the sampler in the right numerical region, especially when the time scale is far from 1.

The individual distribution identities used below are standard. The exact default-prior recipe in `bird` is not claimed to be a named published prior. It is a pragmatic empirical-Bayes construction built from standard survival-distribution relationships:

- Weibull survival, quantile and median formulas are standard Weibull distribution results.
- The relationship between a Weibull survival time and an extreme-value/Gumbel distribution on log-time is the standard accelerated-failure-time representation of the Weibull model.
- The conversion between lognormal mean/variance and Gamma shape/rate is ordinary moment matching.
- The fixed clamps and fallback spreads are package safeguards. They avoid overly tight or extreme defaults in small or homogeneous datasets, but those exact constants are not canonical.

#### Weibull

Stan samples positive `shape` and `scale` parameters. By default, the Weibull model uses a mixed prior specification:

`shape ~ Lognormal(mu_log_shape, sd_log_shape)`

`scale ~ Gamma(scale_prior_shape, scale_prior_rate)`

This keeps the original lognormal shape prior while using the PI's published Gamma prior family for the Weibull scale parameter. The Gamma prior is specified using shape/rate parameterization, so its prior mean is `shape / rate`.

The package still supports lognormal priors for both Weibull parameters as an explicit sensitivity option by setting `prior_family = "lognormal"`. In that case the model uses `log(shape) ~ N(mu_log_shape, sd_log_shape)` and `log(scale) ~ N(mu_log_scale, sd_log_scale)`. The Stan code also enforces generous upper bounds on `shape` and `scale`, based on the data range, to prevent numerical overflow during sampling.

**Adaptive strategy**:
1. Work only with valid observed event times: finite, positive times where `status == 1`.
2. Compute `med_t = median(t_obs)`.
3. Compute the log-times, `log_t = log(t_obs)`, and estimate their spread. The code first uses the median absolute deviation on the log scale and rescales it by `1.4826` to approximate a standard deviation. If that is unavailable, it falls back to `sd(log_t)`. If both fail, it uses `1`.
4. Convert that log-time spread into a rough Weibull shape estimate using the Gumbel approximation:

   `shape_0 = pi / (sqrt(6) * sd(log_t))`

   The result is clamped to `0.25 <= shape_0 <= 5` to avoid extreme prior centers from noisy small samples.

5. Set an internal log-scale prior center and spread: `mu_log_shape = log(shape_0)` and `sd_log_shape = 0.7`.
6. Use the Weibull median identity, `median = scale * (log 2)^(1 / shape)`, to solve for a scale center:

   `mu_log_scale = log(med_t) - (1 / shape_0) * log(log(2))`

   This comes directly from the Weibull survival function `S(t) = exp(-(t / scale)^shape)`. Setting `S(t) = 0.5` gives `(t / scale)^shape = log(2)`, hence `t = scale * (log 2)^(1 / shape)`.

7. Set an internal `sd_log_scale` from the robust log-time spread, clamped to `0.7` to `1.5`.
8. If there are fewer than 5 observed events, widen the priors slightly: `sd_log_shape >= 0.9` and `sd_log_scale >= 1.0`.
9. Keep the shape prior on the lognormal scale. Convert only the internal scale prior into a Gamma prior with the same mean and variance:

   `gamma_shape = mean^2 / variance`

   `gamma_rate = mean / variance`

   where, for a lognormal prior with log-mean `mu` and log-SD `sigma`:

   `mean = exp(mu + sigma^2 / 2)`

   `variance = (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)`

The important point is that the data-adaptive calculation still uses the log-time scale because that is the natural scale for estimating Weibull shape and scale centers. The shape prior is passed to Stan as lognormal, while the scale prior is passed as a Gamma prior matched to the scale center and spread.

**What the Gumbel approximation means**:

This is a distributional shortcut for getting an approximate Weibull shape from the observed spread of log event times.

Start with the Weibull survival function:

`S(t) = exp(-(t / lambda)^alpha)`

If `T` follows that Weibull distribution, then:

`Y = (T / lambda)^alpha`

has a standard exponential distribution. This is useful because it lets us rewrite the event time:

`T = lambda * Y^(1 / alpha)`

Taking logs gives:

`log(T) = log(lambda) + (1 / alpha) * log(Y)`

This is where the earlier "noise" term comes from: `log(Y)` is the random part. Since `Y` is standard exponential, `log(Y)` has a known extreme-value/Gumbel-type distribution. Its standard deviation is:

`sd(log(Y)) = pi / sqrt(6)`

The scale parameter `lambda` only shifts `log(T)` by `log(lambda)`, so it affects the location of the log-times but not their spread. The shape parameter `alpha` rescales the random part, so:

`sd(log(T)) ≈ pi / (alpha * sqrt(6))`

Rearranging gives:

`alpha ≈ pi / (sqrt(6) * sd(log(T)))`

That is the formula used to compute the rough prior center `shape_0`. In practice, this is a moment-matching heuristic: if the observed event times are very spread out on the log scale, the implied Weibull shape is smaller; if they are tightly clustered, the implied shape is larger.

The same relationship appears in standard survival-regression parameterizations. For example, the Weibull accelerated-failure-time model can be written as a location-scale model for `log(T)` with an extreme-value error distribution. In R's `survival::survreg` documentation, the Weibull distribution is described as arising when the logarithm of survival time has an extreme-value distribution; the `survreg` scale parameter maps to `1 / shape`, and the linear predictor maps to `log(scale)`.

**Reviewer notes / rationale**:
- We use **event times only** (not censored) to avoid biasing the adaptive prior toward shorter times.
- The Gumbel approximation is a reasonable way to center the Weibull shape prior because it uses a known relationship between the Weibull distribution and log-time variance. It is a moment-matching center, not a precise estimator and not a claim that the prior itself is published as a complete procedure.
- The clamps on `shape_0` and prior standard deviations are important. Without them, a small or homogeneous set of observed events could create an unrealistically tight or extreme prior.
- The mixed default is appropriate because it preserves the original lognormal shape prior while matching the published Gamma prior family for the Weibull scale parameter. The adaptive Gamma scale hyperparameters are derived from the same weak data-adaptive scale center used by the lognormal option, so the practical fitted behavior should remain close while matching the prior family expected by the original methodology for scale.
- The lognormal prior option remains useful for sensitivity analysis. Development comparisons in `dev/compare_weibull_priors.R` found similar diagnostics, WAIC and imputations across the two families in ordinary scenarios, with differences mainly in fragile/heavily censored tail behavior.
- From a reviewer perspective, this is defensible as a weakly data-adaptive default, especially because the used priors are stored in `fit$priors` and users can supply explicit priors. It would be less appropriate if the package silently used very tight priors or if the priors were not inspectable.
- The main limitation is that censoring is ignored when constructing the adaptive prior. That is intentional to avoid treating censoring times as event times, but in very heavily censored datasets the observed events may not represent the full survival distribution. Sensitivity checks with custom priors are recommended in that setting.
- We enforce **generous upper bounds** on `shape` and `scale` in Stan to prevent numerical overflow during sampling. These are not intended to be informative constraints; they are safety bounds scaled to the observed time range.
- If your time unit is unusual (e.g., hours vs years), consider rescaling time or supplying custom priors for interpretability.

**Formula references**:
- The NIST/SEMATECH e-Handbook Weibull distribution page gives the Weibull CDF, survival function, inverse survival function and median. For the standard Weibull, the median is `(log 2)^(1 / shape)`; adding a scale parameter gives `scale * (log 2)^(1 / shape)`.
- The NIST/SEMATECH e-Handbook extreme-value type I page identifies the extreme-value type I distribution as the Gumbel distribution and gives its standard deviation as `beta * pi / sqrt(6)`. In the log-Weibull representation, `beta = 1 / shape`.
- The R `survival` package documentation for `survreg.distributions` describes the Weibull model as the case where `log(time)` has an extreme-value distribution, and notes the mapping between `survreg`'s scale parameter and Weibull shape.
- Standard survival-analysis texts such as Kalbfleisch and Prentice, Klein and Moeschberger, and Lawless describe the Weibull model through both hazard/survival functions and log-time accelerated-failure-time representations.

**Fallback defaults** (if adaptive fails):
- `prior_family = "gamma_scale"`
- `mu_log_shape = 0, sd_log_shape = 1`
- `scale_prior_shape = 1, scale_prior_rate = 1`

The default prior object also includes lognormal scale fields so users can switch `prior_family` to `"lognormal"` without rebuilding the list.

**Custom example**:
```r
my_priors <- get_default_priors("weibull")
my_priors$mu_log_shape <- log(1.5)
my_priors$sd_log_shape <- 0.5
my_priors$scale_prior_shape <- 3
my_priors$scale_prior_rate <- 0.02
fit <- impute(lung, time = "time", status = "status", distribution = "weibull", priors = my_priors)
```

**Lognormal sensitivity example**:
```r
my_priors <- get_default_priors("weibull")
my_priors$prior_family <- "lognormal"
my_priors$mu_log_shape <- log(1.5)
my_priors$sd_log_shape <- 0.5
my_priors$mu_log_scale <- log(200)
my_priors$sd_log_scale <- 0.8
fit <- impute(lung, time = "time", status = "status", distribution = "weibull", priors = my_priors)
```

#### Exponential

The exponential model has a single parameter - `rate`. The prior is a Gamma distribution specified by `rate_prior_shape` and `rate_prior_rate`.

**Adaptive strategy**:
1. Compute `med_t = median(t_obs)` from observed events.
2. Use the exponential median identity, `median = log(2) / rate`, to get an empirical center:

   `rate_0 = log(2) / med_t`

   This is the exponential special case of the survival identity `S(t) = exp(-rate * t)`. Setting `S(t) = 0.5` gives `rate * t = log(2)`.

3. Set a Gamma prior with shape `2` and rate `2 / rate_0`. Because the mean of `Gamma(shape, rate)` is `shape / rate`, this prior has mean `rate_0`.

**Reviewer notes / rationale**:
- The exponential prior is intentionally simple because the model has only one parameter.
- `shape = 2` gives a moderate amount of prior spread while keeping the prior proper and centered on the observed time scale.
- The median identity is a standard property of the exponential distribution; the adaptive part is using the observed event-time median as an empirical center.
- From a reviewer perspective, this is reasonable as a default centering rule. The bigger question is usually whether the exponential model itself is plausible, because it assumes a constant hazard over time.
- If the time scale is far from 1, adaptive priors or rescaling time are recommended.

**Fallback defaults**:
- `rate_prior_shape = 1, rate_prior_rate = 1` 

**Custom example**:
```r
my_priors <- list(rate_prior_shape = 3, rate_prior_rate = 0.01)
fit <- impute(lung, time = "time", status = "status", distribution = "exponential", priors = my_priors)
```

#### Lognormal

Stan fits `mu` (location) and `sigma` (scale) for the lognormal distribution. Priors expected are:

- `mu_prior_mean`, `mu_prior_sd` (normal prior on mu)
- `sigma_prior_sd` (half-normal prior on sigma)

**Adaptive strategy**:
1. Compute `log_t = log(t_obs)` from observed events.
2. Set `mu_prior_mean = median(log_t)`. This centers the prior median event time near the observed median event time.
3. Set `mu_prior_sd = 1.5`, which is deliberately broad on the log-time scale.
4. Estimate the log-time spread with `mad(log_t)` and use it as the half-normal prior scale for `sigma`.
5. Clamp `sigma_prior_sd` to `0.2` to `2.5` so that very homogeneous or very noisy observed events do not create an unusably narrow or broad prior.

**Reviewer notes / rationale**:
- Priors are placed on `mu` and `sigma` on the log‑time scale, which keeps the model stable and interpretable across time units.
- The lognormal adaptive prior is the most direct of the three because the model is already parameterized on `log(T)`.
- From a reviewer perspective, centering `mu` at the median observed log-time is reasonable and transparent. The fixed `mu_prior_sd = 1.5` is broad enough to avoid overfitting the prior center in most survival-time units.
- As with the Weibull prior, heavy censoring can make observed event times unrepresentative of the full event-time distribution, so sensitivity checks are recommended when the tail drives the analysis.
- The median/log-time relationship is a standard property of the lognormal distribution: if `log(T) ~ Normal(mu, sigma)`, then the median of `T` is `exp(mu)`.

**Fallback defaults**:
- `mu_prior_mean = 0, mu_prior_sd = 2, sigma_prior_sd = 1`.

**Custom example**:
```r
my_priors <- list(
  mu_prior_mean = log(150),
  mu_prior_sd = 0.7,
  sigma_prior_sd = 0.5
)
fit <- impute(lung, time = "time", status = "status", distribution = "lognormal", priors = my_priors)
```

In all cases, you can see what priors were used by inspecting `fit$priors` after the run completes.

### 1.2 Prior References and What They Justify

The references below justify the distributional identities and moment-matching steps used by the adaptive priors. They do not imply that the exact `bird` defaults, clamps or fallback constants are themselves a named published prior.

- Weibull survival/quantile/median formulas: NIST/SEMATECH e-Handbook of Statistical Methods, "Weibull Distribution" (https://www.itl.nist.gov/div898/handbook/eda/section3/eda3668.htm); Lawless (2003), *Statistical Models and Methods for Lifetime Data*; Klein and Moeschberger (2003), *Survival Analysis: Techniques for Censored and Truncated Data*.
- Weibull as a log-time extreme-value/Gumbel model: R `survival::survreg.distributions` documentation (https://rdrr.io/cran/survival/man/survreg.distributions.html); Kalbfleisch and Prentice (2002), *The Statistical Analysis of Failure Time Data*; Lawless (2003).
- Gumbel standard deviation and method-of-moments relationship: NIST/SEMATECH e-Handbook of Statistical Methods, "Extreme Value Type I Distribution" (https://www.itl.nist.gov/div898/handbook/eda/section3/eda366g.htm).
- Robust log-time spread via MAD: Hampel (1974) and Rousseeuw and Croux (1993) are standard references for median absolute deviation and robust scale estimation; R's `stats::mad()` uses the normal-consistency constant approximately equal to `1.4826`.
- Gamma shape/rate moment matching: this follows from the standard Gamma mean and variance formulas, `mean = shape / rate` and `variance = shape / rate^2`, giving `shape = mean^2 / variance` and `rate = mean / variance`.
- Weakly informative/data-scaled defaults: Gelman et al. (2008) discuss weakly informative priors as regularizing defaults. The `bird` adaptive priors follow that philosophy by centering broad priors on the observed event-time scale while retaining user override and sensitivity-analysis options.

Full citations:
- Gelman, A., Jakulin, A., Pittau, M. G., and Su, Y.-S. (2008). A weakly informative default prior distribution for logistic and other regression models. *The Annals of Applied Statistics*, 2(4), 1360-1383.
- Hampel, F. R. (1974). The influence curve and its role in robust estimation. *Journal of the American Statistical Association*, 69(346), 383-393.
- Kalbfleisch, J. D., and Prentice, R. L. (2002). *The Statistical Analysis of Failure Time Data*, 2nd edition. Wiley.
- Klein, J. P., and Moeschberger, M. L. (2003). *Survival Analysis: Techniques for Censored and Truncated Data*, 2nd edition. Springer.
- Lawless, J. F. (2003). *Statistical Models and Methods for Lifetime Data*, 2nd edition. Wiley.
- Rousseeuw, P. J., and Croux, C. (1993). Alternatives to the median absolute deviation. *Journal of the American Statistical Association*, 88(424), 1273-1283.

### 1.3 Nonparametric LDDP Model

The nonparametric engine uses a Dirichlet Process to let the data shape the survival function directly. Priors control how flexibly the model clusters observations and how concentrated the baseline measure is.

When you don't supply a `prior` argument, the code uses these defaults:

```r
prior <- list(
  a0 = 10, b0 = 10,       # Gamma(10, 10) prior on DP concentration alpha
  nu = 4,                 # degrees of freedom for base measure
  m0 = 0,                 # intercept mean
  S0 = matrix(25, 1, 1),  # intercept variance
  psiinv = matrix(1, 1, 1),
  tau1 = 8.01,
  taus1 = 8.01,
  taus2 = 4.01
)
```

These settings encourage moderate clustering and a diffuse baseline that lets the data dominate. We use a more conservative prior on the DP concentration (Gamma(10, 10), mean 1) to reduce extreme tail behavior under heavy right-censoring while keeping alpha learnable. The specific numeric values otherwise align with DPpackage examples and documentation.

**Empirical justification for the `Gamma(10, 10)` default**:
- This was chosen as a regularization change from an earlier `Gamma(10, 1)` setting.
- In repeated sensitivity runs (real `lung` data and simulated right-censored datasets), `Gamma(10, 10)` consistently reduced posterior `alpha` and the average number of occupied clusters (`ncluster`), typically by about 70-80%.
- In most runs this also reduced upper-tail imputation spread (e.g., 95th/99th percentiles), though not in every single random seed under very heavy censoring.
- Therefore, we treat this as a stability-oriented default rather than a guarantee against extreme draws.

**Recommended practice**:
- Keep `Gamma(10, 10)` as the default for routine use.
- For analyses sensitive to upper tails, run prior sensitivity checks (e.g., compare `Gamma(10, 10)` with `Gamma(10, 1)` or fixed `alpha`) and report robustness.

If you want tighter control over the concentration parameter you can fix `alpha` directly:

```r
custom_prior <- list(
  a0 = -1, b0 = -1,   # signals fixed alpha
  alpha = 2,
  nu = 5,
  m0 = 0,
  S0 = matrix(16, 1, 1),
  psiinv = matrix(0.5, 1, 1),
  tau1 = 6.01, taus1 = 6.01, taus2 = 3.01
)
fit_np <- impute(lung, time = "time", status = "status", model = "nonparametric", prior = custom_prior)
```

The LDDP model also accepts an `mcmc` list (`nburn`, `nsave`, `nskip`, `ndisplay`) to control the Fortran sampler. Defaults are `nburn = 1000, nsave = 4000, nskip = 1, ndisplay = 100`.


---

## 2. MCMC Defaults and Convergence Criteria

### 2.1 Parametric models (Stan/rstan)

Parametric models run 4 chains of 1000 warmup + 1000 sampling iterations each, giving 4000 total posterior draws. We chose those numbers to balance speed with reliability. Most datasets converge comfortably within that budget and the multi-chain setup helps detect any sampling pathologies.

We also set `adapt_delta = 0.95` and `max_treedepth = 12` to keep divergences low. Those settings are a bit conservative, which slows things down slightly but avoids reruns.

After each fit, we check that `max(Rhat) <= 1.1` across all parameters, this is the standard convergence threshold. We also look at effective sample sizes (ESS) and warn if bulk or tail ESS drop below 100. Divergent transitions and treedepth hits get reported but don't cause automatic failure. You can inspect traces with `plot(fit, type = "trace")` if something looks strange.

If convergence flags appear, the simplest fix is to double the warmup and sampling iterations. 

### 2.2 Nonparametric model (LDDP via Fortran)

The LDDP engine defaults to 1000 burn-in iterations followed by 4000 saved iterations, with no thinning. Progress messages print every 100 iterations so you can see it's running.

Because the Fortran code doesn't compute Rhat or ESS, we simply check that the run completes without errors and that posterior draws look reasonable. 

---

## 3. Model Fit: WAIC

WAIC (Widely Applicable Information Criterion) is a Bayesian estimate of out-of-sample predictive accuracy. It is used here as a compact model-fit metric for comparing alternative parametric survival models fitted to the same dataset, such as Weibull vs exponential vs lognormal. Lower WAIC values indicate better expected predictive performance after accounting for model flexibility.

The package computes WAIC from posterior draws of the pointwise log likelihood. For each observation and each posterior draw, the Stan models record the appropriate likelihood contribution:

- Observed events use the event-time density, `log f(t)`.
- Right-censored observations use the survival probability, `log S(t)`, because the observed information is that the event time exceeded the censoring time.

This censoring-aware calculation matters. Treating censored times as if they were observed event times would overstate early failures and bias model comparison toward distributions that put too much mass at shorter times.

For each observation `i`, the calculation is:

- `lppd_i = log(mean_s(exp(log_lik[s, i])))`, the log pointwise posterior predictive density.
- `p_waic_i = var_s(log_lik[s, i])`, the effective number of parameters contributed by that observation.
- `WAIC = -2 * (sum_i(lppd_i) - sum_i(p_waic_i))`.

The implementation uses a subtract-the-maximum calculation when exponentiating log likelihoods, which improves numerical stability for small likelihood values.

**How it is used in `bird`**:
- WAIC is computed automatically for Stan-backed parametric fits when pointwise log likelihoods are available.
- Results are stored in `fit$fit_metrics$waic`, including `waic`, `lppd` and `p_waic`.
- `print(fit)` reports WAIC for a single parametric fit.
- `compare_models()` prints WAIC alongside convergence status for each successful model, making it easier to compare parametric candidates on the same data.

**Reviewer notes / rationale**:
- WAIC is intended for comparing models fitted to the same response data. Absolute WAIC values are not meaningful on their own.
- Differences are more important than raw values, and small differences should not be over-interpreted.
- WAIC should be considered alongside convergence diagnostics, posterior predictive behavior, imputation summaries and substantive plausibility.
- WAIC is currently not computed for the nonparametric LDDP engine because that Fortran sampler does not expose the same pointwise log-likelihood draws in the fitted object.
- WAIC can become unstable when individual observations have very high posterior variance in log likelihood, especially in small datasets or heavily censored datasets. In those cases, sensitivity checks and visual diagnostics are more informative than relying on WAIC alone.

---

## 4. Grouped Analysis Design

When you pass a `groups` argument, the package splits the data by that variable and fits a separate model to each subset.

We chose this "separate models" approach because it's straightforward and maximally flexible. The downside is that you need enough data in each group, at least 2 observed events per group, or the fit will fail.

Once you have at least two successful groups, the package computes a pooled Wald test (essentially a log-rank test that pools coefficients across imputed datasets using Rubin's rules). This tells you whether survival differs significantly between groups. You also get descriptive summaries (median survival, event probabilities, survival curves) for each group.

---

## 5. Group Comparisons: Pooled Cox D1 Wald Test

For grouped analyses, `bird` reports a pooled Cox D1 Wald test. This is a global test for whether survival differs between groups after accounting for uncertainty introduced by multiple imputation.

This is based on the established D1 multivariate Wald test for multiply imputed data. The D1 approach was developed for testing several regression parameters jointly after multiple imputation, rather than testing each coefficient separately. This is useful for grouped survival comparisons because a group variable with more than two levels produces multiple Cox model contrasts that should be tested together.

The calculation works by fitting a Cox proportional hazards model separately to each completed dataset:

```r
survival::coxph(Surv(time, status) ~ group, data = completed_dataset)
```

Each fit estimates the log hazard-ratio contrasts for the group variable and their covariance matrix. The package then pools those estimates across imputations using the D1/Rubin framework:

- `bbar` is the average Cox coefficient vector across completed datasets.
- `W` is the average within-imputation covariance matrix.
- `B` is the between-imputation covariance matrix of the coefficient estimates.
- `m` is the number of usable completed datasets.
- `k` is the number of group contrasts.
- `r1 = (1 + 1 / m) * trace(B %*% solve(W)) / k` is the average relative increase in variance due to imputation.
- `T_D1 = (1 + r1) * W` is the D1 covariance approximation used for the joint test.
- The D1 test statistic is `(bbar' * solve(T_D1) * bbar) / k`.
- The p-value is computed from an F distribution with numerator degrees of freedom `k` and D1 denominator degrees of freedom.

The denominator degrees of freedom follow the large-sample D1 formula:

```r
t = k * (m - 1)

if (t > 4) {
  df2 = 4 + (t - 4) * (1 + (1 - 2 / t) / r1)^2
} else {
  df2 = t * (1 + 1 / k) * (1 + 1 / r1)^2 / 2
}
```

If `r1` is zero, there is no estimated between-imputation contribution and the denominator degrees of freedom are treated as infinite. In that limiting case the F-reference test is equivalent to the corresponding large-sample Wald test.

For two groups, there is one contrast, so the pooled coefficient can also be reported as a hazard ratio: `exp(coef)`, with an approximate 95% confidence interval from the pooled standard error. For three or more groups, the main result is the global Wald test because there are multiple contrasts relative to the reference group.

We describe this as "log-rank equivalent" because a Cox model with only a group indicator tests the same broad question as a log-rank comparison: whether the survival experience differs between groups. The Cox formulation is more convenient here because it gives coefficient estimates and covariance matrices that can be pooled cleanly across multiple imputed datasets.

**How it is used in `bird`**:
- `pool_cox_group_test()` fits and pools the Cox models across completed datasets.
- `calculate_log_rank_test()` returns the D1 Wald statistic, numerator degrees of freedom, denominator degrees of freedom and p-value for the global group comparison.
- `calculate_cox_ph_test()` returns the same p-value and, for two-group analyses, a pooled hazard ratio with confidence interval.
- The grouped print summary reports these results under the global group effect and hazard-ratio sections.

**Reviewer notes / rationale**:
- The method is an implementation of the D1 family of multivariate Wald tests for multiple-imputation inference, as described by Li, Raghunathan and Rubin and summarized in Schafer and van Buuren.
- This is more appropriate than a plain chi-square Wald test because the F reference distribution incorporates finite-imputation uncertainty through `r1` and the D1 denominator degrees of freedom.
- This is used for comparing groups within a fitted grouped analysis. It is not the criterion used to choose between imputation model families; model-family comparison uses summaries such as WAIC, convergence diagnostics and clinical/imputation behavior.
- The test uses the completed datasets generated by the imputation model, so the result reflects both the observed data and the imputed event times.
- Fits that produce singular or non-finite Cox estimates are skipped. If fewer imputations are usable, the package reports that note where the hazard-ratio helper is called.
- The usual Cox proportional hazards interpretation still applies. A significant p-value indicates evidence of a group difference in survival, but the hazard-ratio interpretation assumes proportional hazards.
- In small samples, highly separated groups or heavily censored data, the Cox fits can be unstable. In those cases, the pooled Wald result should be interpreted alongside survival curves, event probabilities and median survival summaries.

**References**:
- Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*. Wiley.
- Li, K.-H., Raghunathan, T. E., and Rubin, D. B. (1991). Large-sample significance levels from multiply imputed data using moment-based statistics and an F reference distribution. *Journal of the American Statistical Association*, 86(416), 1065-1073.
- Schafer, J. L. (1997). *Analysis of Incomplete Multivariate Data*. Chapman & Hall/CRC.
- van Buuren, S. (2018). *Flexible Imputation of Missing Data*, 2nd edition. Chapman & Hall/CRC. Section 5.3 describes D1 as the multivariate Wald test for multi-parameter inference after multiple imputation.
- The `mice::D1()` documentation implements the same D1 multivariate Wald-test family for multiply imputed model comparisons. The `mice` project also documents Cox model support for D1-style pooling using `survival::coxph()`, with complete-data degrees of freedom based on the number of events minus the number of parameters.

---

## 6. Clinical Time-Horizon Summaries

In addition to model-fit and group-comparison statistics, `bird` reports clinical-style summaries that are easier to interpret directly: median survival, event probability by fixed time horizons, survival probability differences between groups and the amount of additional time imputed for censored observations.

These summaries are descriptive rather than model-selection criteria. They are designed to answer practical questions such as "what is the probability of the event by 1 year?" or "how much longer are censored observations imputed to survive?"

### 6.1 Event Probability at Fixed Horizons

`event_probability()` estimates the cumulative probability of having experienced the event by one or more time points. This is a CDF-style quantity: `P(T <= t)`.

The defaults depend on the time unit:

- For day-scale data, printed summaries use `30, 90, 180, 365` and `event_probability()` defaults to `180, 365`.
- For month-scale data, printed summaries use `1, 3, 6, 12` and `event_probability()` defaults to `6, 12`.
- For year-scale data, printed summaries use `0.25, 0.5, 1, 2` and `event_probability()` defaults to `0.5, 1`.

For the original data, the function uses observed event times only. For the imputed result, it computes the event probability separately in each completed dataset and then summarizes across imputations with a mean, standard deviation and 2.5% / 97.5% quantiles.

By default, probabilities are estimated with a logspline CDF. If the `logspline` package is unavailable, the function falls back to the empirical CDF. Users can request the empirical method directly with `method = "ecdf"`.

**Reviewer notes / rationale**:
- Fixed horizons make results comparable across models and groups.
- The event probability is reported as a cumulative event probability, not as a survival probability. Survival probability at the same horizon is `1 - event_probability`.
- Logspline smoothing gives a less stepwise estimate than the ECDF, which is useful for small completed datasets, but it is still a descriptive summary and should be checked visually when tails are sparse.
- Horizons must be interpreted in the same time unit as the input data.

### 6.2 Median Survival

For the original right-censored data, median survival is estimated from the Kaplan-Meier curve. For completed datasets, the package computes a Kaplan-Meier median for each imputed dataset and then summarizes those medians across imputations.

The printed imputed median is the mean of the completed-dataset medians. The accompanying interval is based on the empirical 2.5% and 97.5% quantiles across imputations.

**How it is used in `bird`**:
- `calculate_clinical_metrics()` computes original and imputed median survival summaries.
- `print(fit)` reports the original Kaplan-Meier median, the imputed median and the change from original to imputed.
- Grouped and model-comparison summaries use these medians as practical effect-size summaries alongside event probabilities and plots.

**Reviewer notes / rationale**:
- This interval represents variation across imputed completed datasets. It is not the same thing as a single Kaplan-Meier confidence interval.
- If the Kaplan-Meier curve never drops below 0.5, the median may be unavailable for that dataset.
- Median survival is useful for communication, but it can hide tail behavior. It should be read alongside density, survival-curve and event-probability summaries.

### 6.3 Imputation Gain

For censored observations, `bird` keeps the original censoring time and records the completed event time after imputation. The imputation gain is the difference between the imputed event time and the original censoring time.

The package summarizes the original censored times, imputed times and imputation gains using simple descriptive statistics. It also checks whether all gains are positive, which should hold because censored observations represent `T > censoring_time`.

**How it is used in `bird`**:
- `calculate_imputation_summary()` compares original censored times, completed event times and imputation gains.
- `print(fit)` reports the number of censored observations, median original censored time, median imputed time and median gain.
- These summaries are meant as imputation diagnostics rather than inferential tests.

**Reviewer notes / rationale**:
- Imputation gain is a diagnostic for whether the imputation behaves consistently with right-censoring.
- Very large gains can be plausible under heavy censoring, but they should prompt sensitivity checks because they often reflect tail assumptions.

### 6.4 Group Survival Probability Differences

For grouped analyses, the comprehensive summary reports a survival probability difference at a one-year-equivalent horizon:

- `365` for day-scale data
- `12` for month-scale data
- `1` for year-scale data

The calculation first estimates each group's cumulative event probability at that horizon from the completed datasets. It then converts that to survival probability with `1 - event_probability`, identifies the best and worst groups, and reports the absolute difference.

**How it is used in `bird`**:
- `calculate_survival_probability_diff()` computes the best-vs-worst survival probability difference for grouped analyses.
- `print_comprehensive_group_comparison_safe()` reports this as the one-year-equivalent survival probability difference.
- `event_probability_groups()` provides the more general table when users want explicit horizons for every group.

**Reviewer notes / rationale**:
- This is an interpretable effect-size summary, not a formal hypothesis test. The formal global group test is the pooled Cox Wald test.
- The one-year-equivalent default gives a familiar benchmark across common survival-analysis time scales.
- If a study has a more relevant clinical horizon, users should call `event_probability()` or `event_probability_groups()` with explicit times.

---

## 7. Transparency 

Every completed dataset includes columns that let you trace back to what was originally there:

- `time` preserves the original observed/censoring time and `imputed_time` stores the completed event time
- `original_status` preserves the status value before imputation
- `was_censored` flags which observations were censored (so you can see exactly what got imputed)
- `dataset_id` tracks which imputation each row belongs to when you generate multiple datasets
- `.imp` provides labeling for long-format exports

If those columns get in the way of downstream analyses (e.g., you're passing data to external software that doesn't expect them), just use `export(fit, "clean_data", format = "csv", include_original = FALSE)` or drop them manually after calling `complete()`.

---

## 8. Export Design and Limitations

The `export()` function can save completed datasets in a few common formats. Right now, only CSV and RDS are fully implemented:

- **CSV**: One file per dataset (or combined file for long format). Universally readable and easy to share.
- **RDS**: R binary format that preserves all column types and attributes. Best for re-importing into R later.

---

## 9. Legacy Fortran I/O and CRAN NOTE

The nonparametric LDDP implementation vendors legacy DPpackage Fortran code that uses
Fortran file I/O (`open/read/write/rewind/close`) internally. This can trigger the
standard CRAN compiled-code NOTE about `__gfortran_st_*` symbols.

For containment, the package wrapper executes these routines inside a temporary directory
created at runtime (`tempfile("dppackage_")`), restores the original working directory,
and cleans temporary files on exit. This behavior is implemented in
`LDDPsurvival.default` (`R/LDDPsurvival.R`).

This design does not remove the NOTE itself, but it keeps I/O scoped to temporary
locations and avoids persistent file side effects in normal usage.

---

## 10. Known Limitations 

We've prioritised getting the core imputation workflow right (parametric and nonparametric, single cohort and grouped) before branching into more specialised use cases. Further development potentially includes:

**Not yet implemented**:

1. Covariate-adjusted imputation models. Adding covariates like `age` or `sex` directly into the imputation model would require a reworked Stan program and new prior structures.

2. Left truncation and interval censoring. The package assumes right-censoring exclusively.

3. Hierarchical group models. Groups are independent at the moment. A partial-pooling model could borrow strength if you have many small groups, but that would add complexity.

4. Time-varying covariates. Not supported.
