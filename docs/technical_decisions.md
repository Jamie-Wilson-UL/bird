# Technical Decisions/Design Rationale

This document describes key design choices, default settings and assumptions in `bird`.

---

## 1. Prior Specification

When you run a parametric model without supplying priors, the package first tries to learn reasonable hyperparameters from the observed event times. If that's not possible (too few events or degenerate data), it falls back to generic weakly informative settings. Nonparametric models use fixed defaults inherited from DPpackage.

### 1.1 Parametric Models: Default Data-Adaptive Priors 

The parametric models (`weibull`, `exponential`, `lognormal`) aim to let the data speak wherever possible. If you don't supply a `priors` argument:

1. The code pulls out the observed event times (the uncensored observations).
2. It calls `get_adaptive_priors(distribution, t_obs)`, which inspects the median, spread and scale of those times.
3. The function sets hyperparameters that roughly center the prior on plausible parameter ranges for the dataset.
4. If there are issues (e.g., fewer than 2 events or the data are too homogeneous to compute a stable estimate), it falls back to `get_default_priors(distribution)` which gives generic weakly informative settings.

#### Weibull

Stan fits `shape` and `scale` on the log scale to keep the sampler stable. The priors you (or the package) supply are normal distributions on the logs `log(shape) ~ N(mu_log_shape, sd_log_shape)` and `log(scale) ~ N(mu_log_scale, sd_log_scale)`. That means the actual `shape` and `scale` parameters are implicitly log-normal on the positive scale, which keeps them strictly positive without hard bounds. The Stan code also enforces soft upper bounds (generous limits based on the data range) to prevent the sampler from wandering into numerically unstable territory. 

**Adaptive strategy**:
- Estimate a rough shape parameter using the dispersion of log-times (Gumbel approximation: `sd(log T) ≈ π / (shape × sqrt(6))`).
- Center the log-scale prior near the observed median, accounting for the relationship `median = scale × (log 2)^(1/shape)`.
- Use moderate spread to allow the data to override the center if needed.

**Reviewer notes / rationale**:
- We use **event times only** (not censored) to avoid biasing the adaptive prior toward shorter times.
- We enforce **generous upper bounds** on `shape` and `scale` in Stan to prevent numerical overflow during sampling. These are not intended to be informative constraints; they are safety bounds scaled to the observed time range.
- If your time unit is unusual (e.g., hours vs years), consider rescaling time or supplying custom priors for interpretability.

**Fallback defaults** (if adaptive fails):
- `mu_log_shape = 0, sd_log_shape = 1` (shape near 1)
- `mu_log_scale = 0, sd_log_scale = 2` (scale agnostic to units)

**Custom example**:
```r
my_priors <- list(
  mu_log_shape = log(1.5),
  sd_log_shape = 0.5,
  mu_log_scale = log(200),
  sd_log_scale = 0.8
)
fit <- impute(lung, distribution = "weibull", priors = my_priors)
```

#### Exponential

The exponential model has a single parameter - `rate`. The prior is a Gamma distribution specified by `rate_prior_shape` and `rate_prior_rate`.

**Adaptive strategy**:
- Use the observed median to estimate `rate ≈ log(2) / median`.
- Set a Gamma prior centered near that value: `shape = 2, rate = 2 / (log(2) / median)`.

**Reviewer notes / rationale**:
- The exponential prior is intentionally weak; if the time scale is far from 1, adaptive priors or rescaling time is recommended.

**Fallback defaults**:
- `rate_prior_shape = 1, rate_prior_rate = 1` 

**Custom example**:
```r
my_priors <- list(rate_prior_shape = 3, rate_prior_rate = 0.01)
fit <- impute(lung, distribution = "exponential", priors = my_priors)
```

#### Lognormal

Stan fits `mu` (location) and `sigma` (scale) for the lognormal distribution. Priors expected are:

- `mu_prior_mean`, `mu_prior_sd` (normal prior on mu)
- `sigma_prior_sd` (half-normal prior on sigma)

**Adaptive strategy**:
- Set `mu_prior_mean = median(log(t_obs))` to center on the observed scale.
- Use a moderate standard deviation (`mu_prior_sd = 1.5`).
- Let the prior on sigma scale with the median absolute deviation of log-times, clamped to a sensible range (0.2 to 2.5).

**Reviewer notes / rationale**:
- Priors are placed on `mu` and `sigma` on the log‑time scale, which keeps the model stable and interpretable across time units.

**Fallback defaults**:
- `mu_prior_mean = 0, mu_prior_sd = 2, sigma_prior_sd = 1`.

**Custom example**:
```r
my_priors <- list(
  mu_prior_mean = log(150),
  mu_prior_sd = 0.7,
  sigma_prior_sd = 0.5
)
fit <- impute(lung, distribution = "lognormal", priors = my_priors)
```

In all cases, you can see what priors were used by inspecting `fit$priors` after the run completes.

### 1.2 Nonparametric LDDP Model

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
fit_np <- impute(lung, model = "nonparametric", prior = custom_prior)
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

## 3. Grouped Analysis Design

When you pass a `groups` argument, the package splits the data by that variable and fits a separate model to each subset.

We chose this "separate models" approach because it's straightforward and maximally flexible. The downside is that you need enough data in each group, at least 2 observed events per group, or the fit will fail.

Once you have at least two successful groups, the package computes a pooled Wald test (essentially a log-rank test that pools coefficients across imputed datasets using Rubin's rules). This tells you whether survival differs significantly between groups. You also get descriptive summaries (median survival, event probabilities, survival curves) for each group.

---

## 4. Transparency 

Every completed dataset includes columns that let you trace back to what was originally there:

- `time` preserves the original observed/censoring time and `imputed_time` stores the completed event time
- `original_status` preserves the status value before imputation
- `was_censored` flags which observations were censored (so you can see exactly what got imputed)
- `dataset_id` tracks which imputation each row belongs to when you generate multiple datasets
- `.imp` provides labeling for long-format exports

If those columns get in the way of downstream analyses (e.g., you're passing data to external software that doesn't expect them), just use `export(fit, "clean_data", format = "csv", include_original = FALSE)` or drop them manually after calling `complete()`.

---

## 5. Export Design and Limitations

The `export()` function can save completed datasets in a few common formats. Right now, only CSV and RDS are fully implemented:

- **CSV**: One file per dataset (or combined file for long format). Universally readable and easy to share.
- **RDS**: R binary format that preserves all column types and attributes. Best for re-importing into R later.

---

## 6. Legacy Fortran I/O and CRAN NOTE

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

## 7. Known Limitations 

We've prioritised getting the core imputation workflow right (parametric and nonparametric, single cohort and grouped) before branching into more specialised use cases. Further development potentially includes:

**Not yet implemented**:

1. Covariate-adjusted imputation models. Adding covariates like `age` or `sex` directly into the imputation model would require a reworked Stan program and new prior structures.

2. Left truncation and interval censoring. The package assumes right-censoring exclusively.

3. Hierarchical group models. Groups are independent at the moment. A partial-pooling model could borrow strength if you have many small groups, but that would add complexity.

4. Time-varying covariates. Not supported.
