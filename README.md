# bird <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/Jamie-Wilson-UL/bird/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jamie-Wilson-UL/bird/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

`bird` provides Bayesian imputation for right-censored survival data.

It supports:
- Parametric models (`weibull`, `exponential`, `lognormal`) via Stan
- Nonparametric LDDP imputation
- Grouped analyses
- Multi-model comparisons

## Installation

```r
# install.packages("devtools")
devtools::install_github("Jamie-Wilson-UL/bird")

library(bird)
```

If you plan to use parametric engines, `rstan` is used with package-precompiled Stan models.

- CRAN binary installs (Windows/macOS) usually do not need local Stan compilation.
- Source installs still need a working C++ toolchain at install time:
  - Windows: RTools
  - macOS: `xcode-select --install`
  - Linux: `g++` and `make` (e.g. `build-essential`)

Recommended runtime settings:

```r
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## One-Function Workflow

```r
library(survival)
lung <- survival::lung

fit <- impute(
  lung,
  n_imputations = 5,
  distribution = "weibull",
  time_unit = "days"
)
```

Time units:
- `time_unit` defaults to `"days"` if omitted.
- `time_unit` controls labels/reporting and default time-horizon summaries.
- Input `time` values are assumed to already be in that unit (no automatic conversion).

`impute()` auto-detects common time/status columns. If names are unusual, specify them explicitly:

```r
fit <- impute(
  lung,
  time = "time",
  status = "status",
  n_imputations = 5,
  distribution = "lognormal"
)
```

Switching engines is just argument changes:

```r
fit_np  <- impute(lung, model = "nonparametric", n_imputations = 5)
fit_grp <- impute(lung, groups = "sex", distribution = "weibull", n_imputations = 5)
```

Direct functions are also available: `bayesian_impute()` (parametric) and `bayes_np_impute()` (nonparametric).

## Looking at Results

```r
print(fit)
```

Plot helpers:

```r
plot(fit, type = "survival")
plot(fit, type = "survival", n_curves = 20)
plot(fit, type = "completed_dataset_summary")
plot(fit, type = "completed_dataset_summary", dataset_id = 3)
plot(fit, type = "boxplots_comparison")
plot(fit, type = "boxplots_comparison", n_max = 20)
plot(fit, type = "boxplots_comparison", dataset_indices = c(1, 5, 10))
plot(fit, type = "density")
```

Diagnostics:

```r
plot(fit, type = "posterior")
plot(fit, type = "trace")
plot(fit, type = "pairs")
```

The `trace`/`pairs` plots require `bayesplot` (and `posterior`).

Grouped fits use the same plotting API:

```r
print(fit_grp)
plot(fit_grp, type = "survival")
plot(fit_grp, type = "completed_dataset_summary")
plot(fit_grp, type = "survival", combine_groups = FALSE)
plot(fit_grp, type = "boxplots_comparison")
```

For `plot(fit_grp, type = "completed_dataset_summary")`:
- 2 groups: full 4-panel view (histogram, density, survival, boxplot)
- >2 groups: simplified view (survival + boxplot)

### 3-Group Example

```r
lung3 <- survival::lung
lung3$age_group <- cut(
  lung3$age,
  breaks = c(-Inf, 60, 70, Inf),
  labels = c("<60", "60-70", ">70"),
  right = TRUE
)

fit_grp3 <- impute(
  lung3,
  groups = "age_group",
  distribution = "weibull",
  n_imputations = 5,
  time_unit = "days"
)

print(fit_grp3)
plot(fit_grp3, type = "survival")
plot(fit_grp3, type = "survival", combine_groups = FALSE)
plot(fit_grp3, type = "completed_dataset_summary")
```

## Comparing Multiple Models

```r
cmp <- compare_models(
  lung,
  models = c("weibull", "nonparametric"),
  n_imputations = 5
)

print(cmp)
plot(cmp, type = "survival")
plot(cmp, type = "completed_dataset_summary")
```

Accepted `models` values include:
- `"weibull"`
- `"exponential"`
- `"lognormal"`
- `"nonparametric"` (aliases `"np"`, `"lddp"`)

Completed datasets can be extracted combined or per model:

```r
combined <- complete(cmp, dataset = 1, models = "combined")
separate <- complete(cmp, dataset = 1, models = "separate")
weibull_only <- complete(cmp, dataset = 1, models = "weibull")
```

For `plot(cmp, type = "completed_dataset_summary")`:
- 2 models: full 4-panel view (histogram, density, survival, boxplot)
- >2 models: simplified view (survival + boxplot)

## Working with Completed Datasets

```r
complete(fit, dataset = 1)
all_sets <- complete(fit)
long_fmt <- complete(fit, format = "long")
```

Each completed dataset includes transparency columns (`original_time`, `original_status`, `was_censored`).

For grouped analyses:

```r
combined <- complete(fit_grp, dataset = 1, groups = "combined")
separate <- complete(fit_grp, dataset = 1, groups = "separate")
```

## Exporting

```r
export(fit, "weibull_results", format = "csv")
export(fit, "weibull_results", format = "rds")
```

```r
export(fit, "clean_results", format = "csv", include_original = FALSE)
export(fit_grp, "groups", format = "csv", groups = "separate")
```

Current export formats are CSV and RDS.

When exporting multiple datasets to RDS, set `combine = FALSE` to write one file per dataset.

## Advanced Usage

```r
my_priors <- get_default_priors("weibull")
my_priors$mu_log_shape <- 0.2
my_priors$sd_log_shape <- 0.4
my_priors$mu_log_scale <- 0.5
my_priors$sd_log_scale <- 0.3

custom_fit <- impute(
  lung,
  time = "time",
  status = "status",
  distribution = "weibull",
  n_imputations = 20,
  priors = my_priors,
  mcmc_options = list(
    chains = 4,
    iter_warmup = 2000,
    iter_sampling = 2000,
    adapt_delta = 0.9,
    max_treedepth = 12
  ),
  verbose = TRUE
)
```

You can do the same for other distributions (`get_default_priors("lognormal")`, etc.).
Nonparametric runs accept `mcmc = list(nburn = ..., nsave = ..., nskip = ...)` in the same way.

## Troubleshooting

- Columns not found: inspect `colnames(data)` and specify `time`/`status` explicitly.
- Stan compilation errors: confirm C++ toolchain setup.
- MCMC diagnostics: inspect `plot(fit, type = "trace")`, then increase iterations if needed.
- Performance: use fewer imputations or a faster distribution (`"exponential"`).
- Reusing draws: use `generate_complete_datasets()` instead of refitting.

Quick data sanity check:

```r
stopifnot(
  all(lung$time > 0),
  all(lung$status %in% c(0, 1)),
  !any(is.na(lung$time)),
  !any(is.na(lung$status))
)
```

## Summary

- `impute()` is the default entry point.
- Parametric/nonparametric/grouped behavior is controlled by `distribution`, `model`, and `groups`.
- Use `complete()` and `export()` after fitting.
- Plots provide quality checks and diagnostics.
- Custom priors and MCMC controls are available when needed.
