# bird User Guide

---

## Getting Started

Install from GitHub and load the library:

```r
install.packages("devtools")
devtools::install_github("Jamie-Wilson-UL/bird")

library(bird)
```

If you plan to use the parametric engines (Weibull, exponential, lognormal), `rstan` is used with precompiled package models.

- **CRAN binary installs (Windows/macOS)** usually do not need local Stan compilation.
- **Source installs** still need a working C++ toolchain at install time:
  - **Windows**: Install RTools (includes C++ compiler)
  - **macOS**: Install Xcode Command Line Tools: `xcode-select --install`
  - **Linux**: Install `g++` and `make` (e.g. `sudo apt-get install build-essential`)

Recommended runtime settings:

```r
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

The nonparametric LDDP engine doesn't rely on Stan, so it will work even if you skip the toolchain setup.

Throughout the examples below, we use the classic `survival::lung` dataset:

```r
library(survival)
lung <- survival::lung
```

---

## One-Function Workflow

For most situations, you can just use `impute()`:

```r
fit <- impute(
  lung,
  n_imputations = 5,
  distribution = "weibull",    # choose "weibull", "exponential", or "lognormal"
  time_unit = "days"
)
```

Time units:
- `time_unit` defaults to `"days"` if omitted.
- `time_unit` controls labels/reporting and default time-horizon summaries.
- The input `time` values are assumed to already be in that unit (no automatic conversion is performed).

By default, `impute()` tries to detect the time and status columns (looking for common names like `time`, `status`, `duration`, etc.). If your column names are unusual it is advisable to specify them explicitly:

```r
fit <- impute(
  lung,
  time = "time",           # you can also pass indices, e.g. time = 3
  status = "status",
  n_imputations = 5,
  distribution = "lognormal"
)
```

Switching engines (distributions/groups) is just a change of arguments:

```r
fit_np  <- impute(lung, model = "nonparametric", n_imputations = 5)
fit_grp <- impute(lung, groups = "sex", distribution = "weibull", n_imputations = 5)
```

If you need direct access to the underlying functions, they are available as `bayesian_impute()` (parametric) and `bayes_np_impute()` (nonparametric). They behave exactly the same as `impute()` once you pass columns or `prepare_survival_data(...)` output.

---

## Looking at the Results

The print method gives a basic statistical summary with sensible next steps:

```r
print(fit)
```

You can check your imputations with the plotting helpers:

```r
plot(fit, type = "survival")                  # observed vs imputed survival curves
plot(fit, type = "survival", n_curves = 20)  # show more imputed survival curves
plot(fit, type = "completed_dataset_summary") # overview of one random completed dataset
plot(fit, type = "completed_dataset_summary", dataset_id = 3)  # pick a specific dataset
plot(fit, type = "boxplots_comparison")       # compares up to 10 random datasets by default
plot(fit, type = "boxplots_comparison", n_max = 20)            # increase the displayed limit
plot(fit, type = "boxplots_comparison", dataset_indices = c(1, 5, 10)) # choose exact datasets
plot(fit, type = "density")                   # mean imputed density + pointwise 95% band
```

When you want to dive deeper into model diagnostics, these are available too:

```r
plot(fit, type = "posterior")  # (random) posterior distributions
plot(fit, type = "trace")      # MCMC convergence (requires bayesplot package)
plot(fit, type = "pairs")      # parameter relationships (requires bayesplot package)
```

Note: The `trace` and `pairs` plot types rely on the `bayesplot` package. If you see a message asking you to install it, run `install.packages(c("bayesplot", "posterior"))` and the plots will work.

All of the group-specific objects (`impute(..., groups = ...)`) respond to the same plotting API. For example:

```r
print(fit_grp)
plot(fit_grp, type = "survival")
plot(fit_grp, type = "completed_dataset_summary") 
plot(fit_grp, type = "survival", combine_groups = FALSE)  # separate panel per group
plot(fit_grp, type = "boxplots_comparison")
```

For `plot(fit_grp, type = "completed_dataset_summary")`:
- 2 groups: full 4-panel view (histogram, density, survival, boxplot),
- >2 groups: simplified view (survival + boxplot).

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
plot(fit_grp3, type = "completed_dataset_summary")  # simplified (survival + boxplot) for >2 groups
```

---

## Comparing Multiple Models

If you want to compare parametric and nonparametric imputations directly, use `compare_models()`:

```r
cmp <- compare_models(
  lung,
  models = c("weibull", "nonparametric"),
  n_imputations = 5
)
```

You can then use a comparison-oriented print/plot workflow:

```r
print(cmp)
plot(cmp, type = "survival")
plot(cmp, type = "completed_dataset_summary")
```

`models` can include any combination of:
- `"weibull"`
- `"exponential"`
- `"lognormal"`
- `"nonparametric"` (also accepts aliases `"np"`, `"lddp"`)

Completed datasets are available in combined or model-specific form:

```r
# Stack all selected models into one dataset (adds .model)
combined <- complete(cmp, dataset = 1, models = "combined")

# Return one object per model
separate <- complete(cmp, dataset = 1, models = "separate")

# Extract only one model from the comparison
weibull_only <- complete(cmp, dataset = 1, models = "weibull")
```

For `plot(cmp, type = "completed_dataset_summary")`:
- 2 models: full 4-panel view (histogram, density, survival, boxplot),
- >2 models: simplified view (survival + boxplot) to avoid overcrowded plots.

---

## Working with Completed Datasets

The package stores completed datasets in a way that is deliberately transparent. You can pull out one dataset or many, in either wide, long or list form:

```r
complete(fit, dataset = 1)            # a single completed dataset
all_sets <- complete(fit)             # list of datasets
long_fmt <- complete(fit, format = "long")
```

Each dataset carries the original time/status columns alongside the imputed versions (`original_time`, `original_status`, `was_censored`), making it easy to see exactly what changed.

For grouped analyses you can choose to combine groups or keep them separate:

```r
combined <- complete(fit_grp, dataset = 1, groups = "combined")
separate <- complete(fit_grp, dataset = 1, groups = "separate")
```

---

## Exporting

When you are satisfied with a fit, export the data in the format that suits your downstream work:

```r
export(fit, "weibull_results", format = "csv")
export(fit, "weibull_results", format = "rds")
```

There are helper arguments if you want to drop the transparency columns or focus on particular datasets:

```r
export(fit, "clean_results", format = "csv", include_original = FALSE)
export(fit_grp, "groups", format = "csv", groups = "separate")
```

At the moment only CSV and RDS exports are enabled in the package. (Support for Excel/Stata/SPSS is planned, but not yet implemented.)

When exporting multiple datasets to RDS, set `combine = FALSE` to write one file per dataset.

---

## Advanced Usage

The defaults should be sufficient for most applications, but it is also possible to specify custom priors and MCMC settings. A good pattern is to start from the defaults and tweak the pieces you care about:

```r
# Start from the default Weibull priors and modify them
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

You can do the same for other distributions (`get_default_priors("lognormal")`, etc.). Nonparametric runs accept `mcmc = list(nburn = ..., nsave = ..., nskip = ...)` in the same way.

---

## Troubleshooting

A few quick checks:

- **Columns not found**: inspect with `colnames(lung)` and specify `time = "..."`, `status = "..."` explicitly (indices work too).
- **Stan compilation errors**: confirm your C++ toolchain is installed (RTools / Xcode CLT / build-essential).
- **MCMC diagnostics**: call `plot(fit, type = "trace")` to inspect convergence; increase iterations if needed.
- **Performance**: use fewer imputations (e.g., 5 instead of 50) or choose a faster distribution (`"exponential"`).
- **Reusing draws**: if you want 100 datasets but already ran 20, use `generate_complete_datasets()` on the stored posterior draws rather than refitting.

A quick sanity check to confirm your data meet the minimum requirements:

```r
stopifnot(
  all(lung$time > 0),
  all(lung$status %in% c(0, 1)),
  !any(is.na(lung$time)),
  !any(is.na(lung$status))
)
```

---

## Summary

- `impute()` is the default entry point.
- Parametric vs nonparametric vs grouped behaviour is controlled with `distribution`, `model`, and `groups` arguments.
- `complete()` and `export()` once you have a fit.
- Plots provide quality checks and diagnostics.
- Custom priors/MCMC settings are available if desired.
