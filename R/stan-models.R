# Stan Model Compilation and Management
# This file handles pre-compilation and caching of Stan models

# Internal helper for looking up model names in generated stanmodels lists.
stan_model_lookup_candidates <- function(model_name) {
  unique(c(
    model_name,
    paste0(model_name, "_imputation"),
    paste0("bird_", model_name),
    paste0("bird_", model_name, "_imputation"),
    paste0("bayessurvival_", model_name),
    paste0("bayessurvival_", model_name, "_imputation")
  ))
}

# Returns a precompiled model from `stanmodels` (when present), otherwise NULL.
get_precompiled_stan_model <- function(model_name) {
  ns <- asNamespace("bird")
  if (!exists("stanmodels", envir = ns, inherits = FALSE)) {
    return(NULL)
  }
  models <- get("stanmodels", envir = ns, inherits = FALSE)
  if (!is.list(models)) {
    # rstantools-generated stanmodels may not always be a plain list
    # (e.g., simplified by sapply). Coerce while preserving names.
    models <- as.list(models)
  }
  if (length(models) == 0) {
    return(NULL)
  }

  hits <- intersect(stan_model_lookup_candidates(model_name), names(models))
  if (length(hits) == 0) {
    return(NULL)
  }

  model <- models[[hits[[1]]]]
  if (inherits(model, "stanmodel")) model else NULL
}

# Resolve Stan backend option with backward compatibility for old key.
get_stan_backend_option <- function() {
  backend_new <- getOption("bird.stan_backend", NULL)
  if (!is.null(backend_new)) {
    return(backend_new)
  }
  getOption("bayessurvival.stan_backend", "auto")
}

# Select backend with option bird.stan_backend = auto|precompiled|runtime
select_stan_backend <- function(model_name) {
  backend <- get_stan_backend_option()
  valid <- c("auto", "precompiled", "runtime")
  if (!backend %in% valid) {
    warning(
      "Invalid option bird.stan_backend='", backend,
      "'. Falling back to 'auto'."
    )
    backend <- "auto"
  }

  if (backend == "auto") {
    return(if (!is.null(get_precompiled_stan_model(model_name))) "precompiled" else "runtime")
  }
  backend
}

#' Compile and cache Stan models
#' @keywords internal
compile_stan_models <- function() {
  # Define available models
  models <- list(
    weibull = "weibull_imputation.stan",
    exponential = "exponential_imputation.stan",
    lognormal = "lognormal_imputation.stan"
  )
  
  # Compile each model
  for (model_name in names(models)) {
    backend <- select_stan_backend(model_name)
    if (backend == "precompiled") {
      precompiled <- get_precompiled_stan_model(model_name)
      if (!is.null(precompiled)) {
        .bird_env[[paste0(model_name, "_model")]] <- precompiled
        message("Using precompiled Stan model: ", model_name)
        next
      }
      if (identical(get_stan_backend_option(), "precompiled")) {
        stop(
          "Requested precompiled Stan backend but model not found for '",
          model_name, "'."
        )
      }
    }

    check_rstan()
    stan_file <- system.file("stan", models[[model_name]], package = "bird")
    
    if (!file.exists(stan_file)) {
      stop("Stan model file not found: ", stan_file)
    }
    
    message("Compiling Stan model: ", model_name)
    
    # Compile model
    model <- rstan::stan_model(file = stan_file, model_name = paste0("bird_", model_name))
    
    # Store in package environment
    .bird_env[[paste0(model_name, "_model")]] <- model
  }
  
  invisible(TRUE)
}

#' Get compiled Stan model
#' @param model_name Name of the model (e.g., "weibull")
#' @return Compiled Stan model object
#' @keywords internal
get_stan_model <- function(model_name) {
  model_key <- paste0(model_name, "_model")

  # Return cached model if already compiled
  if (exists(model_key, envir = .bird_env, inherits = FALSE)) {
    return(.bird_env[[model_key]])
  }

  backend <- select_stan_backend(model_name)
  if (backend == "precompiled") {
    model <- get_precompiled_stan_model(model_name)
    if (is.null(model)) {
      stop(
        "Precompiled Stan model not found for '", model_name, "'. ",
        "If you are developing, set options(bird.stan_backend = 'runtime') ",
        "or generate package stanmodels via rstantools."
      )
    }
    .bird_env[[model_key]] <- model
    return(model)
  }

  # Compile on demand 
  check_rstan()
  stan_filename <- paste0(model_name, "_imputation.stan")
  stan_file <- system.file("stan", stan_filename, package = "bird")
  if (!nzchar(stan_file) || !file.exists(stan_file)) {
    stop("Stan model file not found: ", stan_filename)
  }

  message("Compiling Stan model: ", model_name)
  model <- rstan::stan_model(file = stan_file, model_name = paste0("bird_", model_name))
  .bird_env[[model_key]] <- model
  model
}

#' Check if Stan models are compiled
#' @return Logical indicating if models are ready
#' @keywords internal
stan_models_ready <- function() {
  required_models <- c("weibull_model", "exponential_model", "lognormal_model")
  all(sapply(required_models, function(x) exists(x, envir = .bird_env)))
}

#' Get default priors for a given distribution
#' @param distribution Distribution name ("weibull", "exponential", "lognormal")
#' @return List of default prior parameters
#' @export
get_default_priors <- function(distribution = "weibull") {
  switch(distribution,
    "weibull" = list(
      # Log-parameterisation defaults 
      mu_log_shape = 0,           # centered at shape ~ 1
      sd_log_shape = 1.0,         # weakly-informative on log-shape
      mu_log_scale = 0,           # unit-agnostic location
      sd_log_scale = 2.0          # weakly-informative on log-scale
    ),
    "exponential" = list(
      # Rate parameter (gamma) weak prior centered near 1
      rate_prior_shape = 1.0,
      rate_prior_rate = 1.0
    ),
    "lognormal" = list(
      # Location parameter (normal)
      mu_prior_mean = 0,
      mu_prior_sd = 2,           # slightly wider default
      # Scale parameter (half-normal)
      sigma_prior_sd = 1
    ),
    stop("Unknown distribution: ", distribution, ". Supported distributions: weibull, exponential, lognormal")
  )
}

#' Compute data-adaptive priors from observed events
#'
#' Adaptive Weibull priors use the dispersion of log-times to set a rough center
#' for log-shape via a Gumbel approximation (sd(log T) ~= pi/(alpha*sqrt(6))), and
#' center log-scale to match the observed median using the relation
#' median = scale * (log 2)^(1/shape).
#'
#' @param distribution Distribution name
#' @param t_obs Numeric vector of observed event times (>0)
#' @return List of prior hyperparameters matching prepare_stan_data expectations
#' @export
get_adaptive_priors <- function(distribution = "weibull", t_obs) {
  t_obs <- t_obs[is.finite(t_obs) & t_obs > 0]
  if (length(t_obs) < 2) {
    # Fallback to defaults if insufficient data
    return(get_default_priors(distribution))
  }
  med_t <- stats::median(t_obs)
  if (!is.finite(med_t) || med_t <= 0) {
    return(get_default_priors(distribution))
  }
  
  if (distribution == "weibull") {
    # Use dispersion of log-times to set a rough alpha0 via Gumbel approximation
    log_t <- log(t_obs)
    mad_log <- stats::mad(log_t, constant = 1)
    sd_log <- if (is.finite(mad_log) && mad_log > 0) 1.4826 * mad_log else stats::sd(log_t)
    if (!is.finite(sd_log) || sd_log <= 0) sd_log <- 1
    alpha0 <- pi / (sqrt(6) * sd_log)
    alpha0 <- max(0.25, min(5, alpha0))

    mu_log_shape <- log(alpha0)
    sd_log_shape <- 0.7

    mu_log_scale <- log(med_t) - (1/alpha0) * log(log(2))
    sd_log_scale <- max(0.7, min(1.5, 1.4826 * mad_log))

    # Guards for tiny/homogeneous data
    if (length(t_obs) < 5) {
      sd_log_shape <- max(sd_log_shape, 0.9)
      sd_log_scale <- max(sd_log_scale, 1.0)
    }

    list(
      mu_log_shape = mu_log_shape,
      sd_log_shape = sd_log_shape,
      mu_log_scale = mu_log_scale,
      sd_log_scale = sd_log_scale
    )
  } else if (distribution == "exponential") {
    # For Exponential, median = log(2)/lambda => lambda ≈ log(2)/median
    emp_rate <- log(2) / med_t
    shape <- 2.0
    rate <- shape / emp_rate
    list(
      rate_prior_shape = shape,
      rate_prior_rate = rate
    )
  } else if (distribution == "lognormal") {
    log_t <- log(t_obs)
    mu_mean <- stats::median(log_t)
    # Use robust spread for sigma prior scale
    mad_log <- stats::mad(log_t, constant = 1)
    if (!is.finite(mad_log) || mad_log <= 0) mad_log <- 1
    list(
      mu_prior_mean = mu_mean,
      mu_prior_sd = 1.5,
      sigma_prior_sd = max(0.2, min(2.5, mad_log))
    )
  } else {
    get_default_priors(distribution)
  }
} 
