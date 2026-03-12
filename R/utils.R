# Utility functions for bird package

#' Check if required packages are available
#' @keywords internal
check_required_packages <- function() {
  required <- c("rstan", "posterior", "survival", "ggplot2")
  missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
  
  if (length(missing) > 0) {
    stop("Required packages not available: ", paste(missing, collapse = ", "),
         "\nPlease install with: install.packages(c(", 
         paste0("'", missing, "'", collapse = ", "), "))")
  }
}

#' Check if rstan is installed
#' @keywords internal
check_rstan <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop(
      "The 'rstan' package is required but not installed. ",
      "Install with: install.packages('rstan')",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Format time for display
#' @param time_seconds Time in seconds
#' @return Formatted time string
#' @keywords internal
format_time <- function(time_seconds) {
  if (time_seconds < 60) {
    return(sprintf("%.1f seconds", time_seconds))
  } else if (time_seconds < 3600) {
    return(sprintf("%.1f minutes", time_seconds / 60))
  } else {
    return(sprintf("%.1f hours", time_seconds / 3600))
  }
}

#' Standardize completed dataset column order
#'
#' Reorders completed datasets so original/new censor indicators and times are
#' grouped together for easier inspection:
#' `original_status`, `status`, `was_censored`, `original_time`, `time`.
#'
#' Metadata columns (when present) are kept at the front:
#' `.imp`, `.id`, `.model`, `dataset_id`.
#'
#' @param data Completed dataset (data.frame)
#' @param time_col Name of imputed/current time column
#' @param status_col Name of imputed/current status column
#' @return Reordered data.frame
#' @keywords internal
standardize_complete_column_order <- function(data, time_col = "time", status_col = "status") {
  if (!is.data.frame(data) || ncol(data) == 0) {
    return(data)
  }

  cols <- names(data)

  meta_set <- c(".imp", ".id", ".model", "dataset_id")
  compare_block <- c("original_status", status_col, "was_censored", "original_time", time_col)

  # Always place metadata at the front in a consistent order.
  ordered_front <- intersect(c(".imp", ".id", ".model", "dataset_id"), cols)
  ordered_compare <- intersect(compare_block, cols)

  ordered <- unique(c(ordered_front, ordered_compare))
  rest <- setdiff(cols, ordered)

  data[, c(ordered, rest), drop = FALSE]
}

#' Safe evaluation with error handling
#' @param expr Expression to evaluate
#' @param error_msg Custom error message
#' @return Result of expression or error
#' @keywords internal
safe_eval <- function(expr, error_msg = "An error occurred") {
  tryCatch(
    expr,
    error = function(e) {
      stop(error_msg, ": ", e$message, call. = FALSE)
    }
  )
}


# Get data summary
get_data_summary <- function(data, time_col, status_col, verbose = TRUE) {
  
  n_total <- nrow(data)
  n_observed <- sum(data[[status_col]] == 1)
  n_censored <- sum(data[[status_col]] == 0)
  censoring_rate <- n_censored / n_total
  
  if (verbose) {
    cat("  Total observations:", n_total, "\n")
    cat("  Observed events:", n_observed, "\n") 
    cat("  Censored observations:", n_censored, "\n")
    cat("  Censoring rate:", round(censoring_rate * 100, 1), "%\n")
  }
  
  return(list(
    n_total = n_total,
    n_observed = n_observed,
    n_censored = n_censored,
    censoring_rate = censoring_rate
  ))
}


# Get default MCMC options
get_default_mcmc_options <- function() {
  return(list(
    iter_warmup = 1000,
    iter_sampling = 1000,
    chains = 4,
    adapt_delta = 0.95,
    max_treedepth = 12
  ))
}

normalize_mcmc_options <- function(mcmc_options) {
  defaults <- get_default_mcmc_options()

  if (is.null(mcmc_options)) {
    return(defaults)
  }

  if (!is.list(mcmc_options)) {
    stop("'mcmc_options' must be a list.", call. = FALSE)
  }

  mcmc_options <- mcmc_options[!vapply(mcmc_options, is.null, logical(1))]

  merged <- utils::modifyList(defaults, mcmc_options)

  merged$iter_warmup <- as.integer(merged$iter_warmup)
  merged$iter_sampling <- as.integer(merged$iter_sampling)
  merged$chains <- as.integer(merged$chains)
  merged$adapt_delta <- as.numeric(merged$adapt_delta)
  merged$max_treedepth <- as.integer(merged$max_treedepth)

  merged
}

build_rstan_control <- function(mcmc_options) {
  mcmc_options <- normalize_mcmc_options(mcmc_options)

  control <- list()
  if (is.finite(mcmc_options$adapt_delta)) {
    control$adapt_delta <- mcmc_options$adapt_delta
  }
  if (is.finite(mcmc_options$max_treedepth)) {
    control$max_treedepth <- mcmc_options$max_treedepth
  }
  control
}


# Extract MCMC summary (distribution-agnostic)
extract_mcmc_summary <- function(stan_fit, distribution = "weibull") {
  
  # Get parameter names based on distribution
  param_names <- get_distribution_params(distribution)
  
  # Add log_lik_total 
  all_vars <- c(param_names, "log_lik_total")
  
  draws <- posterior::subset_draws(
    posterior::as_draws_df(stan_fit),
    variable = all_vars
  )
  param_summary <- posterior::summarise_draws(draws)

  sampler_params <- tryCatch(
    rstan::get_sampler_params(stan_fit, inc_warmup = FALSE),
    error = function(e) NULL
  )
  diagnostics <- NULL
  if (!is.null(sampler_params) && length(sampler_params) > 0) {
    diagnostics <- data.frame(
      chain = seq_along(sampler_params),
      num_divergent = vapply(sampler_params, function(m) {
        if (!("divergent__" %in% colnames(m))) return(NA_real_)
        sum(m[, "divergent__"] > 0)
      }, numeric(1)),
      stringsAsFactors = FALSE
    )
  }
  
  return(list(
    parameters = param_summary,
    diagnostics = diagnostics
  ))
}

#' Get parameter names for a distribution
#' @param distribution Distribution name ("weibull", "exponential", "lognormal")
#' @return Character vector of parameter names
#' @keywords internal
get_distribution_params <- function(distribution) {
  switch(distribution,
    "weibull" = c("shape", "scale"),
    "exponential" = "rate",
    "lognormal" = c("mu", "sigma"),
    stop("Unknown distribution: ", distribution, ". Supported distributions: weibull, exponential, lognormal")
  )
}

#' Get parameter labels for plotting
#' @param distribution Distribution name ("weibull", "exponential", "lognormal")
#' @return Character vector of parameter labels for plots
#' @keywords internal
get_distribution_param_labels <- function(distribution) {
  switch(distribution,
    "weibull" = c("Shape Parameter", "Scale Parameter"),
    "exponential" = "Rate Parameter",
    "lognormal" = c("Location Parameter (mu)", "Scale Parameter (sigma)"),
    stop("Unknown distribution: ", distribution, ". Supported distributions: weibull, exponential, lognormal")
  )
} 
