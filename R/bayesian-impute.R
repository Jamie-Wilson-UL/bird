# Core Bayesian Imputation Function

#' Bayesian Imputation for Right-Censored Survival Data
#'
#' This function implements the Bayesian imputation methodology from Moghaddam et al. (2022)
#' for handling right-censored survival data. It treats censored observations as missing
#' data and uses Bayesian methods to generate posterior distributions for each censored 
#' observation.
#'
#' @param data A data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column (1 = event, 0 = censored)
#' @param groups Name of the group variable for group comparison (optional)
#' @param n_imputations Number of imputed complete datasets to generate (default: 10)
#' @param distribution Distribution to use ("weibull", "exponential", "lognormal")
#' @param priors List of prior parameters for the distribution parameters
#' @param mcmc_options List of MCMC sampling options (chains, iterations, etc.)
#' @param verbose Print progress messages (default: TRUE)
#' @param time_unit Optional label for the time unit used in outputs (default: "days")
#'
#' @return A bayesian_imputation object containing:
#' \describe{
#'   \item{original_data}{The original input data}
#'   \item{imputed_datasets}{List of n_imputations complete datasets}
#'   \item{posterior_imputations}{Matrix of posterior draws for each censored observation}
#'   \item{posterior_samples}{MCMC samples for distribution parameters}
#'   \item{model_info}{Information about the fitted model}
#'   \item{diagnostics}{MCMC convergence diagnostics}
#'   \item{mcmc_summary}{Summary statistics for all parameters}
#' }
#'
#' @details
#' Completed datasets preserve the original observations alongside their imputations. Each
#' completed dataset retains all original covariates and adds transparency columns, allowing
#' users to distinguish between original and imputed values when analysing results.
#'
#' The methodology works as follows:
#' \enumerate{
#'   \item Fit a Bayesian survival model using MCMC
#'   \item Generate one imputation per censored observation per MCMC iteration
#'   \item Store the full posterior distribution for each censored observation
#'   \item Generate complete datasets by independently sampling from each censored observation's posterior
#' }
#'
#' When groups is specified, the function splits the data by group and fits separate
#' models to each group.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(lung, package = "survival")
#' lung_clean <- na.omit(lung[, c("time", "status")])
#' 
#' # Basic imputation (single group)
#' result <- bayesian_impute(
#'   data = lung_clean,
#'   time_col = "time",  
#'   status_col = "status",  
#'   n_imputations = 20
#' )
#' 
#' # Group comparison
#' result_groups <- bayesian_impute(
#'   data = lung_clean,
#'   time_col = "time",
#'   status_col = "status",
#'   groups = "sex",
#'   n_imputations = 20
#' )
#' 
#' # View results
#' print(result)
#' plot(result, type = "survival")
#' 
#' # Generate additional datasets 
#' more_datasets <- generate_complete_datasets(
#'   result$original_data, 
#'   result$time_col, 
#'   result$status_col,
#'   result$posterior_imputations, 
#'   50  # Generate 50 more datasets
#' )
#' }
#' @export
bayesian_impute <- function(data, time_col = NULL, status_col = NULL, 
                           groups = NULL,  # Optional groups parameter
                           n_imputations = 10,
                           distribution = "weibull",
                           priors = NULL,
                           mcmc_options = NULL,
                           verbose = TRUE,
                           time_unit = NULL) {
  
  # Handle survival_data objects (auto-detection)
  if (inherits(data, "survival_data")) {
    if (is.null(time_col)) {
      time_col <- attr(data, "time_col")
    }
    if (is.null(status_col)) {
      status_col <- attr(data, "status_col")
    }
  }
  
  # If columns are not provided or not present, auto-prepare to detect them
  if (is.null(time_col) || is.null(status_col) ||
      !(time_col %in% names(data)) || !(status_col %in% names(data))) {
    if (verbose) message("Auto-preparing survival data (detecting time/status columns)...")
    prepared <- prepare_survival_data(data, verbose = FALSE, time_unit = time_unit %||% attr(data, "time_unit") %||% "days")
    data <- prepared
    time_col <- attr(prepared, "time_col")
    status_col <- attr(prepared, "status_col")
    if (verbose && !is.null(time_col) && !is.null(status_col)) {
      cat(sprintf("Detected columns: time='%s', status='%s'\n", time_col, status_col))
    }
  }
  
  # Check that time_col and status_col are provided
  if (is.null(time_col) || is.null(status_col)) {
    stop("You must specify 'time_col' and 'status_col'. ",
         "Alternatively, use prepare_survival_data() first for auto-detection.")
  }
  
  if (is.null(groups)) {
    return(bayesian_impute_single(data, time_col, status_col, 
                                 n_imputations, distribution, 
                                 priors, mcmc_options, verbose, time_unit = time_unit))
  } else {
    return(bayesian_impute_groups_safe(data, time_col, status_col, groups,
                                      n_imputations, distribution, 
                                      priors, mcmc_options, verbose, time_unit = time_unit))
  }
}

#' Bayesian Imputation for Single Group (Internal)
#'
#' This is the bayesian_impute logic, separated for single-group analysis.
#' This function is called internally by bayesian_impute() when no groups are specified.
#'
#' @param data A data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column (1 = event, 0 = censored)
#' @param n_imputations Number of imputed complete datasets to generate (default: 10)
#' @param distribution Distribution to use ("weibull", "exponential", "lognormal")
#' @param priors List of prior parameters for the distribution parameters
#' @param mcmc_options List of MCMC sampling options (chains, iterations, etc.)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A bayesian_imputation object
#' @keywords internal
bayesian_impute_single <- function(data, time_col = NULL, status_col = NULL, 
                                  n_imputations = 10,
                                  distribution = "weibull",
                                  priors = NULL,
                                  mcmc_options = NULL,
                                  verbose = TRUE,
                                  time_unit = NULL) {
  
  if (verbose) {
    cat("Starting Bayesian imputation for censored survival data...\n")
    cat("Distribution:", distribution, "\n")
    cat("Number of imputations:", n_imputations, "\n")
  }
  
  # Validate distribution
  valid_distributions <- c("weibull", "exponential", "lognormal")
  if (!distribution %in% valid_distributions) {
    stop("Invalid distribution: ", distribution, ". Supported distributions: ", 
         paste(valid_distributions, collapse = ", "))
  }
  
  # Validate inputs (skip if already validated)
  if (!inherits(data, "survival_data")) {
    if (verbose) cat("Validating input data...\n")
    validation_result <- validate_survival_data(data, time_col, status_col)
    data <- validation_result$data  # Use the potentially converted data
  } else {
    if (verbose) cat("Data already validated, skipping validation...\n")
  }
  
  # Validate sample size requirements (consistent with group validation)
  n_total <- nrow(data)
  n_observed <- sum(data[[status_col]] == 1, na.rm = TRUE)
  n_censored <- sum(data[[status_col]] == 0, na.rm = TRUE)
  
  if (n_total < 5) {
    stop("Dataset has only ", n_total, " observations (need >= 5)")
  }
  if (n_observed < 2) {
    stop("Dataset has only ", n_observed, " events (need >= 2)")
  }
  if (n_censored == 0) {
    warning("Dataset has no censored observations")
  }
  
  # Prepare data
  if (verbose) cat("Data summary:\n")
  data_summary <- get_data_summary(data, time_col, status_col, verbose)
  
  # Set priors if not provided: use data-adaptive priors by default
  if (verbose) cat("Validating parameters...\n")
  if (is.null(priors)) {
    observed_mask <- data[[status_col]] == 1
    t_obs <- data[[time_col]][observed_mask]
    priors <- get_adaptive_priors(distribution, t_obs)
    if (verbose) {
      cat("Using data-adaptive priors based on observed events\n")
    }
  }
  
  # Merge user-provided MCMC options with sensible defaults
  mcmc_options <- normalize_mcmc_options(mcmc_options)
  
  # Prepare Stan data
  stan_data <- prepare_stan_data(data, time_col, status_col, distribution, priors)
  
  # Load and compile Stan model
  if (verbose) cat("Loading Stan model...\n")
  stan_model <- get_stan_model(distribution)
  
  # Initial values: Currently lets Stan auto-initialise 
  
  # Run MCMC sampling
  if (verbose) {
    cat("Running MCMC sampling...\n")
    cat("  Warmup iterations:", mcmc_options$iter_warmup, "\n")
    cat("  Sampling iterations:", mcmc_options$iter_sampling, "\n")
    cat("  Number of chains:", mcmc_options$chains, "\n")
  }
  
  start_time <- Sys.time()
  
  stan_fit <- rstan::sampling(
    object = stan_model,
    data = stan_data,
    iter = mcmc_options$iter_warmup + mcmc_options$iter_sampling,
    warmup = mcmc_options$iter_warmup,
    chains = mcmc_options$chains,
    cores = min(mcmc_options$chains, 4),
    refresh = if (verbose) max(mcmc_options$iter_sampling %/% 10, 1) else 0,
    control = build_rstan_control(mcmc_options)
  )
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  if (verbose) cat("Processing results...\n")
  
  # Extract MCMC results
  mcmc_summary <- extract_mcmc_summary(stan_fit, distribution)
  
  # Extract posterior samples for parameters (distribution-specific)
  param_names <- get_distribution_params(distribution)
  posterior_samples <- posterior::subset_draws(
    posterior::as_draws_df(stan_fit),
    variable = param_names
  )
  
  # Extract posterior distributions for each censored observation
  posterior_imputations <- extract_posterior_distributions(stan_fit, data_summary$n_censored)
  
  # Compute WAIC from pointwise log-likelihoods if available
  fit_metrics <- list()
  waic_res <- tryCatch({
    compute_waic_from_stan(stan_fit, data, time_col, status_col, distribution)
  }, error = function(e) NULL)
  if (!is.null(waic_res)) fit_metrics$waic <- waic_res
  
  # Generate multiple complete datasets by sampling from posteriors
  imputed_datasets <- generate_complete_datasets(data, time_col, status_col, 
                                               posterior_imputations, n_imputations)
  
  if (verbose) {
    cat("Imputation completed.\n")
    cat("Runtime:", round(runtime, 1), "seconds\n")
  }
  
  # Create result object with model_info structure expected by methods
  result <- structure(list(
    original_data = data,
    time_col = time_col,
    status_col = status_col,
    n_imputations = n_imputations,
    distribution = distribution,
    priors = priors,
    mcmc_options = mcmc_options,
    data_summary = data_summary,
    stan_fit = stan_fit,
    mcmc_summary = mcmc_summary,
    posterior_samples = posterior_samples,  # Parameter posterior draws
    imputed_datasets = imputed_datasets,
    posterior_imputations = posterior_imputations,  # Full posteriors for each censored obs
    runtime = runtime,
    model_info = list(
      distribution = distribution,
      n_imputations = n_imputations,
      runtime = runtime,
      data_summary = data_summary,
      mcmc_options = mcmc_options,
      time_unit = time_unit %||% attr(data, "time_unit") %||% "days"
    ),
    diagnostics = compute_diagnostics(stan_fit, distribution, mcmc_options = mcmc_options),
    fit_metrics = fit_metrics
  ), class = "bayesian_imputation")
  
  return(result)
}

# Helper function to extract posterior distributions for each censored observation
extract_posterior_distributions <- function(stan_fit, n_censored) {
  # Early exit if no censored observations 
  if (is.null(n_censored) || n_censored <= 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  
  out <- rstan::extract(stan_fit, pars = "t_imputed", permuted = TRUE)
  t_imputed <- out$t_imputed
  if (is.null(t_imputed)) {
    stop("t_imputed not found in Stan output")
  }

  # Expected shape: draws x n_censored
  if (is.vector(t_imputed)) {
    t_imputed <- matrix(as.numeric(t_imputed), ncol = 1)
  }
  if (!is.matrix(t_imputed)) {
    t_imputed <- as.matrix(t_imputed)
  }
  if (ncol(t_imputed) != n_censored) {
    stop("Unexpected t_imputed dimensions: expected ", n_censored, " columns, got ", ncol(t_imputed))
  }

  t_imputed
}

#' Generate Complete Datasets from Posterior Distributions
#'
#' This function generates multiple complete datasets by sampling from the posterior
#' distributions of censored observations.
#'
#' @param data Original data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column (1 = event, 0 = censored)
#' @param posterior_imputations Matrix of posterior draws for each censored observation
#'   (rows = MCMC draws, columns = censored observations)
#' @param n_imputations Number of complete datasets to generate
#'
#' @return List of n_imputations complete datasets
#' 
#' @details
#' For each censored observation in each dataset, the function independently samples 
#' from its posterior distribution. The imputed values are constrained to be >= the 
#' original censoring time.
#'
#' @examples
#' \dontrun{
#' # After running bayesian_impute
#' result <- bayesian_impute(data, "time", "status")
#' 
#' # Generate 50 additional datasets
#' more_datasets <- generate_complete_datasets(
#'   result$original_data,
#'   result$time_col,
#'   result$status_col,
#'   result$posterior_imputations,
#'   50
#' )
#' }
#' @export
generate_complete_datasets <- function(data, time_col, status_col, 
                                     posterior_imputations, n_imputations) {
  
  censored_indices <- which(data[[status_col]] == 0)
  n_censored <- length(censored_indices)
  
  # Early exit: if there are no censored observations, just replicate original data
  if (n_censored == 0) {
    datasets <- vector("list", n_imputations)
    for (i in seq_len(n_imputations)) {
      imputed_data <- data
      imputed_data$original_time <- data[[time_col]]
      imputed_data$original_status <- data[[status_col]]
      imputed_data$was_censored <- data[[status_col]] == 0
      imputed_data$dataset_id <- i
      imputed_data <- standardize_complete_column_order(
        imputed_data,
        time_col = time_col,
        status_col = status_col
      )
      datasets[[i]] <- imputed_data
    }
    return(datasets)
  }
  
  n_posterior_draws <- nrow(posterior_imputations)
  
  datasets <- list()
  
  for (i in 1:n_imputations) {
    # Create a copy of the original data
    imputed_data <- data
    
    # Add original value columns for transparency
    imputed_data$original_time <- data[[time_col]]
    imputed_data$original_status <- data[[status_col]]
    imputed_data$was_censored <- data[[status_col]] == 0
    
    # Each censored observation samples independently 
    # from its posterior 
    for (j in 1:n_censored) {
      # Sample a random MCMC draw for each censored observation
      random_draw_idx <- sample(1:n_posterior_draws, 1)
      
      # Extract imputed value from this observation's randomly chosen draw
      imputed_value <- posterior_imputations[random_draw_idx, j]
      
      # Ensure imputed time >= censoring time
      obs_idx <- censored_indices[j]
      imputed_value <- max(imputed_value, data[[time_col]][obs_idx])
      
      # Replace the censored time with the imputed time
      imputed_data[[time_col]][obs_idx] <- imputed_value
      # Mark as observed (event)
      imputed_data[[status_col]][obs_idx] <- 1
    }
    
    # Add dataset identifier
    imputed_data$dataset_id <- i
    imputed_data <- standardize_complete_column_order(
      imputed_data,
      time_col = time_col,
      status_col = status_col
    )
    datasets[[i]] <- imputed_data
  }
  
  return(datasets)
}

# Prepare Stan data for different distributions
prepare_stan_data <- function(data, time_col, status_col, distribution, priors) {
  
  # Separate observed and censored data
  observed_mask <- data[[status_col]] == 1
  censored_mask <- data[[status_col]] == 0
  
  t_obs <- data[[time_col]][observed_mask]
  t_cens <- data[[time_col]][censored_mask]
  
  # Base Stan data (common to all distributions)
  stan_data <- list(
    n_obs = length(t_obs),
    n_cens = length(t_cens),
    t_obs = t_obs,
    t_cens = t_cens
  )
  
  # Add distribution-specific prior parameters
  if (distribution == "weibull") {
    # Log-parameterisation priors
    stan_data$mu_log_shape <- priors$mu_log_shape
    stan_data$sd_log_shape <- priors$sd_log_shape
    stan_data$mu_log_scale <- priors$mu_log_scale
    stan_data$sd_log_scale <- priors$sd_log_scale
    # Weak but finite upper bounds
    t_all <- c(t_obs, t_cens)
    t_max <- max(t_all, na.rm = TRUE)
    shape_upper <- 50
    scale_upper <- 1e4 * t_max
    stan_data$shape_upper <- shape_upper
    stan_data$scale_upper <- scale_upper
  } else if (distribution == "exponential") {
    stan_data$rate_prior_shape <- priors$rate_prior_shape
    stan_data$rate_prior_rate <- priors$rate_prior_rate
  } else if (distribution == "lognormal") {
    stan_data$mu_prior_mean <- priors$mu_prior_mean
    stan_data$mu_prior_sd <- priors$mu_prior_sd
    stan_data$sigma_prior_sd <- priors$sigma_prior_sd
  } else {
    stop("Unknown distribution: ", distribution)
  }
  
  return(stan_data)
}

#' Compute MCMC diagnostics
#' @param mcmc_result Stan sampling result
#' @param distribution Distribution name to determine which parameters to check
#' @param mcmc_options Optional MCMC options list (used for max treedepth hits)
#' @return List of diagnostic information
#' @keywords internal
compute_diagnostics <- function(mcmc_result, distribution = "weibull", mcmc_options = NULL) {
  param_names <- get_distribution_params(distribution)

  draws <- posterior::subset_draws(
    posterior::as_draws_df(mcmc_result),
    variable = param_names
  )

  summ <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")

  max_rhat <- max(summ$rhat, na.rm = TRUE)
  min_ess_bulk <- min(summ$ess_bulk, na.rm = TRUE)
  min_ess_tail <- min(summ$ess_tail, na.rm = TRUE)

  sampler_params <- tryCatch(
    rstan::get_sampler_params(mcmc_result, inc_warmup = FALSE),
    error = function(e) NULL
  )

  num_divergent <- NA_integer_
  num_max_treedepth <- NA_integer_

  if (!is.null(sampler_params) && length(sampler_params) > 0) {
    num_divergent <- sum(vapply(sampler_params, function(m) {
      if (!("divergent__" %in% colnames(m))) return(0L)
      sum(m[, "divergent__"] > 0)
    }, integer(1)))

    if (!is.null(mcmc_options) && !is.null(mcmc_options$max_treedepth)) {
      max_td <- as.integer(mcmc_options$max_treedepth)
      num_max_treedepth <- sum(vapply(sampler_params, function(m) {
        if (!("treedepth__" %in% colnames(m))) return(0L)
        sum(m[, "treedepth__"] >= max_td)
      }, integer(1)))
    }
  }

  list(
    num_divergent = num_divergent,
    any_divergent = isTRUE(num_divergent > 0),
    num_max_treedepth = num_max_treedepth,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    convergence_ok = isTRUE(max_rhat <= 1.1) && isTRUE(min_ess_bulk >= 100)
  )
} 

#' Safe Group Validation
#'
#' Validates group variables and data structure before processing
#'
#' @param data Data frame containing survival data
#' @param groups Name of the group variable
#' @param time_col Name of the time column
#' @param status_col Name of the status column
#' @param verbose Print progress messages
#' @return List with validation results
#' @keywords internal
validate_groups_safe <- function(data, groups, time_col, status_col, verbose = TRUE) {
  
  # Check group variable exists
  if (!groups %in% names(data)) {
    return(list(valid = FALSE, error = paste("Group variable '", groups, "' not found")))
  }
  
  # Get group values safely
  groups_data <- data[[groups]]
  if (all(is.na(groups_data))) {
    return(list(valid = FALSE, error = "Group variable contains only NA values"))
  }
  
  # Handle different data types
  if (is.factor(groups_data)) {
    group_values <- levels(groups_data)
  } else {
    group_values <- sort(unique(groups_data[!is.na(groups_data)]))
    group_values <- as.character(group_values)
  }
  
  if (length(group_values) < 2) {
    return(list(valid = FALSE, error = "Need at least 2 groups for comparison"))
  }
  
  if (verbose) {
    cat("Group validation:\n")
    cat("  Variable:", groups, "\n")
    cat("  Groups:", paste(group_values, collapse = ", "), "\n")
  }
  
  return(list(
    valid = TRUE,
    group_values = group_values,
    groups_data = groups_data
  ))
}

#' Single Group Validation
#'
#' Validates a single group's data before processing
#'
#' @param group_data Data frame for a single group
#' @param time_col Name of the time column
#' @param status_col Name of the status column
#' @param group_name Name of the group being validated
#' @return TRUE if valid, stops with error if invalid
#' @keywords internal
validate_single_group_safe <- function(group_data, time_col, status_col, group_name) {
  
  n_total <- nrow(group_data)
  n_censored <- sum(group_data[[status_col]] == 0, na.rm = TRUE)
  n_observed <- sum(group_data[[status_col]] == 1, na.rm = TRUE)
  
  # Minimum requirements
  if (n_total < 5) {
    stop("Group '", group_name, "' has only ", n_total, " observations (need >= 5)")
  }
  if (n_observed < 2) {
    stop("Group '", group_name, "' has only ", n_observed, " events (need >= 2)")
  }
  if (n_censored == 0) {
    warning("Group '", group_name, "' has no censored observations")
  }
  
  # Check for NA values in key columns
  if (any(is.na(group_data[[time_col]]))) {
    stop("Group '", group_name, "' has NA values in time column")
  }
  if (any(is.na(group_data[[status_col]]))) {
    stop("Group '", group_name, "' has NA values in status column")
  }
  
  return(TRUE)
}

#' Safe Group Processing Function
#'
#' Processes groups individually with error isolation 
#'
#' @param data Data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column
#' @param groups Name of the group variable
#' @param n_imputations Number of imputations
#' @param distribution Distribution to use
#' @param priors Prior parameters
#' @param mcmc_options MCMC options
#' @param verbose Print progress messages
#' @return Group comparison object
#' @keywords internal
bayesian_impute_groups_safe <- function(data, time_col, status_col, groups,
                                      n_imputations, distribution, 
                                      priors, mcmc_options, verbose,
                                      time_unit = NULL) {
  
  # Validate groups before processing
  group_validation <- validate_groups_safe(data, groups, time_col, status_col, verbose)
  if (!group_validation$valid) {
    stop("Group validation failed: ", group_validation$error)
  }
  
  # Process each group individually
  group_results <- list()
  group_errors <- list()
  
  for (group_val in group_validation$group_values) {
    if (verbose) cat("Processing group:", group_val, "\n")
    
    # Each group gets its own tryCatch
    group_result <- tryCatch({
      
      # Extract group data
      group_data <- data[data[[groups]] == group_val, ]
      
      # Apply same status conversion as bayesian_impute_single 
      # (skip if already validated)
      if (!inherits(group_data, "survival_data")) {
        validation_result <- validate_survival_data(group_data, time_col, status_col)
        group_data_converted <- validation_result$data
      } else {
        group_data_converted <- group_data
      }
      
      # Validate the converted data
      validate_single_group_safe(group_data_converted, time_col, status_col, group_val)
      
      # Fit model using single-group logic
      group_result <- bayesian_impute_single(group_data_converted, time_col, status_col,
                                            n_imputations, distribution, 
                                            priors, mcmc_options, verbose = FALSE)  
      
      if (verbose) cat("Group", group_val, "completed successfully\n")
      group_result
      
    }, error = function(e) {
      # Don't let one group crash run
      error_msg <- paste("Error in group '", group_val, "': ", e$message, sep = "")
      if (verbose) cat("ERROR:", error_msg, "\n")
      return(list(error = error_msg, group = group_val))
    })
    
    # Store result
    if ("error" %in% names(group_result)) {
      group_errors[[group_val]] <- group_result$error
      if (verbose) cat("Stored error for group", group_val, "\n")
    } else {
      group_results[[group_val]] <- group_result
      if (verbose) cat("Stored result for group", group_val, "\n")
    }
  }
  
  # Check results
  if (verbose) {
    cat("Loop completed. Group results:", length(group_results), "Group errors:", length(group_errors), "\n")
    cat("Group names:", paste(names(group_results), collapse = ", "), "\n")
  }
  
  if (length(group_results) == 0) {
    stop("All groups failed. Errors:\n", 
         paste(unlist(group_errors), collapse = "\n"))
  }
  
  if (length(group_errors) > 0) {
    warning("Some groups failed:\n", 
            paste(unlist(group_errors), collapse = "\n"))
  }
  
  # Create group comparison object
  result <- structure(
    list(
      group_results = group_results,
      group_errors = group_errors,
      groups = groups,
      group_names = names(group_results),
      original_data = data,
      validation = group_validation
    ),
    class = c("bayesian_imputation_groups", "bayesian_imputation")
  )
  
  if (verbose) cat("Created group comparison object with class:", class(result), "\n")
  
  return(result)
} 


# Compute WAIC from Stan fit (requires pointwise log-likelihoods in Stan)
#' @keywords internal
compute_waic_from_stan <- function(stan_fit, data, time_col, status_col, distribution) {
  observed_mask <- data[[status_col]] == 1
  censored_mask <- data[[status_col]] == 0
  n <- nrow(data)
  
  ll_obs <- NULL
  ll_cens <- NULL

  # Extract pointwise ll arrays as draws x count matrices (if present)
  if (sum(observed_mask) > 0) {
    ll_obs <- tryCatch(
      rstan::extract(stan_fit, pars = "ll_obs", permuted = TRUE)$ll_obs,
      error = function(e) NULL
    )
  }
  if (sum(censored_mask) > 0) {
    ll_cens <- tryCatch(
      rstan::extract(stan_fit, pars = "ll_cens", permuted = TRUE)$ll_cens,
      error = function(e) NULL
    )
  }

  if (is.null(ll_obs) && is.null(ll_cens)) {
    stop("Pointwise log-likelihoods not found in Stan output")
  }

  # If one class is absent and the corresponding pointwise log-lik
  # is also unavailable from Stan, skip WAIC 
  if (sum(observed_mask) == 0 && is.null(ll_cens)) return(NULL)
  if (sum(censored_mask) == 0 && is.null(ll_obs)) return(NULL)
  
  # Combine into a draws x N matrix aligned to original rows
  S <- nrow(ll_obs) %||% nrow(ll_cens)
  if (is.null(S)) stop("No log-likelihood draws found")
  log_lik <- matrix(0, nrow = S, ncol = n)
  
  # Fill observed positions
  if (!is.null(ll_obs)) {
    log_lik[, observed_mask] <- ll_obs
  }
  # Fill censored positions
  if (!is.null(ll_cens)) {
    log_lik[, censored_mask] <- ll_cens
  }
  
  # WAIC calculations
  # lppd_i = log(mean_s exp(log_lik[s,i]))
  # p_waic = sum_i var_s(log_lik[s,i])
  # WAIC = -2 * (sum_i lppd_i - p_waic)
  maxlog <- apply(log_lik, 2, max)
  # stabilise exp via subtract-max 
  exp_centered <- exp(sweep(log_lik, 2, maxlog, FUN = "-"))
  mean_exp <- colMeans(exp_centered)
  lppd_i <- log(mean_exp) + maxlog
  p_waic_i <- apply(log_lik, 2, stats::var)
  lppd <- sum(lppd_i)
  p_waic <- sum(p_waic_i)
  waic <- -2 * (lppd - p_waic)
  
  list(
    waic = waic,
    lppd = lppd,
    p_waic = p_waic
  )
} 
