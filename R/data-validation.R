# Data Validation and Processing Functions

#' Validate survival data input
#' @param data Data frame containing survival data
#' @param time_var Name of time variable
#' @param status_var Name of status variable
#' @param verbose Logical; print status-coding and censoring-pattern messages
#' @return List with validated and processed data
#' @keywords internal
validate_survival_data <- function(data, time_var, status_var, verbose = TRUE) {
  
  # Check basic inputs
  checkmate::assert_data_frame(data, min.rows = 1)
  checkmate::assert_string(time_var)
  checkmate::assert_string(status_var)
  
  # Check required columns exist
  if (!time_var %in% names(data)) {
    stop("Time variable '", time_var, "' not found in data")
  }
  
  if (!status_var %in% names(data)) {
    stop("Status variable '", status_var, "' not found in data")
  }
  
  # Extract variables
  time <- data[[time_var]]
  status <- data[[status_var]]
  
  # Validate time variable
  if (!is.numeric(time)) {
    stop("Time variable must be numeric")
  }
  
  if (any(time <= 0, na.rm = TRUE)) {
    stop("All survival times must be positive")
  }
  
  if (any(is.na(time))) {
    stop("Missing values in time variable not allowed")
  }
  
  # Validate status variable
  if (!is.numeric(status) && !is.logical(status)) {
    stop("Status variable must be numeric or logical")
  }
  
  # Convert status to numeric if logical
  if (is.logical(status)) {
    status <- as.numeric(status)
  }
  
  # Check status values and handle different coding schemes
  unique_status <- unique(status[!is.na(status)])
  
  if (any(is.na(status))) {
    stop("Missing values in status variable not allowed")
  }
  
  # Handle different coding conventions
  if (all(unique_status %in% c(0, 1))) {
    # Standard 0/1 coding - no conversion needed
    if (verbose) {
      message("Status coding detected: 0 = censored, 1 = event")
    }
  } else if (all(unique_status %in% c(1, 2))) {
    # Common 1/2 coding - convert to 0/1
    if (verbose) {
      message("Status coding detected: 1 = censored, 2 = event")
      message("Converting to standard coding: 0 = censored, 1 = event")
    }
    status <- status - 1  # Convert 1,2 to 0,1
    data[[status_var]] <- status  # Update the original data
  } else {
    stop("Status variable must contain only:\n",
         "  - 0/1 coding (0 = censored, 1 = event), or\n", 
         "  - 1/2 coding (1 = censored, 2 = event)\n",
         "Found values: ", paste(sort(unique_status), collapse = ", "))
  }
  
  # Check for minimum events and censoring
  n_events <- sum(status == 1)
  n_censored <- sum(status == 0)
  
  if (n_events < 2) {
    stop("At least 2 observed events required for estimation")
  }
  
  if (n_censored < 1) {
    warning("No censored observations found. Consider using standard survival analysis methods.")
  }
  
  # Separate observed and censored times
  observed_times <- time[status == 1]
  censored_times <- time[status == 0]
  
  # Validate censoring pattern
  max_observed <- ifelse(length(observed_times) > 0, max(observed_times), 0)
  min_censored <- ifelse(length(censored_times) > 0, min(censored_times), Inf)
  
  if (min_censored <= max_observed) {
    if (verbose) {
      message("Note: Some censoring times occur before the largest observed time.")
    }
  }
  
  return(list(
    data = data,  # Return the potentially modified data
    n_total = length(time),
    n_observed = n_events,
    n_censored = n_censored,
    observed_times = observed_times,
    censored_times = censored_times,
    all_times = time,
    all_status = status,
    censoring_rate = n_censored / length(time)
  ))
}

#' Validate imputation parameters
#' @param n_imputations Number of imputations
#' @param distribution Distribution name
#' @param mcmc_options MCMC options list
#' @return Validated parameters list
#' @keywords internal
validate_imputation_params <- function(n_imputations, distribution, mcmc_options) {
  
  # Validate n_imputations
  checkmate::assert_int(n_imputations, lower = 1, upper = 1000)
  
  if (n_imputations < 10) {
    warning("Using fewer than 10 imputations may give unreliable results")
  }
  
  # Validate distribution
  checkmate::assert_string(distribution)
  supported_distributions <- c("weibull", "exponential", "lognormal")
  
  if (!distribution %in% supported_distributions) {
    stop("Distribution '", distribution, "' not supported. ",
         "Supported: ", paste(supported_distributions, collapse = ", "))
  }
  
  mcmc_options <- normalize_mcmc_options(mcmc_options)
  
  # Validate specific MCMC parameters
  checkmate::assert_int(mcmc_options$iter_warmup, lower = 100)
  checkmate::assert_int(mcmc_options$iter_sampling, lower = 100)
  checkmate::assert_int(mcmc_options$chains, lower = 1, upper = 8)
  checkmate::assert_number(mcmc_options$adapt_delta, lower = 0.5, upper = 0.99)
  checkmate::assert_int(mcmc_options$max_treedepth, lower = 5, upper = 20)
  
  return(list(
    n_imputations = n_imputations,
    distribution = distribution,
    mcmc_options = mcmc_options
  ))
}
