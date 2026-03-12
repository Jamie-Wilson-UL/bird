#' bayes_np_impute
#' ----------------
#' Non-parametric Bayesian imputation for right-censored survival data using a
#' *Linear Dependent Dirichlet-Process* model (originally from **DPpackage**).
#'
#' This is a convenience wrapper around the internal
#' [`dp_np_LDDPsurvival()`] function (a vendored copy of
#' `DPpackage::LDDPsurvival()`), making it as easy to call as
#' [`bayesian_impute()`].
#'
#' @param data A `data.frame` with survival information.
#' @param time_col Name of the survival time column.
#' @param status_col Name of the status column (1 = event, 0 = right-censored).
#' @param groups Optional grouping variable for stratified analysis
#' @param n_imputations Number of imputed complete datasets to generate (default: 10)
#' @param prior A named list of prior hyper-parameters (same structure as the
#'   original `LDDPsurvival` prior list).  If `NULL`, sensible defaults are
#'   constructed automatically. The default alpha prior is Gamma(10, 10), making
#'   the DP concentration parameter learnable.
#' @param mcmc A named list with MCMC settings (`nburn`, `nsave`, `nskip`,
#'   `ndisplay`).  Defaults to `nburn = 1000`, `nsave = 4000`, `nskip = 1`,
#'   `ndisplay = 100`.
#' @param grid Numeric vector of grid points where the posterior survival &
#'   density are evaluated.  Defaults to `seq(0.01, max(time)*1.1, length = 100)`.
#' @param verbose Logical; print progress?  (Default `TRUE`).
#' @param time_unit Optional label for the time unit stored in the result (default: "days")
#'
#' @return A `bayesian_imputation` object containing:
#' \describe{
#'   \item{original_data}{The original input data}
#'   \item{imputed_datasets}{List of n_imputations complete datasets}
#'   \item{posterior_imputations}{Matrix of posterior draws for each censored observation}
#'   \item{posterior_samples}{MCMC samples for model parameters}
#'   \item{model_info}{Information about the fitted model}
#'   \item{diagnostics}{MCMC convergence diagnostics}
#'   \item{fit}{The raw `LDDPsurvival` object for advanced inspection}
#' }
#'
#' @export
#' @examples
#' set.seed(1)
#' n  <- 50
#' t  <- rexp(n, 0.05)               # true times
#' c  <- rexp(n, 0.02)               # censor times
#' time   <- pmin(t, c)
#' status <- as.integer(t <= c)
#' df <- data.frame(time, status)
#'
#' # Default: alpha ~ Gamma(10, 10) (learnable)
#' res <- bayes_np_impute(df, "time", "status", n_imputations = 10, verbose = FALSE,
#'                        mcmc = list(nburn = 100, nsave = 200, nskip = 2, ndisplay = 100))
#' 
#' # View results
#' print(res)
#' plot(res)
#' 
#' # Access completed datasets
#' head(res$imputed_datasets[[1]])
#'
#' # Custom prior with fixed alpha = 1
#' custom_prior <- list(a0 = -1, b0 = -1, nu = 4, m0 = 0, 
#'                      S0 = matrix(100, 1, 1), psiinv = matrix(1, 1, 1),
#'                      tau1 = 6.01, taus1 = 6.01, taus2 = 2.01)
#' res_fixed <- bayes_np_impute(df, "time", "status", n_imputations = 5, 
#'                              prior = custom_prior, verbose = FALSE,
#'                              mcmc = list(nburn = 100, nsave = 200, nskip = 2, ndisplay = 100))
#'
#' @references Jara, A., Hanson, T., Quintana, F. A., Müller, P., & Rosner, G. L. (2011). DPpackage: Bayesian Semi- and Nonparametric Modeling in R. Journal of Statistical Software, 40(5), 1–30. https://doi.org/10.18637/jss.v040.i05
#'
#' @export
bayes_np_impute <- function(data,
                            time_col = NULL,
                            status_col = NULL,
                            groups = NULL,  # Optional groups parameter
                            n_imputations = 10,
                            prior = NULL,
                            mcmc = list(nburn = 1000, nsave = 4000, nskip = 1, ndisplay = 100),
                            grid = NULL,
                            verbose = TRUE,
                            time_unit = NULL) {
  
  # Track whether user explicitly passed columns
  arg_missing_time <- missing(time_col)
  arg_missing_status <- missing(status_col)
  
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
  } else {
    # No auto-prepare path; if user didn't pass columns and we're using defaults
    if (verbose && arg_missing_time && arg_missing_status &&
        (time_col %in% names(data)) && (status_col %in% names(data))) {
      cat(sprintf("Detected columns: time='%s', status='%s'\n", time_col, status_col))
    }
  }
  
  # Check that time_col and status_col are provided
  if (is.null(time_col) || is.null(status_col)) {
    stop("You must specify 'time_col' and 'status_col'. ",
         "Alternatively, use prepare_survival_data() first for auto-detection.")
  }
  
  if (is.null(groups)) {
    return(bayes_np_impute_single(data, time_col, status_col, 
                                 n_imputations, prior, mcmc, grid, verbose, time_unit = time_unit))
  } else {
    return(bayes_np_impute_groups_safe(data, time_col, status_col, groups,
                                      n_imputations, prior, mcmc, grid, verbose, time_unit = time_unit))
  }
}

#' Bayesian Nonparametric Imputation for Single Group (Internal)
#'
#' Original bayes_np_impute logic, now separated for single-group analysis.
#' This function is called internally by bayes_np_impute() when no groups are specified.
#'
#' @param data A data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column (1 = event, 0 = censored)
#' @param n_imputations Number of imputed complete datasets to generate (default: 10)
#' @param prior List of prior parameters for the LDDP model
#' @param mcmc List of MCMC sampling options (nburn, nsave, nskip, ndisplay)
#' @param grid Prediction grid for survival curves
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A bayesian_imputation object
#' @keywords internal
bayes_np_impute_single <- function(data,
                                   time_col = "time",
                                   status_col = "status",
                                   n_imputations = 10,
                                   prior = NULL,
                                   mcmc = list(nburn = 1000, nsave = 4000, nskip = 1, ndisplay = 100),
                                   grid = NULL,
                                   verbose = TRUE,
                                   time_unit = NULL) {
  
  if (verbose) {
    cat("Starting Bayesian nonparametric imputation for censored survival data...\n")
    cat("Model: Linear Dependent Dirichlet Process (LDDP)\n")
    cat("Number of imputations:", n_imputations, "\n")
  }
  
  # Validate inputs and prepare data (same as parametric version)
  # (skip if already validated)
  if (!inherits(data, "survival_data")) {
    if (verbose) cat("Validating input data...\n")
    validation_result <- validate_survival_data(data, time_col, status_col)
    data <- validation_result$data  # Use the potentially converted data
  } else {
    if (verbose) cat("Data already validated, skipping validation...\n")
  }
  
  stopifnot(is.data.frame(data))
  if (!all(c(time_col, status_col) %in% names(data))) {
    stop("time_col or status_col not found in data")
  }
  
  time   <- data[[time_col]]
  status <- data[[status_col]]
  if (any(time <= 0)) stop("Times must be positive for log-scale transformation")

  # Validate sample size requirements (consistent with group validation)
  n_total <- nrow(data)
  n_observed <- sum(status == 1, na.rm = TRUE)
  n_censored <- sum(status == 0, na.rm = TRUE)
  
  if (n_total < 5) {
    stop("Dataset has only ", n_total, " observations (need >= 5)")
  }
  if (n_observed < 2) {
    stop("Dataset has only ", n_observed, " events (need >= 2)")
  }
  if (n_censored == 0) {
    warning("Dataset has no censored observations")
  }

  # Data summary (same as parametric version)
  if (verbose) {
    cat("Data summary:\n")
    cat("  Total observations:", nrow(data), "\n")
    cat("  Events:", sum(status == 1), "\n")
    cat("  Censored:", sum(status == 0), "\n")
  }

  # --- build interval-censored matrix ---------------------------------------
  left  <- time
  right <- ifelse(status == 1, time, -999)   # right-censored coded as -999
  ymat  <- cbind(left, right)
  z     <- matrix(1, nrow(data), 1)          

  # --- default prior --------------------------------------------------------
  if (is.null(prior)) {
    prior <- list(
      a0 = 10, b0 = 10,               # alpha ~ Gamma(10, 10)
      nu = 4,
      m0 = 0,
      S0 = matrix(25, 1, 1),          
      psiinv = matrix(1, 1, 1),
      tau1 = 8.01, taus1 = 8.01, taus2 = 4.01  
    )
  }

  # --- prediction grid ------------------------------------------------------
  if (is.null(grid)) {
    grid <- seq(0.01, max(time) * 1.1, length.out = 100)
  }
  zpred <- matrix(1, 1, 1)  # survival for baseline covariate pattern

  if (verbose) {
    cat("Running MCMC sampling...\n")
    cat("  Burn-in iterations:", mcmc$nburn, "\n")
    cat("  Sampling iterations:", mcmc$nsave, "\n")
    cat("  Thinning interval:", mcmc$nskip, "\n")
  }

  start_time <- Sys.time()

  fit <- dp_np_LDDPsurvival(ymat ~ 1,
                            zpred   = zpred,
                            prior   = prior,
                            mcmc    = mcmc,
                            state   = NULL,
                            status  = TRUE,
                            grid    = grid)

  runtime <- difftime(Sys.time(), start_time, units = "secs")

  if (verbose) {
    cat("MCMC completed in", round(runtime, 2), "seconds\n")
    cat("Generating", n_imputations, "completed datasets...\n")
  }

  # Get the full posterior of imputed log-times
  y_posterior <- fit$save.state$ysave  # nsave x nrec matrix
  
  # Extract posterior distributions for censored observations only
  censored_indices <- which(status == 0)
  
  if (length(censored_indices) == 0) {
    # No censored observations - return original data for all imputations
    if (verbose) {
      cat("No censored observations found - returning original data for all imputations\n")
    }
    imputed_datasets <- replicate(n_imputations, data, simplify = FALSE)
    posterior_imputations <- matrix(numeric(0), nrow = nrow(y_posterior), ncol = 0)
  } else {
    posterior_imputations <- y_posterior[, censored_indices, drop = FALSE]
    
    # Generate multiple completed datasets
    imputed_datasets <- generate_complete_datasets_np(
      data, time_col, status_col, posterior_imputations, censored_indices, n_imputations
    )
  }
  
  # Create data summary
  data_summary <- list(
    n_total = nrow(data),
    n_observed = sum(status == 1),
    n_censored = sum(status == 0),
    time_range = range(time),
    censoring_rate = mean(status == 0)
  )
  
  # Create result object matching parametric approach structure
  result <- structure(list(
    original_data = data,
    time_col = time_col,
    status_col = status_col,
    posterior_samples = list(
      y_posterior = y_posterior,
      censored_indices = censored_indices
    ),
    imputed_datasets = imputed_datasets,
    posterior_imputations = posterior_imputations,
    runtime = runtime,
    model_info = list(
      distribution = "nonparametric_lddp",
      n_imputations = n_imputations,
      runtime = runtime,
      data_summary = data_summary,
      mcmc_options = list(
        chains = 1,  # DP approach uses single chain
        iter_warmup = mcmc$nburn,
        iter_sampling = mcmc$nsave
      ),
      prior = prior,
      time_unit = time_unit %||% attr(data, "time_unit") %||% "days"
    ),
    diagnostics = list(
      n_iterations = nrow(y_posterior),
      n_censored = length(censored_indices),
      convergence_ok = TRUE,
      convergence = "MCMC completed successfully"
    ),
    fit = fit  # Keep the original fit for advanced users
  ), class = "bayesian_imputation")
  
  if (verbose) {
    cat("Imputation completed successfully!\n")
  }
  
  return(result)
} 

#' Generate Complete Datasets from Nonparametric DP Posterior
#' 
#' This function generates multiple completed datasets by sampling from the 
#' posterior distribution of imputed times for censored observations.
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
#' @keywords internal
generate_complete_datasets_np <- function(data, time_col, status_col, 
                                        posterior_imputations, censored_indices, n_imputations) {
  
  n_censored <- length(censored_indices)
  
  # Early exit: no censored observations -> replicate original data
  if (n_censored == 0) {
    datasets <- vector("list", n_imputations)
    for (i in seq_len(n_imputations)) {
      imputed_data <- data
      imputed_data$original_time <- data[[time_col]]
      imputed_data$original_status <- data[[status_col]]
      imputed_data$was_censored <- data[[status_col]] == 0
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
    
    # Sample posterior 
    joint_draw_idx <- sample(1:n_posterior_draws, 1)
    
    # Impute censored observations 
    for (j in 1:n_censored) {
      
      # Get the imputed log-time and convert to actual time
      log_time_imputed <- posterior_imputations[joint_draw_idx, j]
      time_imputed <- exp(log_time_imputed)
      
      # Ensure imputed time is >= original censoring time
      original_censoring_time <- data[[time_col]][censored_indices[j]]
      time_imputed <- max(time_imputed, original_censoring_time)
      
      # Update the dataset
      imputed_data[[time_col]][censored_indices[j]] <- time_imputed
      imputed_data[[status_col]][censored_indices[j]] <- 1L
    }

    imputed_data <- standardize_complete_column_order(
      imputed_data,
      time_col = time_col,
      status_col = status_col
    )
    datasets[[i]] <- imputed_data
  }
  
  return(datasets)
} 

#' Bayesian Nonparametric Imputation for Multiple Groups (Internal)
#'
#' This function handles group-based nonparametric imputation by fitting separate
#' LDDP models for each group and combining the results.
#'
#' @param data A data frame containing survival data
#' @param time_col Name of the time column
#' @param status_col Name of the status column (1 = event, 0 = censored)
#' @param groups Name of the grouping variable
#' @param n_imputations Number of imputed complete datasets to generate (default: 10)
#' @param prior List of prior parameters for the LDDP model
#' @param mcmc List of MCMC sampling options (nburn, nsave, nskip, ndisplay)
#' @param grid Prediction grid for survival curves
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return A bayesian_imputation_groups object
#' @keywords internal
bayes_np_impute_groups_safe <- function(data, time_col, status_col, groups,
                                       n_imputations, prior, mcmc, grid, verbose,
                                       time_unit = NULL) {
  
  # Validate groups before any processing
  group_validation <- validate_groups_safe(data, groups, time_col, status_col, verbose)
  if (!group_validation$valid) {
    stop("Group validation failed: ", group_validation$error)
  }
  
  # Process each group with error isolation
  group_results <- list()
  group_errors <- list()
  
  for (group_val in group_validation$group_values) {
    if (verbose) cat("Processing group:", group_val, "\n")
    
    # Each group gets its own tryCatch
    group_result <- tryCatch({
      
      # Extract group data
      group_data <- data[data[[groups]] == group_val, ]
      
      # First, apply the same status conversion logic that bayes_np_impute_single uses
      # (skip if already validated)
      if (!inherits(group_data, "survival_data")) {
        validation_result <- validate_survival_data(group_data, time_col, status_col)
        group_data_converted <- validation_result$data
      } else {
        group_data_converted <- group_data
      }
      
      # Now validate the converted data
      validate_single_group_safe(group_data_converted, time_col, status_col, group_val)
      
      # Fit model using single-group logic
      group_result <- bayes_np_impute_single(group_data_converted, time_col, status_col,
                                            n_imputations, prior, mcmc, grid, verbose = FALSE,
                                            time_unit = time_unit)  # Reduce verbosity
      
      if (verbose) cat("Group", group_val, "completed successfully\n")
      group_result
      
    }, error = function(e) {
      # Don't let one group crash everything
      error_msg <- paste("Error in group '", group_val, "': ", e$message, sep = "")
      if (verbose) cat("ERROR:", error_msg, "\n")
      return(list(error = error_msg, group = group_val))
    })
    
    # Store result even if it's an error
    if ("error" %in% names(group_result)) {
      group_errors[[group_val]] <- group_result$error
      if (verbose) cat("Stored error for group", group_val, "\n")
    } else {
      group_results[[group_val]] <- group_result
      if (verbose) cat("Stored result for group", group_val, "\n")
    }
  }
  
  # Check if we have any successful results
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
