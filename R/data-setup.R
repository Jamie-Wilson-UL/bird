# Data Setup for Survival Analysis
# Functions to prepare survival data for Bayesian imputation

#' Prepare Survival Data for Bayesian Imputation
#'
#' This function prepares survival data for use with bayesian_impute().
#' It validates explicitly specified time and status columns and creates a
#' properly formatted survival dataset object.
#'
#' @param data A data frame containing survival data
#' @param time Column name or index for survival times.
#' @param status Column name or index for event status (1=event, 0=censored).
#' @param validate Logical; whether to run data validation (default: TRUE)
#' @param verbose Logical; print validation information (default: TRUE)
#' @param time_unit Character string describing the unit used for time values (default "days")
#'
#' @return A survival_data object ready for bayesian_impute()
#'
#' @details
#' Automatic column detection is not yet supported. Users must specify both
#' `time` and `status` so the package does not guess the wrong columns.
#'
#' The resulting object retains all original columns but adds metadata about
#' which columns represent time and status for use by bayesian_impute().
#'
#' @examples
#' \dontrun{
#' library(survival)
#' lung_data <- prepare_survival_data(lung, time = "time", status = "status")
#'
#' my_data <- prepare_survival_data(mydata, time = "duration", status = "died")
#'
#' # Then use directly with bayesian_impute
#' result <- bayesian_impute(lung_data, n_imputations = 20)
#' }
#'
#' @seealso \code{\link{bayesian_impute}}, \code{\link{validate_survival_data}}
#' @export
prepare_survival_data <- function(data, 
                                 time,
                                 status,
                                 validate = TRUE,
                                 verbose = TRUE,
                                 time_unit = "days") {
  
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame")
  }
  
  if (nrow(data) == 0) {
    stop("Data cannot be empty")
  }
  
  if (missing(time) || missing(status) || is.null(time) || is.null(status)) {
    stop("Please specify both 'time' and 'status'; auto-detection is not yet supported.\n",
         "Available columns: ", paste(names(data), collapse = ", "))
  }
  
  # Convert column names to character if they were provided as indices
  if (is.numeric(time)) {
    time <- names(data)[time]
  }
  if (is.numeric(status)) {
    status <- names(data)[status]
  }
  
  # Validate columns exist
  if (!time %in% names(data)) {
    stop("Time column '", time, "' not found in data")
  }
  if (!status %in% names(data)) {
    stop("Status column '", status, "' not found in data")
  }
  
  # Run validation if requested
  if (validate) {
    validation_result <- validate_survival_data(data, time, status, verbose = verbose)
    data <- validation_result$data  # Use the potentially converted data
  }
  
  # Create enhanced data object
  # Normalise time_unit to a simple label
  if (is.null(time_unit) || !nzchar(time_unit)) time_unit <- "days"
  time_unit <- tolower(as.character(time_unit))

  result <- structure(
    data,
    time_col = time,
    status_col = status,
    n_total = nrow(data),
    n_censored = sum(data[[status]] == 0),
    censoring_rate = mean(data[[status]] == 0),
    time_unit = time_unit,
    class = c("survival_data", class(data))
  )
  
  if (verbose) {
    message("Successfully prepared survival data with ", nrow(data), " observations")
    message("Censoring rate: ", round(attr(result, "censoring_rate") * 100, 1), "%")
    message("Time unit: ", attr(result, "time_unit"))
  }
  
  return(result)
}


#' Print method for survival_data objects, primarily an older debugging/validation tool
#' @param x A survival_data object
#' @param ... Additional arguments (ignored)
#' @export
print.survival_data <- function(x, ...) {
  cat("Survival Data\n")
  cat("=============\n")
  cat("Observations:", attr(x, "n_total"), "\n")
  cat("Time column: '", attr(x, "time_col"), "'\n", sep = "")
  cat("Status column: '", attr(x, "status_col"), "'\n", sep = "")
  cat("Censored:", attr(x, "n_censored"), 
      "(", round(attr(x, "censoring_rate") * 100, 1), "%)\n", sep = "")
  cat("\nFirst few rows:\n")
  print(head(as.data.frame(x)))
  
  invisible(x)
}

#' Quick Survival Analysis Setup
#'
#' One-liner function for users who want to go from raw data to imputation
#' with sensible defaults. Useful for beginners or quick analyses.
#'
#' @param data Raw survival data
#' @param groups Name of the group variable for group comparison (optional)
#' @param n_imputations Number of imputed datasets (default: 10)
#' @param distribution Survival distribution to use (default: "weibull")
#' @param time_unit Label describing the time unit to display (default "days")
#' @param ... Additional arguments passed to `impute()`. For raw data, this
#'   must include `time` and `status`.
#'
#' @return bayesian_imputation object
#' @export
quick_impute <- function(data, groups = NULL, n_imputations = 10, distribution = "weibull", time_unit = "days", ...) {
  
  # Delegate to unified wrapper
  result <- impute(
    data = data,
    groups = groups,
    n_imputations = n_imputations,
    distribution = distribution,
    time_unit = time_unit,
    ...
  )
  
  return(result)
}
