#' Extract completed datasets
#'
#' Generic function to extract completed datasets from a
#' `bayesian_imputation` or `bayesian_imputation_groups` object. See the
#' method documentation for details.
#'
#' @param object Object produced by [bayesian_impute()] or related helpers.
#' @param ... Passed to methods.
#' @export
complete <- function(object, ...) {
  UseMethod("complete")
}

# Complete Method for Bayesian Imputation Objects

#' Extract Complete Datasets from Bayesian Imputation
#'
#' @description
#' This method extracts complete datasets from a `bayesian_imputation`
#' object, providing simple access to imputed datasets without requiring
#' additional processing.
#'
#' @param object A bayesian_imputation object from `bayesian_impute()`
#' @param dataset Which dataset(s) to extract:
#'   \itemize{
#'     \item \code{"all"} (default): Return all imputed datasets as a list
#'     \item Integer (e.g., \code{1}): Return specific dataset number
#'     \item Vector (e.g., \code{c(1,3,5)}): Return specific datasets as a list
#'   }
#' @param format Output format:
#'   \itemize{
#'     \item \code{"wide"} (default): Each dataset as separate data.frame
#'     \item \code{"long"}: All datasets stacked with .imp and .id columns
#'     \item \code{"list"}: Always return as list (same as "wide" for multiple datasets)
#'   }
#' @param include.original Logical; include original data with missing values
#'   as dataset 0 (default: FALSE)
#' @param ... Additional arguments passed to downstream helpers
#'
#' @return 
#' \itemize{
#'   \item Single dataset: data.frame with imputed values and original value columns
#'   \item Multiple datasets: list of data.frames or single long data.frame
#' }
#'
#' @details
#' The completed datasets contain the following columns:
#' \itemize{
#'   \item \code{time}: Imputed survival time (or original time if was an event)
#'   \item \code{status}: Always 1 (all observations are events in completed datasets)
#'   \item \code{original_time}: Original time value (censoring time for censored observations)
#'   \item \code{original_status}: Original status (0 = censored, 1 = event)
#'   \item \code{was_censored}: Logical flag indicating originally censored observations
#'   \item \code{dataset_id}: Dataset identifier (in wide format)
#'   \item \code{.imp}: Imputation number (in long format)
#'   \item \code{.id}: Row identifier (in long format)
#'   \item All original covariates: Preserved unchanged
#' }
#'
#' For group analyses, the \code{groups} parameter controls how group data is returned:
#' \itemize{
#'   \item \code{"combined"} (default): Merges all groups into a single dataset with group variable preserved
#'   \item \code{"separate"}: Returns separate datasets for each group as a list
#'   \item \code{group_name}: Returns data for a specific group only
#' }
#'
#' @details
#' This function provides an intuitive interface for extracting completed datasets,
#' eliminating the need for the more complex `generate_complete_datasets()` function
#' for simple use cases.
#'
#' @examples
#' \dontrun{
#' # Basic imputation
#' imp <- bayesian_impute(lung, "time", "status", n_imputations = 5)
#' 
#' # Extract first dataset
#' dataset1 <- complete(imp, dataset = 1)
#' 
#' # Extract all datasets
#' all_datasets <- complete(imp)  # or complete(imp, "all")
#' 
#' # Extract specific datasets
#' some_datasets <- complete(imp, dataset = c(1, 3, 5))
#' 
#' # Get long format (all datasets stacked)
#' long_data <- complete(imp, format = "long")
#' 
#' # Include original data with missing values
#' with_original <- complete(imp, include.original = TRUE)
#' }
#'
#' @seealso \code{\link{bayesian_impute}}, \code{\link{generate_complete_datasets}}
#' @export
complete.bayesian_imputation <- function(object, 
                                       dataset = "all", 
                                       format = c("wide", "long", "list"), 
                                       include.original = FALSE,
                                       ...) {
  
  format <- match.arg(format)
  
  # Validate inputs
  if (!inherits(object, "bayesian_imputation")) {
    stop("Object must be of class 'bayesian_imputation'")
  }
  
  if (is.null(object$imputed_datasets) || length(object$imputed_datasets) == 0) {
    stop("No imputed datasets found in object. Run bayesian_impute() first.")
  }

  time_col <- if (!is.null(object$time_col)) object$time_col else "time"
  status_col <- if (!is.null(object$status_col)) object$status_col else "status"
  
  n_imputations <- length(object$imputed_datasets)
  
  # Handle dataset selection
  if (is.character(dataset) && dataset == "all") {
    selected_datasets <- object$imputed_datasets
    dataset_numbers <- seq_len(n_imputations)
  } else if (is.numeric(dataset)) {
    # Validate dataset numbers
    if (any(dataset < 1 | dataset > n_imputations)) {
      stop("Dataset numbers must be between 1 and ", n_imputations)
    }
    selected_datasets <- object$imputed_datasets[dataset]
    dataset_numbers <- dataset
  } else {
    stop("'dataset' must be 'all' or a numeric vector")
  }
  
  # Add original data if requested
  if (include.original) {
    original_with_id <- object$original_data
    original_with_id$.imp <- 0
    original_with_id$.id <- seq_len(nrow(original_with_id))
    
    # Ensure original data has same columns as imputed datasets
    if (length(selected_datasets) > 0) {
      # Get columns from first imputed dataset
      imputed_cols <- names(selected_datasets[[1]])
      original_cols <- names(original_with_id)
      
      # Add missing columns to original data
      missing_cols <- setdiff(imputed_cols, original_cols)
      for (col in missing_cols) {
        original_with_id[[col]] <- NA  # Fill with NA for original data
      }
      
      # Reorder columns to match imputed datasets
      original_with_id <- original_with_id[, imputed_cols]
      original_with_id <- standardize_complete_column_order(
        original_with_id,
        time_col = time_col,
        status_col = status_col
      )
    }
    
    selected_datasets <- c(list(original_with_id), selected_datasets)
    dataset_numbers <- c(0, dataset_numbers)
  }
  
  # Handle different output formats
  if (format == "long") {
    # Stack all datasets with .imp and .id columns
    long_data <- do.call(rbind, lapply(seq_along(selected_datasets), function(i) {
      d <- selected_datasets[[i]]
      if (!(".imp" %in% names(d))) {
        d$.imp <- dataset_numbers[i]
      }
      if (!(".id" %in% names(d))) {
        d$.id <- seq_len(nrow(d))
      }
      return(d)
    }))
    
    # Reorder columns to put .imp and .id first
    long_data <- standardize_complete_column_order(
      long_data,
      time_col = time_col,
      status_col = status_col
    )
    
    return(long_data)
    
  } else if (format %in% c("wide", "list")) {
    # Return as list or single dataset
    
    # Add .imp column to each dataset for consistency
    for (i in seq_along(selected_datasets)) {
      if (!(".imp" %in% names(selected_datasets[[i]]))) {
        selected_datasets[[i]]$.imp <- dataset_numbers[i]
      }
      selected_datasets[[i]] <- standardize_complete_column_order(
        selected_datasets[[i]],
        time_col = time_col,
        status_col = status_col
      )
    }
    
    # If single dataset requested, return data.frame; otherwise return list
    if (length(selected_datasets) == 1 && length(dataset) == 1 && is.numeric(dataset)) {
      return(selected_datasets[[1]])
    } else {
      names(selected_datasets) <- paste0("dataset_", dataset_numbers)
      return(selected_datasets)
    }
  }
}

 
