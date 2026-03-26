#' Export Imputed Datasets
#'
#' Export completed datasets from Bayesian imputation results in CSV or RDS
#' formats. Provides user-friendly options for saving imputed data for future
#' analysis.
#'
#' @param object A bayesian_imputation object (from bayesian_impute or bayes_np_impute)
#' @param file_path Base file path 
#' @param format Output format for saving ("csv", "rds")
#' @param datasets Which datasets to save ("all", "first", or specific numbers)
#' @param groups How to handle groups ("combined", "separate", or specific group name)
#' @param include_original Whether to include original value columns
#' @param combine If TRUE, write a single combined file across imputations
#' @param ... Additional arguments passed to writing functions
#'
#' @return Invisible list with file paths of saved files
#'
#' @examples
#' \dontrun{
#' # Load and impute data
#' data(lung, package = "survival")
#' result <- bayesian_impute(lung, n_imputations = 5)
#'
#' # Export all datasets as CSV files
#' export(result, "lung_imputed", format = "csv")
#'
#' # Export as RDS for future R analysis
#' export(result, "lung_imputed", format = "rds")
#'
#' # Export only first dataset
#' export(result, "lung_imputed", format = "csv", datasets = "first")
#'
#' # For group analysis
#' result_groups <- bayesian_impute(data, groups = "treatment")
#' export(result_groups, "group_imputed", format = "csv", groups = "combined")
#' }
#'
#' @export
export <- function(object, file_path, 
                  format = c("csv", "rds"),
                  datasets = "all",
                  groups = "combined",
                  include_original = TRUE,
                  combine = FALSE,
                  ...) {
  
  format <- match.arg(format)
  
  # Validate input
  if (!inherits(object, "bayesian_imputation")) {
    stop("object must be a bayesian_imputation object")
  }
  
  # Determine available dataset count 
  if (inherits(object, "bayesian_imputation_groups")) {
    # For group objects, check the first successful group
    successful_groups <- names(object$group_results)
    if (length(successful_groups) == 0) {
      stop("No successful groups found. Check group_errors for details.")
    }
    max_dataset <- length(object$group_results[[successful_groups[1]]]$imputed_datasets)
  } else {
    max_dataset <- length(object$imputed_datasets)
  }

  # Determine which datasets to save 
  if (identical(datasets, "all")) {
    dataset_numbers <- if (max_dataset > 0) seq_len(max_dataset) else integer(0)
  } else if (identical(datasets, "first")) {
    dataset_numbers <- 1
  } else if (is.numeric(datasets)) {
    dataset_numbers <- datasets
  } else {
    stop("datasets must be 'all', 'first', or a numeric vector")
  }

  if (length(dataset_numbers) == 0) {
    stop("No imputed datasets available to export.")
  }
  if (max(dataset_numbers) > max_dataset) {
    stop("Requested dataset ", max(dataset_numbers), " but only ", max_dataset, " available")
  }
  
  # Handle group objects
  if (inherits(object, "bayesian_imputation_groups")) {
    return(export_group_imputed_data(object, file_path, format, datasets, groups, include_original, combine, ...))
  }
  
  # Get completed datasets
  if (length(dataset_numbers) == 1) {
    # Single dataset
    completed_data <- complete(object, dataset = dataset_numbers, format = "wide")
    if (!include_original) {
      completed_data <- remove_original_columns(completed_data)
    }
  } else {
    # Multiple datasets
    completed_data <- complete(object, dataset = dataset_numbers, format = "long")
    if (!include_original) {
      completed_data <- remove_original_columns(completed_data)
    }
  }
  
  # Export based on format
  exported_files <- switch(format,
    "csv" = export_csv(completed_data, file_path, dataset_numbers, combine = combine, ...),
    "rds" = export_rds(completed_data, file_path, dataset_numbers, combine = combine, ...)
  )
  
  # Print summary
  cat("Exported imputed datasets:\n")
  for (file in exported_files) {
    cat("  ", file, "\n")
  }
  
  invisible(exported_files)
}

sanitize_filename_component <- function(x) {
  safe <- gsub("[^A-Za-z0-9._-]+", "_", x)
  safe <- gsub("^_+|_+$", "", safe)
  if (!nzchar(safe)) "group" else safe
}

resolve_dataset_numbers <- function(data, dataset_numbers) {
  if (is.character(dataset_numbers)) {
    if (identical(dataset_numbers, "first")) return(1L)
    if (identical(dataset_numbers, "all")) {
      if (".imp" %in% names(data)) {
        return(sort(unique(data$.imp)))
      }
      return(1L)
    }
    stop("datasets must be 'all', 'first', or a numeric vector")
  }
  if (is.numeric(dataset_numbers)) {
    return(dataset_numbers)
  }
  stop("datasets must be 'all', 'first', or a numeric vector")
}

#' Export Group Imputed Data
#' @keywords internal
export_group_imputed_data <- function(object, file_path, format, datasets, groups, include_original, combine, ...) {
  
  if (identical(groups, "combined")) {
    # Get combined data
    if (identical(datasets, "first")) {
      completed_data <- complete(object, dataset = 1, format = "wide")
    } else {
      completed_data <- complete(object, format = "long")
    }
    
    if (!include_original) {
      completed_data <- remove_original_columns(completed_data)
    }
    
    # Export combined data
    exported_files <- switch(format,
      "csv" = export_csv(completed_data, file_path, datasets, suffix = "_combined", combine = combine, ...),
      "rds" = export_rds(completed_data, file_path, datasets, suffix = "_combined", combine = combine, ...)
    )
    
  } else if (identical(groups, "separate")) {
    # Export each group separately
    exported_files <- character()
    
    for (group_name in object$group_names) {
      group_safe <- sanitize_filename_component(group_name)
      group_file_path <- paste0(file_path, "_", group_safe)
      
      if (identical(datasets, "first")) {
        group_data <- complete(object, dataset = 1, groups = group_name)
      } else {
        group_data <- complete(object, groups = group_name, format = "long")
      }
      
      if (!include_original) {
        group_data <- remove_original_columns(group_data)
      }
      
      group_files <- switch(format,
        "csv" = export_csv(group_data, group_file_path, datasets, suffix = paste0("_", group_safe), combine = combine, ...),
        "rds" = export_rds(group_data, group_file_path, datasets, suffix = paste0("_", group_safe), combine = combine, ...)
      )
      
      exported_files <- c(exported_files, group_files)
    }
    
  } else {
    # Specific group
    if (!groups %in% object$group_names) {
      stop("Group '", groups, "' not found. Available: ", paste(object$group_names, collapse = ", "))
    }
    
    group_safe <- sanitize_filename_component(groups)
    group_file_path <- paste0(file_path, "_", group_safe)
    
    if (identical(datasets, "first")) {
      group_data <- complete(object, dataset = 1, groups = groups)
    } else {
      group_data <- complete(object, groups = groups, format = "long")
    }
    
    if (!include_original) {
      group_data <- remove_original_columns(group_data)
    }
    
    exported_files <- switch(format,
      "csv" = export_csv(group_data, group_file_path, datasets, suffix = paste0("_", group_safe), combine = combine, ...),
      "rds" = export_rds(group_data, group_file_path, datasets, suffix = paste0("_", group_safe), combine = combine, ...)
    )
  }
  
  invisible(exported_files)
}

#' Remove Original Value Columns
#' @keywords internal
remove_original_columns <- function(data) {
  cols_to_remove <- if ("imputed_time" %in% names(data)) {
    c("time", "original_time", "original_status", "was_censored")
  } else {
    c("original_time", "original_status", "was_censored")
  }
  cols_to_remove <- cols_to_remove[cols_to_remove %in% names(data)]
  if (length(cols_to_remove) > 0) {
    data <- data[, !names(data) %in% cols_to_remove, drop = FALSE]
  }
  data
}

#' Export as CSV
#' @keywords internal
export_csv <- function(data, file_path, dataset_numbers, suffix = "", combine = FALSE, ...) {
  if (length(dataset_numbers) == 1 || isTRUE(combine)) {
    # Single combined file (either one dataset or combine=TRUE)
    file_name <- paste0(file_path, suffix, ".csv")
    utils::write.csv(data, file_name, row.names = FALSE, ...)
    return(file_name)
  } else {
    # Multiple files
    file_names <- character()
    for (i in dataset_numbers) {
      dataset_data <- data[data$.imp == i, ]
      file_name <- paste0(file_path, "_dataset_", i, suffix, ".csv")
      utils::write.csv(dataset_data, file_name, row.names = FALSE, ...)
      file_names <- c(file_names, file_name)
    }
    return(file_names)
  }
}

#' Export as RDS
#' @keywords internal
export_rds <- function(data, file_path, dataset_numbers, suffix = "", combine = FALSE, ...) {
  dataset_numbers <- resolve_dataset_numbers(data, dataset_numbers)

  if (length(dataset_numbers) == 1 || isTRUE(combine)) {
    # Single file
    file_name <- paste0(file_path, suffix, ".rds")
    saveRDS(data, file_name, ...)
    return(file_name)
  } else {
    # Multiple files (one per dataset)
    if (!(".imp" %in% names(data))) {
      stop("Cannot split RDS output without a '.imp' column")
    }
    file_names <- character()
    for (i in dataset_numbers) {
      dataset_data <- data[data$.imp == i, , drop = FALSE]
      file_name <- paste0(file_path, "_dataset_", i, suffix, ".rds")
      saveRDS(dataset_data, file_name, ...)
      file_names <- c(file_names, file_name)
    }
    return(file_names)
  }
}

#' Quick Export Function
#'
#' Convenience function for quick export of imputed datasets.
#'
#' @param object A bayesian_imputation object
#' @param file_path Base file path
#' @param format Output format ("csv", "rds")
#' @param time_unit Unused; reserved for backward compatibility
#' @param ... Additional arguments forwarded to `export()`
#'
#' @return Invisible file path
#'
#' @examples
#' \dontrun{
#' result <- bayesian_impute(data, n_imputations = 5)
#' quick_export(result, "my_imputed_data")
#' }
#'
#' @export
quick_export <- function(object, file_path, format = "csv", time_unit = "days", ...) {
  unused <- time_unit
  export(
    object = object,
    file_path = file_path,
    format = format,
    ...
  )
}
