#' Calculate Clinical Survival Metrics
#' 
#' Helper functions to compute clinically meaningful survival statistics
#' from bayesian_imputation objects
#' 
#' @keywords internal

resolve_time_unit_label <- function(x, fallback = "days") {
  if (!is.null(x) && !is.null(x$model_info) && !is.null(x$model_info$time_unit)) {
    return(x$model_info$time_unit)
  }
  if (!is.null(x) && inherits(x, "survival_data")) {
    tu <- attr(x, "time_unit")
    if (!is.null(tu) && nzchar(tu)) return(tu)
  }
  fallback
}

default_clinical_time_points <- function(time_unit = "days") {
  unit <- tolower(time_unit %||% "days")
  if (unit %in% c("month", "months")) return(c(12, 24, 60))
  if (unit %in% c("year", "years")) return(c(1, 2, 5))
  c(365, 730, 1825)
}

default_event_probability_times <- function(time_unit = "days") {
  unit <- tolower(time_unit %||% "days")
  if (unit %in% c("month", "months")) return(c(6, 12))
  if (unit %in% c("year", "years")) return(c(0.5, 1))
  c(180, 365)
}

#' Calculate imputation summary comparing censored vs imputed values
#' @param x bayesian_imputation object
#' @return List with imputation metrics
#' @keywords internal
calculate_imputation_summary <- function(x) {
  
  # Input validation
  if (!inherits(x, "bayesian_imputation")) {
    stop("Input must be a bayesian_imputation object")
  }
  
  time_var <- x$time_col
  status_var <- x$status_col
  
  # Original data
  original_data <- x$original_data
  
  # Find censored observations in original data
  censored_mask <- original_data[[status_var]] == 0
  original_censored_times <- original_data[[time_var]][censored_mask]
  n_censored <- sum(censored_mask)
  
  if (n_censored == 0) {
    return(list(
      n_censored = 0,
      has_imputations = FALSE,
      message = "No censored observations to impute"
    ))
  }
  
  # Check if imputed datasets exist
  if (is.null(x$imputed_datasets) || length(x$imputed_datasets) == 0) {
    return(list(
      n_censored = n_censored,
      has_imputations = FALSE,
      message = "No imputed datasets available"
    ))
  }
  
  # Get imputed times from completed datasets
  n_imputations <- length(x$imputed_datasets)
  imputed_times_matrix <- matrix(NA, nrow = n_censored, ncol = n_imputations)
  
  for (i in 1:n_imputations) {
    # Get times for censored observations in this completed dataset
    imputed_times_matrix[, i] <- x$imputed_datasets[[i]][[time_var]][censored_mask]
  }
  
  # Calculate summary statistics
  imputed_times_mean <- rowMeans(imputed_times_matrix)  # Mean across imputations for each observation
  imputation_gains <- imputed_times_mean - original_censored_times  # How much longer they lived
  
  # Overall statistics
  original_censored_summary <- c(
    mean = mean(original_censored_times),
    median = median(original_censored_times),
    min = min(original_censored_times),
    max = max(original_censored_times)
  )
  
  imputed_times_summary <- c(
    mean = mean(imputed_times_mean),
    median = median(imputed_times_mean), 
    min = min(imputed_times_mean),
    max = max(imputed_times_mean)
  )
  
  imputation_gains_summary <- c(
    mean = mean(imputation_gains),
    median = median(imputation_gains),
    min = min(imputation_gains),
    max = max(imputation_gains)
  )
  
  # Quality checks
  all_gains_positive <- all(imputation_gains > 0)
  sample_imputed_values <- if (n_censored >= 5) imputed_times_mean[1:5] else imputed_times_mean
  sample_original_values <- if (n_censored >= 5) original_censored_times[1:5] else original_censored_times
  sample_gains <- if (n_censored >= 5) imputation_gains[1:5] else imputation_gains
  
  return(list(
    n_censored = n_censored,
    has_imputations = TRUE,
    original_censored_summary = original_censored_summary,
    imputed_times_summary = imputed_times_summary,
    imputation_gains_summary = imputation_gains_summary,
    all_gains_positive = all_gains_positive,
    sample_original_values = sample_original_values,
    sample_imputed_values = sample_imputed_values,
    sample_gains = sample_gains,
    n_imputations = n_imputations
  ))
}

#' Calculate survival statistics for original and imputed data
#' @param x bayesian_imputation object
#' @param time_points Optional vector of time points for survival probability
#'   estimates. If `NULL`, defaults are chosen using `time_unit`.
#' @return List with clinical metrics
#' @keywords internal
calculate_clinical_metrics <- function(x, time_points = NULL) {
  
  # Input validation
  if (!inherits(x, "bayesian_imputation")) {
    stop("Input must be a bayesian_imputation object")
  }
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for clinical metrics")
  }

  if (is.null(time_points)) {
    time_unit <- resolve_time_unit_label(x, fallback = "days")
    time_points <- default_clinical_time_points(time_unit)
  }
  
  time_var <- x$time_col
  status_var <- x$status_col
  
  # Original data and observed events only 
  original_data <- x$original_data
  original_observed <- original_data[original_data[[status_var]] == 1, ]
  
  # Calculate original survival statistics using KM estimates
  if (nrow(original_observed) > 0) {
    # KM survival curve
    km_original <- survival::survfit(
      formula = as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
      data = original_data
    )
    
    # KM median survival
    km_med <- suppressWarnings(as.numeric(km_original$table["median"]))
    if (length(km_med) == 0 || is.na(km_med) || !is.finite(km_med)) {
      # Fallback: first time KM drops to <= 0.5
      idx <- which(km_original$surv <= 0.5)[1]
      km_med <- if (!is.na(idx)) km_original$time[idx] else NA_real_
    }
    original_median <- km_med
    
    # KM-based mean and quantiles
    # Mean: area under the KM curve
    original_mean <- tryCatch({
      if (length(km_original$time) > 1) {
        # Integrate the survival function to get mean
        times <- c(0, km_original$time)
        surv <- c(1, km_original$surv)
        # Use trapezoidal rule
        sum(diff(times) * (surv[-1] + surv[-length(surv)]) / 2)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)
    
    # Quantiles: find times where survival drops to 0.75, 0.25
    original_q25 <- tryCatch({
      idx <- which(km_original$surv <= 0.75)[1]
      if (!is.na(idx)) km_original$time[idx] else NA_real_
    }, error = function(e) NA_real_)
    
    original_q75 <- tryCatch({
      idx <- which(km_original$surv <= 0.25)[1]
      if (!is.na(idx)) km_original$time[idx] else NA_real_
    }, error = function(e) NA_real_)
    
    # Survival probabilities at time points
    original_surv_probs <- summary(km_original, times = time_points, extend = TRUE)$surv
  } else {
    original_median <- original_mean <- original_q25 <- original_q75 <- NA
    original_surv_probs <- rep(NA, length(time_points))
  }
  
  # Check if imputed datasets exist
  if (is.null(x$imputed_datasets) || length(x$imputed_datasets) == 0) {
    return(list(
      original = list(
        median = original_median,
        mean = original_mean,
        q25 = original_q25,
        q75 = original_q75,
        survival_probs = original_surv_probs
      ),
      imputed = list(
        median_summary = c(mean = NA, sd = NA, q25 = NA, q75 = NA),
        survival_summary = matrix(NA, nrow = 4, ncol = length(time_points)),
        individual_medians = numeric(0)
      ),
      time_points = time_points,
      interpretation = list(message = "No imputed datasets available"),
      change = list(
        median_change = NA,
        median_change_percent = NA
      )
    ))
  }
  
  # Imputed data statistics (from completed datasets)
  n_imputations <- length(x$imputed_datasets)
  imputed_medians <- numeric(n_imputations)
  imputed_means <- numeric(n_imputations)
  imputed_surv_probs <- matrix(NA, nrow = n_imputations, ncol = length(time_points))
  
  for (i in 1:n_imputations) {
    dataset <- x$imputed_datasets[[i]]
    
    # KM-based median survival per completed dataset
    km_imputed <- survival::survfit(
      formula = as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
      data = dataset
    )
    med_idx <- which(km_imputed$surv <= 0.5)[1]
    imputed_medians[i] <- if (!is.na(med_idx)) km_imputed$time[med_idx] else NA_real_
    
    # Mean time 
    imputed_means[i] <- mean(dataset[[time_var]])
    
    # Survival probabilities at time points
    surv_summary <- summary(km_imputed, times = time_points, extend = TRUE)
    imputed_surv_probs[i, ] <- surv_summary$surv
  }
  
  # Summary statistics across imputations 
  good <- is.finite(imputed_medians)
  if (sum(good) == 0) {
    imputed_median_summary <- c(
      mean = NA_real_, sd = NA_real_,
      q025 = NA_real_, q975 = NA_real_,
      q25 = NA_real_, q75 = NA_real_
    )
  } else {
    imputed_median_summary <- c(
      mean = mean(imputed_medians[good], na.rm = TRUE),
      sd   = sd(imputed_medians[good],  na.rm = TRUE),
      q025 = unname(stats::quantile(imputed_medians[good], 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(imputed_medians[good], 0.975, na.rm = TRUE)),
      q25  = unname(stats::quantile(imputed_medians[good], 0.25,  na.rm = TRUE)),
      q75  = unname(stats::quantile(imputed_medians[good], 0.75,  na.rm = TRUE))
    )
  }
  
  imputed_surv_summary <- apply(imputed_surv_probs, 2, function(x) {
    c(mean = mean(x), sd = sd(x), q25 = quantile(x, 0.25), q75 = quantile(x, 0.75))
  })
  
  return(list(
    original = list(
      median = original_median,
      mean = original_mean,
      q25 = original_q25,
      q75 = original_q75,
      survival_probs = original_surv_probs
    ),
    imputed = list(
      median_summary = imputed_median_summary,
      survival_summary = imputed_surv_summary,
      individual_medians = imputed_medians
    ),
    time_points = time_points,
    change = list(
      median_change = imputed_median_summary["mean"] - original_median,
      median_change_percent = if (!is.na(original_median) && original_median > 0) {
        ((imputed_median_summary["mean"] - original_median) / original_median) * 100
      } else NA
    )
  ))
}


#' Format survival probability for display
#' @param prob Survival probability
#' @return Formatted string with percentage
#' @keywords internal
format_survival_prob <- function(prob) {
  if (is.na(prob) || is.null(prob)) return("N/A")
  if (prob < 0 || prob > 1) return("Invalid")
  paste0(round(prob * 100, 1), "%")
}

#' Format survival time with appropriate units
#' @param time Time value
#' @param decimals Number of decimal places
#' @param unit Optional unit label; defaults to inferred unit from `x`
#' @param x Optional object used to infer the time unit
#' @return Formatted string
#' @keywords internal
format_survival_time <- function(time, decimals = 1, unit = NULL, x = NULL) {
  if (is.na(time) || is.null(time)) return("N/A")
  if (!is.numeric(time) || time < 0) return("Invalid")

  # Resolve unit: explicit arg > object info > default
  resolved_unit <- unit
  if (is.null(resolved_unit) && !is.null(x) && !is.null(x$model_info) && !is.null(x$model_info$time_unit)) {
    resolved_unit <- x$model_info$time_unit
  }
  if (is.null(resolved_unit) && inherits(x, "survival_data")) {
    resolved_unit <- attr(x, "time_unit")
  }
  if (is.null(resolved_unit)) resolved_unit <- "days"

  paste(round(time, decimals), resolved_unit)
} 

#' Event probability by clinical time horizons
#'
#' Computes the probability of experiencing the event by specified time horizons
#' for original observed events and across completed imputed datasets.
#'
#' @param x A `bayesian_imputation` object
#' @param times Optional numeric vector of time horizons (in the same units as
#'   the data). If `NULL`, defaults are chosen using `time_unit`.
#' @param method Estimation method: "logspline" (default) or "ecdf" fallback
#' @return A data.frame with columns: `time`, `original_prob`, `imputed_mean`,
#'   `imputed_sd`, `imputed_q025`, `imputed_q975`, `n_imputations`
#' @examples
#' # event_probability(result, times = c(180, 365))
#' @export
event_probability <- function(x, times = NULL, method = c("logspline", "ecdf")) {
  method <- match.arg(method)

  if (!inherits(x, "bayesian_imputation")) {
    stop("x must be a bayesian_imputation object")
  }

  if (is.null(times)) {
    time_unit <- resolve_time_unit_label(x, fallback = "days")
    times <- default_event_probability_times(time_unit)
  }

  time_var <- x$time_col
  status_var <- x$status_col

  # Original observed events only
  original_times <- x$original_data[[time_var]][x$original_data[[status_var]] == 1]

  # Function to get CDF at times for a vector of times
  get_cdf <- function(times_vec, eval_times) {
    if (method == "logspline") {
      if (!requireNamespace("logspline", quietly = TRUE)) {
        warning("logspline not available; falling back to ECDF")
        Fhat <- stats::ecdf(times_vec)
        return(sapply(eval_times, Fhat))
      }
      fit <- logspline::logspline(times_vec)
      return(logspline::plogspline(eval_times, fit))
    } else {
      Fhat <- stats::ecdf(times_vec)
      return(sapply(eval_times, Fhat))
    }
  }

  # Original probabilities
  original_prob <- if (length(original_times) > 0) get_cdf(original_times, times) else rep(NA_real_, length(times))

  # Imputed probabilities across completed datasets
  if (is.null(x$imputed_datasets) || length(x$imputed_datasets) == 0) {
    imputed_mat <- matrix(NA_real_, nrow = 0, ncol = length(times))
  } else {
    imputed_mat <- matrix(NA_real_, nrow = length(x$imputed_datasets), ncol = length(times))
    for (i in seq_along(x$imputed_datasets)) {
      ds_times <- x$imputed_datasets[[i]][[time_var]]
      imputed_mat[i, ] <- get_cdf(ds_times, times)
    }
  }

  imputed_mean <- if (nrow(imputed_mat) > 0) colMeans(imputed_mat, na.rm = TRUE) else rep(NA_real_, length(times))
  imputed_sd <- if (nrow(imputed_mat) > 1) apply(imputed_mat, 2, stats::sd, na.rm = TRUE) else rep(NA_real_, length(times))
  imputed_q025 <- if (nrow(imputed_mat) > 0) apply(imputed_mat, 2, stats::quantile, 0.025, na.rm = TRUE) else rep(NA_real_, length(times))
  imputed_q975 <- if (nrow(imputed_mat) > 0) apply(imputed_mat, 2, stats::quantile, 0.975, na.rm = TRUE) else rep(NA_real_, length(times))

  out <- data.frame(
    time = times,
    original_prob = original_prob,
    imputed_mean = imputed_mean,
    imputed_sd = imputed_sd,
    imputed_q025 = imputed_q025,
    imputed_q975 = imputed_q975,
    n_imputations = length(x$imputed_datasets),
    stringsAsFactors = FALSE
  )

  # Add formatted columns (keep raw numeric columns)
  out$original <- ifelse(is.na(out$original_prob), "N/A", sprintf("%.1f%%", 100 * out$original_prob))
  out$imputed  <- ifelse(is.na(out$imputed_mean),  "N/A",
                         sprintf("%.1f%% (%.1f-%.1f)", 100 * out$imputed_mean, 100 * out$imputed_q025, 100 * out$imputed_q975))

  # Reorder columns 
  out <- out[, c("time", "original", "imputed", "n_imputations",
                 "original_prob", "imputed_mean", "imputed_sd", "imputed_q025", "imputed_q975")]

  return(out)
}

#' Event probability by clinical time horizons (groups)
#'
#' Computes the event probability for each group at specified
#' time horizons, summarising across completed datasets. This uses the
#' cumulative distribution function.
#'
#' @param x A `bayesian_imputation_groups` object
#' @param times Optional numeric vector of time horizons. If `NULL`, defaults
#'   are chosen using the first group's `time_unit`.
#' @param method Estimation method: "logspline" (default) or "ecdf"
#' @return A data.frame with columns: `group`, `time`, `original_prob`,
#'   `imputed_mean`, `imputed_sd`, `imputed_q025`, `imputed_q975`, `n_imputations`
#' @examples
#' # event_probability_groups(result_groups, times = c(180, 365))
#' @export
event_probability_groups <- function(x, times = NULL, method = c("logspline", "ecdf")) {
  method <- match.arg(method)

  if (!inherits(x, "bayesian_imputation_groups")) {
    stop("x must be a bayesian_imputation_groups object")
  }

  if (is.null(times)) {
    ref <- x$group_results[[x$group_names[1]]]
    time_unit <- resolve_time_unit_label(ref, fallback = "days")
    times <- default_event_probability_times(time_unit)
  }

  rows <- list()
  for (g in x$group_names) {
    res_g <- event_probability(x$group_results[[g]], times = times, method = method)
    res_g$group <- g
    rows[[length(rows) + 1]] <- res_g[, c("group", "time", "original_prob", "imputed_mean", "imputed_sd", "imputed_q025", "imputed_q975", "n_imputations")]
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  # Add formatted columns 
  out$original <- ifelse(is.na(out$original_prob), "N/A", sprintf("%.1f%%", 100 * out$original_prob))
  out$imputed  <- ifelse(is.na(out$imputed_mean),  "N/A",
                         sprintf("%.1f%% (%.1f-%.1f)", 100 * out$imputed_mean, 100 * out$imputed_q025, 100 * out$imputed_q975))

  # Reorder columns to surface the friendly display first
  out <- out[, c("group", "time", "original", "imputed", "n_imputations",
                 "original_prob", "imputed_mean", "imputed_sd", "imputed_q025", "imputed_q975")]

  return(out)
} 
