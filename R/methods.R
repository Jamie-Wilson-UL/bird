# S3 Methods for bayesian_imputation objects

# Default reporting horizons for cumulative event probabilities.
default_event_time_points <- function(time_unit = "days") {
  unit <- tolower(time_unit %||% "days")
  if (unit %in% c("month", "months")) return(c(1, 3, 6, 12))
  if (unit %in% c("year", "years")) return(c(0.25, 0.5, 1, 2))
  c(30, 90, 180, 365)
}

format_time_horizon_label <- function(time_point, time_unit = "days") {
  value <- as.numeric(time_point)
  if (!is.finite(value)) return(as.character(time_point))
  value_label <- if (abs(value - round(value)) < .Machine$double.eps^0.5) {
    format(round(value), trim = TRUE, scientific = FALSE)
  } else {
    format(round(value, 2), trim = TRUE, scientific = FALSE)
  }
  paste(value_label, time_unit %||% "days")
}

default_one_year_equivalent <- function(time_unit = "days") {
  unit <- tolower(time_unit %||% "days")
  if (unit %in% c("month", "months")) return(12)
  if (unit %in% c("year", "years")) return(1)
  365
}

#' Print method for bayesian_imputation objects
#' @param x A bayesian_imputation object
#' @param ... Additional arguments (unused)
#' @export
print.bayesian_imputation <- function(x, ...) {
  # Optional banner (suppressible when nested)
  if (!isTRUE(attr(x, "suppress_banner"))) {
    cat("\n")
    cat("==============================================================\n")
    cat("   BAYESIAN IMPUTATION FOR CENSORED SURVIVAL DATA\n")
    cat("==============================================================\n")
    cat("\n")
  }
  
  # Model information compact format
  cat("MODEL INFORMATION\n")
  cat("------------------\n")
  cat("Distribution:", toupper(x$model_info$distribution), "\n")
  cat("Imputations:", x$model_info$n_imputations, "\n")
  if (!is.null(x$fit_metrics$waic)) {
    cat(sprintf("Model fit (WAIC, lower is better): %.2f  [p_waic=%.2f]", x$fit_metrics$waic$waic, x$fit_metrics$waic$p_waic), "\n")
  }
  cat("Runtime:", format_time(x$model_info$runtime), "\n")
  cat("\n")
  
  # Data summary 
  cat("DATA SUMMARY\n")
  cat("------------\n")
  cat("Total observations:", x$model_info$data_summary$n_total, "\n")
  cat("Observed events:", x$model_info$data_summary$n_observed, "\n")
  cat("Censored:", x$model_info$data_summary$n_censored, 
      paste0("(", round(x$model_info$data_summary$censoring_rate * 100, 1), "%)"), "\n")
  cat("\n")
  
  # Clinical survival metrics 
  if (x$model_info$data_summary$n_observed > 0) {
    clinical_metrics <- calculate_clinical_metrics(x)
    
    cat("CLINICAL SUMMARY\n")
    cat("----------------\n")
    
    # Median survival comparison 
    orig_median <- format_survival_time(clinical_metrics$original$median, x = x)
    imp_median <- format_survival_time(clinical_metrics$imputed$median_summary["mean"], x = x)
    imp_ci_low <- format_survival_time(clinical_metrics$imputed$median_summary["q025"], x = x)
    imp_ci_high <- format_survival_time(clinical_metrics$imputed$median_summary["q975"], x = x)
    
    cat("Median Survival:\n")
    cat("  Original (observed):", orig_median, "\n")
    cat("  After imputation:   ", imp_median, "\n")
    cat("  95% CI:            [", imp_ci_low, "-", imp_ci_high, "]\n")
    
    # Change in survival 
    if (!is.na(clinical_metrics$change$median_change)) {
      sign_chr <- if (clinical_metrics$change$median_change >= 0) "+" else "-"
      cat("  Change:             ", sign_chr, " ", 
          format_survival_time(abs(clinical_metrics$change$median_change), x = x),
          " (", round(abs(clinical_metrics$change$median_change_percent), 1), "%)\n", sep = "")
    }
    cat("\n")
    
    # Cumulative event probabilities at key time points (imputed mean)
    time_unit <- x$model_info$time_unit %||% "days"
    time_points <- default_event_time_points(time_unit)
    ep <- tryCatch(event_probability(x, times = time_points), error = function(e) NULL)
    if (!is.null(ep)) {
      cat("Cumulative event probability (imputed mean):\n")
      for (i in seq_along(time_points)) {
        val <- ep$imputed_mean[i]
        if (is.finite(val)) {
          cat(sprintf("  %s: %.1f%%\n", format_time_horizon_label(time_points[i], time_unit), 100 * val))
        }
      }
      cat("\n")
    }
  }
  
  # Imputation summary 
  if (x$model_info$data_summary$n_censored > 0) {
    imputation_summary <- calculate_imputation_summary(x)
    
    if (imputation_summary$has_imputations) {
      cat("IMPUTATION SUMMARY\n")
      cat("------------------\n")
      cat("Censored observations:", imputation_summary$n_censored, "\n")
      
      orig_med <- format_survival_time(imputation_summary$original_censored_summary["median"], x = x)
      imp_med <- format_survival_time(imputation_summary$imputed_times_summary["median"], x = x)
      gain_med <- format_survival_time(imputation_summary$imputation_gains_summary["median"], x = x)
      
      cat("Original censoring (median):", orig_med, "\n")
      cat("Imputed times (median):     ", imp_med, "\n")
      cat("Median extension beyond censoring:", gain_med, 
          " (imputed ", imp_med, " vs censoring ", orig_med, ")\n", sep = "")
      
      # Quality indicator
      if (imputation_summary$all_gains_positive) {
        cat("Quality: [OK] All imputations valid\n")
      } else {
        cat("Quality: [WARNING] Some imputations may be invalid\n")
      }
      cat("\n")
    }
  }
  
  # MCMC diagnostics 
  cat("MCMC DIAGNOSTICS\n")
  cat("----------------\n")
  cat("Chains:", x$model_info$mcmc_options$chains, "\n")
  cat("Iterations:", x$model_info$mcmc_options$iter_warmup, "warmup +", 
      x$model_info$mcmc_options$iter_sampling, "sampling\n")
  
  # Diagnostics for parametric models
  if (!identical(x$model_info$distribution, "nonparametric_lddp")) {
    dx <- x$diagnostics
    if (!is.null(dx)) {
      # Print if values are available
      if (!is.null(dx$max_rhat) || !is.null(dx$min_ess_bulk) || !is.null(dx$min_ess_tail)) {
        cat(sprintf("Max Rhat: %s | Min ESS(bulk): %s | Min ESS(tail): %s\n",
          if (!is.null(dx$max_rhat)) sprintf("%.3f", as.numeric(dx$max_rhat)) else "NA",
          if (!is.null(dx$min_ess_bulk)) sprintf("%.0f", as.numeric(dx$min_ess_bulk)) else "NA",
          if (!is.null(dx$min_ess_tail)) sprintf("%.0f", as.numeric(dx$min_ess_tail)) else "NA"
        ))
      }
      if (!is.null(dx$num_divergent) || !is.null(dx$num_max_treedepth)) {
        cat(sprintf("Divergences: %s | Max treedepth hits: %s\n",
          if (!is.null(dx$num_divergent)) sprintf("%d", as.integer(dx$num_divergent)) else "NA",
          if (!is.null(dx$num_max_treedepth)) sprintf("%d", as.integer(dx$num_max_treedepth)) else "NA"
        ))
      }
    }
  }
  
  # Convergence status with indicator
  if (x$diagnostics$convergence_ok) {
    cat("Convergence: [OK] Good\n")
  } else {
    cat("Convergence: [WARNING] Check diagnostics\n")
  }
  cat("\n")
  
  # Next steps shown only when not suppressed
  if (!isTRUE(attr(x, "suppress_next_steps"))) {
    cat("NEXT STEPS\n")
    cat("----------\n")
    cat("- plot(result)                    # View survival curves\n")
    cat("- complete(result, dataset = 1)   # Extract completed dataset\n")
    if (x$model_info$n_imputations <= 10) {
      cat("- plot(result, type = 'boxplots_comparison')  # Compare datasets\n")
      cat("- plot(result, type = 'density')              # Density comparison\n")
    } else {
      cat("- plot(result, type = 'boxplots_comparison')  # Compare first 10 datasets (", x$model_info$n_imputations, " total)\n", sep = "")
      cat("- plot(result, type = 'density')              # Density comparison\n")
    }
    cat("\n")
  }
  
  invisible(x)
}

#' Plot method for bayesian_imputation objects
#' @param x A bayesian_imputation object
#' @param type Type of plot: "survival", "trace", "pairs", "posterior", "completed_dataset_summary", "boxplots_comparison", "density"
#' @param n_curves For survival plots: number of imputed curves to display (default: 10)
#' @param alpha For survival plots: transparency for imputed curves (default: 0.3)
#' @param show_original For survival plots: whether to show original Kaplan-Meier curve (default: TRUE)
#' @param n_max Backward-compatible alias:
#'   for `type = "survival"`, treated as `n_curves`;
#'   for `type = "boxplots_comparison"`, maximum datasets shown.
#' @param ... Additional arguments passed to plotting functions.
#'   For `type = "completed_dataset_summary"` and `type = "density"`,
#'   you can pass `dataset_id = <integer>` to select a specific imputed dataset
#'   (random if omitted).
#'   For `type = "completed_dataset_summary"`, you can also pass
#'   `panels = c("hist", "density", "survival", "boxplot")`
#'   (or `"auto"` for the default layout).
#' @export
plot.bayesian_imputation <- function(x, type = "survival", n_curves = 10, alpha = 0.3, show_original = TRUE, n_max = NULL, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  if (!is.null(n_max) && type == "survival") {
    n_curves <- n_max
  }
  
  switch(type,
    "survival" = plot_survival_curves(x, n_curves = n_curves, alpha = alpha, show_original = show_original, ...),
    "trace" = plot_trace_plots(x, ...),
    "pairs" = plot_pairs(x, ...),
    "posterior" = plot_posterior_censored(x, ...),
    "completed_dataset_summary" = plot_completed_dataset_summary(x, ...),
    "boxplots_comparison" = plot_boxplots_comparison(x, n_max = n_max %||% 10, ...),
    "density" = plot_density_comparison(x, ...),
    stop("Unknown plot type: ", type, ". Options: 'survival', 'trace', 'pairs', 'posterior', 'completed_dataset_summary', 'boxplots_comparison', 'density'")
  )
}

# S3 Methods for bayesian_imputation_groups objects

#' Print method for bayesian_imputation_groups objects
#' @param x A bayesian_imputation_groups object
#' @param ... Additional arguments (unused)
#' @export
print.bayesian_imputation_groups <- function(x, ...) {
  # Top-level banner printed once for group output
  if (!isTRUE(attr(x, "suppress_banner"))) {
    cat("\n")
    cat("==============================================================\n")
    cat("   BAYESIAN IMPUTATION GROUP SUMMARY\n")
    cat("==============================================================\n")
    cat("\n")
  }
  
  # Combined overview header and context
  cat("GROUP OVERVIEW:\n")
  cat("--------------\n")
  cat("Groups:", paste(x$group_names, collapse = ", "), "\n")
  cat("Group variable:", x$groups, "\n")
  cat("Total groups analyzed:", length(x$group_names), "\n")
  
  # Show any errors first
  if (length(x$group_errors) > 0) {
    cat("\nWARNINGS:\n")
    cat("---------\n")
    for (group_name in names(x$group_errors)) {
      cat("  Group '", group_name, "': ", x$group_errors[[group_name]], "\n", sep = "")
    }
  }
  
  # Group lines
  if (length(x$group_names) > 0) {
    print_group_overview_safe(x)
  }
  
  # Detailed results for each group (suppress banner and NEXT STEPS in nested prints)
  if (length(x$group_names) > 0) {
    cat("\nDETAILED GROUP RESULTS:\n")
    cat("---------------------\n")
    for (group_name in x$group_names) {
      cat("\nGROUP:", group_name, "\n")
      cat("----------------------------------------------\n")
      gr <- x$group_results[[group_name]]
      attr(gr, "suppress_banner") <- TRUE
      attr(gr, "suppress_next_steps") <- TRUE
      print(gr)
    }
  }
  
  # Add comprehensive comparison if we have at least 2 successful groups
  if (length(x$group_names) >= 2) {
    cat("\nGROUP COMPARISONS\n")
    cat("-----------------\n")
    print_comprehensive_group_comparison_safe(x)
  }
}

#' Plot method for bayesian_imputation_groups objects
#' @param x A bayesian_imputation_groups object
#' @param type Type of plot ("survival", "trace", "pairs", "posterior", "completed_dataset_summary", "boxplots_comparison", "density")
#' @param combine_groups For survival plots: whether to show all groups on same plot (TRUE) or separate plots (FALSE)
#' @param n_curves For survival plots: number of imputed curves to show per group (default: 10)
#' @param alpha For survival plots: transparency for imputed curves
#' @param show_original For survival plots: whether to show original Kaplan-Meier curves
#' @param n_max Backward-compatible alias:
#'   for `type = "survival"`, treated as `n_curves`;
#'   for `type = "boxplots_comparison"`, maximum datasets shown.
#' @param ... Additional arguments passed to plotting functions.
#'   For `type = "completed_dataset_summary"`, you can pass:
#'   `dataset_id = <integer>` and
#'   `panels = c("hist", "density", "survival", "boxplot")`
#'   (or `"auto"` for the default layout).
#' @export
plot.bayesian_imputation_groups <- function(x, type = "survival", combine_groups = TRUE, 
                                           n_curves = 10, alpha = 0.3, show_original = TRUE, n_max = NULL, ...) {
  
  # Check if we have enough groups to plot
  if (length(x$group_names) < 1) {
    stop("No successful groups to plot")
  }

  if (!is.null(n_max) && type == "survival") {
    n_curves <- n_max
  }
  
  switch(type,
    "survival" = plot_survival_curves_by_group_safe(x, n_curves, alpha, show_original, combine_groups, ...),
    "trace" = plot_individual_groups_safe(x, "trace", ...),
    "pairs" = plot_individual_groups_safe(x, "pairs", ...),
    "posterior" = plot_posterior_censored_groups(x, ...),
    "completed_dataset_summary" = plot_completed_dataset_summary_groups(x, ...),
    "boxplots_comparison" = plot_boxplots_comparison_groups(x, n_max = n_max %||% 10, ...),
    "density" = plot_density_comparison_groups(x, ...),
    stop("Unknown plot type: ", type, ". Options: 'survival', 'trace', 'pairs', 'posterior', 'completed_dataset_summary', 'boxplots_comparison', 'density'")
  )
}

#' Complete method for bayesian_imputation_groups objects
#' @param object A bayesian_imputation_groups object
#' @param dataset Dataset number to extract (NULL for all)
#' @param format Output format ("wide", "long", "list")
#' @param groups How to handle groups ("combined", "separate", or specific group name)
#' @param ... Additional arguments
#' @export
complete.bayesian_imputation_groups <- function(object, dataset = NULL, 
                                              format = c("wide", "long", "list"),
                                              groups = "combined", ...) {
  
  format <- match.arg(format)
  
  if (groups == "combined") {
    # Combine all groups into single datasets (default)
    combined_datasets <- combine_group_datasets_safe(object, dataset, format)
    return(combined_datasets)
    
  } else if (groups == "separate") {
    # Return separate datasets for each group
    separate_datasets <- list()
    for (group_name in object$group_names) {
      separate_datasets[[group_name]] <- complete(
        object$group_results[[group_name]], dataset, format
      )
    }
    return(separate_datasets)
    
  } else {
    # Return specific group
    if (groups %in% object$group_names) {
      return(complete(object$group_results[[groups]], dataset, format))
    } else {
      stop("Group '", groups, "' not found. Available groups: ", 
           paste(object$group_names, collapse = ", "))
    }
  }
}

# Helper functions for group methods

#' Print group overview with key metrics
#' @param x bayesian_imputation_groups object
#' @keywords internal
print_group_overview_safe <- function(x) {
  # Create a summary table of key metrics for each group
  overview_data <- data.frame(
    Group = character(),
    N_Total = integer(),
    N_Censored = integer(),
    Censoring_Rate = numeric(),
    Median_Mean_Days = numeric(),
    Median_Low_Days = numeric(),
    Median_High_Days = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (group_name in x$group_names) {
    group_result <- x$group_results[[group_name]]
    
    # Get clinical metrics
    clinical_metrics <- calculate_clinical_metrics(group_result)
    median_mean_num <- as.numeric(clinical_metrics$imputed$median_summary["mean"])
    median_low_num  <- as.numeric(clinical_metrics$imputed$median_summary["q25"])
    median_high_num <- as.numeric(clinical_metrics$imputed$median_summary["q75"])
    
    # Add to overview data
    overview_data <- rbind(overview_data, data.frame(
      Group = group_name,
      N_Total = group_result$model_info$data_summary$n_total,
      N_Censored = group_result$model_info$data_summary$n_censored,
      Censoring_Rate = round(group_result$model_info$data_summary$censoring_rate * 100, 1),
      Median_Mean_Days = median_mean_num,
      Median_Low_Days = median_low_num,
      Median_High_Days = median_high_num,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print overview lines 
  for (i in seq_len(nrow(overview_data))) {
    row <- overview_data[i, ]
    cat(sprintf(
      "  %s: %d obs (%d censored, %.1f%%), median survival: %s [%s - %s]\n",
      row$Group,
      row$N_Total,
      row$N_Censored,
      row$Censoring_Rate,
      format_survival_time(row$Median_Mean_Days, x = x$group_results[[row$Group]]),
      format_survival_time(row$Median_Low_Days, x = x$group_results[[row$Group]]),
      format_survival_time(row$Median_High_Days, x = x$group_results[[row$Group]])
    ))
  }
}


#' Print survival comparison between groups
#' @param x bayesian_imputation_groups object
#' @param group_metrics List of clinical metrics for each group
#' @keywords internal
print_survival_comparison_safe <- function(x, group_metrics) {
  # Compare median survival times
  medians <- sapply(group_metrics, function(gm) gm$imputed$median_summary["mean"])
  names(medians) <- x$group_names
  
  # Find best and worst performing groups
  best_group <- names(medians)[which.max(medians)]
  worst_group <- names(medians)[which.min(medians)]
  
  cat("Median survival times:\n")
  for (group_name in x$group_names) {
    median_time <- format_survival_time(medians[group_name], x = x$group_results[[group_name]])
    cat(sprintf("  %s: %s\n", group_name, median_time))
  }
  
  # Calculate difference between best and worst
  if (length(medians) >= 2) {
    best_median <- medians[best_group]
    worst_median <- medians[worst_group]
    difference <- best_median - worst_median
    
    cat(sprintf("\nSurvival difference (best - worst): %s\n", 
                format_survival_time(difference, x = x$group_results[[best_group]])))
    cat(sprintf("Best performing group: %s\n", best_group))
    cat(sprintf("Worst performing group: %s\n", worst_group))
  }
  
  # Event probabilities at key time points (from imputed datasets)
  time_unit <- x$group_results[[1]]$model_info$time_unit %||% "days"
  time_points <- default_event_time_points(time_unit)
  has_any <- FALSE
  rows <- list()
  for (tp in time_points) {
    row_vals <- c()
    for (group_name in x$group_names) {
      gr <- x$group_results[[group_name]]
      # Use the user-facing helper to compute event probabilities across imputations
      ep <- tryCatch(event_probability(gr, times = tp), error = function(e) NULL)
      val <- if (!is.null(ep)) ep$imputed_mean[1] else NA_real_
      row_vals[group_name] <- val
      if (is.finite(val)) has_any <- TRUE
    }
    rows[[length(rows)+1]] <- list(time = tp, vals = row_vals)
  }
  if (has_any) {
    cat("\nCumulative event probability at key time points (imputed mean):\n")
    for (r in rows) {
      cat(sprintf("  %s:\n", format_time_horizon_label(r$time, time_unit)))
      for (group_name in x$group_names) {
        val <- r$vals[group_name]
        if (is.finite(val)) {
          cat(sprintf("    %s: %.1f%%\n", group_name, 100*val))
        }
      }
    }
  }
}


#' Calculate pooled Cox Wald test (log-rank equivalent) for group comparison
#' @param object bayesian_imputation_groups object
#' @keywords internal
calculate_log_rank_test <- function(object) {
  tryCatch({
    pooled <- pool_cox_group_test(object)
    if (is.null(pooled)) return(NULL)
    list(
      p_value = pooled$p_value,
      chi_square = pooled$wald,
      df = pooled$df
    )
  }, error = function(e) {
    return(NULL)
  })
}

#' Calculate pooled Cox hazard-ratio summary for group comparison
#' @param object bayesian_imputation_groups object
#' @keywords internal
calculate_cox_ph_test <- function(object) {
  tryCatch({
    pooled <- pool_cox_group_test(object)
    if (is.null(pooled)) return(NULL)
    
    # Surface if fewer imputations were usable
    total_m <- length(object$group_results[[1]]$imputed_datasets)
    if (!is.null(pooled$m_used) && pooled$m_used < total_m) {
      cat(sprintf("Note: pooled Cox used %d/%d imputations due to singular fits.\n", pooled$m_used, total_m))
    }

    # If there is only one contrast (two groups), provide HR and CI
    hr <- NA_real_
    hr_ci_lower <- NA_real_
    hr_ci_upper <- NA_real_
    if (length(pooled$coef) == 1 && is.matrix(pooled$vcov)) {
      se <- sqrt(diag(pooled$vcov))[1]
      hr <- exp(pooled$coef[1])
      hr_ci_lower <- exp(pooled$coef[1] - 1.96 * se)
      hr_ci_upper <- exp(pooled$coef[1] + 1.96 * se)
    }
    
    list(
      p_value = pooled$p_value,
      hazard_ratio = hr,
      hr_ci_lower = hr_ci_lower,
      hr_ci_upper = hr_ci_upper
    )
  }, error = function(e) {
    return(NULL)
  })
}

#' Pool Cox PH group comparison across imputations (Rubin's rules)
#' @param object bayesian_imputation_groups object
#' @keywords internal
pool_cox_group_test <- function(object) {
  tryCatch({
    m_total <- length(object$group_results[[1]]$imputed_datasets)
    
    coefs <- list()
    covs <- list()
    used <- 0L
    
    for (i in seq_len(m_total)) {
      dat <- data.frame()
      for (g in object$group_names) {
        d <- object$group_results[[g]]$imputed_datasets[[i]]
        d$group <- g
        dat <- rbind(dat, d)
      }
      dat$group <- factor(dat$group, levels = object$group_names)
      
      fit <- tryCatch({
        survival::coxph(
          survival::Surv(
            dat[[object$group_results[[1]]$time_col]],
            dat[[object$group_results[[1]]$status_col]]
          ) ~ group,
          data = dat
        )
      }, error = function(e) NULL)
      
      if (!is.null(fit)) {
        cf <- tryCatch(stats::coef(fit), error = function(e) NULL)
        vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
        if (!is.null(cf) && !is.null(vc) && all(is.finite(cf)) && all(is.finite(vc))) {
          used <- used + 1L
          coefs[[used]] <- cf
          covs[[used]] <- vc
        }
      }
    }
    
    if (used == 0L) return(NULL)
    
    # Align coefficient vectors in case of dropped levels
    ref_levels <- names(coefs[[1]])
    for (j in seq_len(used)) {
      if (!identical(names(coefs[[j]]), ref_levels)) {
        tmp <- rep(NA_real_, length(ref_levels))
        names(tmp) <- ref_levels
        tmp[names(coefs[[j]])] <- coefs[[j]]
        coefs[[j]] <- tmp
        # Expand covariance accordingly
        Vtmp <- matrix(NA_real_, nrow = length(ref_levels), ncol = length(ref_levels))
        rownames(Vtmp) <- colnames(Vtmp) <- ref_levels
        Vsrc <- covs[[j]]
        Vtmp[rownames(Vsrc), colnames(Vsrc)] <- Vsrc
        covs[[j]] <- Vtmp
      }
    }
    
    coef_mat <- do.call(rbind, coefs)
    p <- ncol(coef_mat)
    if (p == 0) return(NULL)
    
    # Within- and between-imputation variance (stabilise with tiny ridge)
    eps <- 1e-10
    base_dim <- ncol(covs[[1]])
    W <- Reduce("+", covs) / used + diag(eps, base_dim)
    if (used >= 2) {
      B <- stats::cov(coef_mat)
    } else {
      B <- matrix(0, nrow = p, ncol = p)
      rownames(B) <- colnames(B) <- colnames(coef_mat)
    }
    bbar <- colMeans(coef_mat)
    Tmat <- W + (1 + 1/used) * B
    
    # Stabilize in case of near-singularity
    solve_T <- tryCatch(solve(Tmat), error = function(e) NULL)
    if (is.null(solve_T)) {
      ridge <- 1e-8
      Tmat <- Tmat + diag(ridge, nrow(Tmat))
      solve_T <- tryCatch(solve(Tmat), error = function(e) NULL)
    }
    if (is.null(solve_T)) return(NULL)
    
    wald <- as.numeric(t(bbar) %*% solve_T %*% bbar)
    df <- length(bbar)
    p_value <- stats::pchisq(wald, df = df, lower.tail = FALSE)
    
    list(coef = bbar, vcov = Tmat, wald = wald, df = df, p_value = p_value, m_used = used)
  }, error = function(e) {
    NULL
  })
}

#' Calculate median survival difference between groups
#' @param object bayesian_imputation_groups object
#' @keywords internal
calculate_median_survival_diff <- function(object) {
  tryCatch({
    medians <- numeric(length(object$group_names))
    names(medians) <- object$group_names
    for (i in seq_along(object$group_names)) {
      group_name <- object$group_names[i]
      group_result <- object$group_results[[group_name]]
      clinical_metrics <- calculate_clinical_metrics(group_result)
      medians[i] <- clinical_metrics$imputed$median_summary["mean"]
    }
    best_group <- names(medians)[which.max(medians)]
    worst_group <- names(medians)[which.min(medians)]
    difference <- medians[best_group] - medians[worst_group]
    list(difference = difference, best_group = best_group, worst_group = worst_group)
  }, error = function(e) {
    NULL
  })
}

#' Calculate survival probability difference at specific time point
#' @param object bayesian_imputation_groups object
#' @param time_point Time point in the same units as the input data
#' @keywords internal
calculate_survival_probability_diff <- function(object, time_point) {
  tryCatch({
    probs <- numeric(length(object$group_names))
    names(probs) <- object$group_names
    for (i in seq_along(object$group_names)) {
      gname <- object$group_names[i]
      gr <- object$group_results[[gname]]
      ep <- event_probability(gr, times = time_point)
      if (nrow(ep) >= 1 && is.finite(ep$imputed_mean[1])) {
        probs[i] <- 1 - ep$imputed_mean[1]
      } else {
        probs[i] <- NA_real_
      }
    }
    if (all(is.finite(probs))) {
      best_group <- names(probs)[which.max(probs)]
      worst_group <- names(probs)[which.min(probs)]
      difference <- probs[best_group] - probs[worst_group]
      list(difference = difference, best_group = best_group, worst_group = worst_group, time_point = time_point)
    } else {
      NULL
    }
  }, error = function(e) {
    NULL
  })
}

#' Print comprehensive group comparison (with statistical tests)
#' @param object bayesian_imputation_groups object
#' @keywords internal
print_comprehensive_group_comparison_safe <- function(object) {
  if (length(object$group_names) < 2) {
    cat("Need at least 2 groups for comparison\n")
    return()
  }
  
  # Collect clinical metrics for all groups
  group_metrics <- list()
  for (group_name in object$group_names) {
    group_result <- object$group_results[[group_name]]
    group_metrics[[group_name]] <- calculate_clinical_metrics(group_result)
  }
  
  # Show survival comparison
  cat("SURVIVAL COMPARISON:\n")
  cat("-------------------\n")
  print_survival_comparison_safe(object, group_metrics)
  
  # Statistical tests
  cat("\n1. Global Group Effect (Pooled Cox Wald Test):\n")
  cat("----------------------------------------------\n")
  log_rank_results <- calculate_log_rank_test(object)
  if (!is.null(log_rank_results$p_value)) {
    cat(sprintf("  Pooled Cox (Wald) P-value: %.4f\n", log_rank_results$p_value))
    if (log_rank_results$p_value < 0.05) {
      cat("  (Statistically significant difference in survival)\n")
    } else {
      cat("  (No statistically significant difference in survival)\n")
    }
  } else {
    cat("  Pooled Cox Wald test could not be calculated.\n")
  }
  cat("\n2. Hazard Ratio Summary (from pooled Cox model):\n")
  cat("------------------------------------------------\n")
  cox_results <- calculate_cox_ph_test(object)
  if (!is.null(cox_results)) {
    if (!is.na(cox_results$hazard_ratio)) {
      cat(sprintf("  Hazard ratio: %.3f\n", cox_results$hazard_ratio))
      if (!is.na(cox_results$hr_ci_lower) && !is.na(cox_results$hr_ci_upper)) {
        cat(sprintf("  95%% CI: [%.3f, %.3f]\n", cox_results$hr_ci_lower, cox_results$hr_ci_upper))
      }
    } else {
      cat("  Multiple groups detected: a single HR is not shown in this summary.\n")
    }
  } else {
    cat("  Hazard-ratio summary could not be calculated.\n")
  }
  cat("\n3. Median Survival Difference:\n")
  cat("------------------------------\n")
  median_diff_results <- calculate_median_survival_diff(object)
  if (!is.null(median_diff_results$difference)) {
    unit_source <- object$group_results[[median_diff_results$best_group]]
    cat(sprintf("  Median Survival Difference: %s\n", format_survival_time(median_diff_results$difference, x = unit_source)))
    cat(sprintf("  Best group: %s | Worst group: %s\n", median_diff_results$best_group, median_diff_results$worst_group))
  } else {
    cat("  Median Survival Difference could not be calculated.\n")
  }
  time_unit <- object$group_results[[1]]$model_info$time_unit %||% "days"
  one_year_equiv <- default_one_year_equivalent(time_unit)
  cat("\n4. Overall Survival Probability at 1-Year Equivalent:\n")
  cat("--------------------------------------\n")
  survival_prob_results <- calculate_survival_probability_diff(object, one_year_equiv)
  if (!is.null(survival_prob_results$difference)) {
    cat(sprintf(
      "  Survival Probability Difference at %s: %.1f%%\n",
      format_time_horizon_label(one_year_equiv, time_unit),
      survival_prob_results$difference * 100
    ))
    cat(sprintf("  Best group: %s | Worst group: %s\n", survival_prob_results$best_group, survival_prob_results$worst_group))
  } else {
    cat(sprintf(
      "  Survival Probability Difference at %s could not be calculated.\n",
      format_time_horizon_label(one_year_equiv, time_unit)
    ))
  }
}

# Group plotting functions moved to R/plotting-helpers.R

#' Combine group datasets safely
#' @param object bayesian_imputation_groups object
#' @param dataset Dataset number
#' @param format Output format
#' @keywords internal
combine_group_datasets_safe <- function(object, dataset, format) {
  # Combine datasets from all groups
  time_col <- object$group_results[[1]]$time_col
  status_col <- object$group_results[[1]]$status_col
  
  if (format == "list") {
    # Return list of combined datasets
    combined_list <- list()
    n_datasets <- length(object$group_results[[1]]$imputed_datasets)
    
    for (i in seq_len(n_datasets)) {
      combined_data <- data.frame()
      for (group_name in object$group_names) {
        group_data <- object$group_results[[group_name]]$imputed_datasets[[i]]
        combined_data <- rbind(combined_data, group_data)
      }
      combined_data <- format_completed_dataset_output(
        combined_data,
        time_col = time_col
      )
      combined_data <- standardize_complete_column_order(
        combined_data,
        time_col = "imputed_time",
        status_col = status_col
      )
      combined_list[[i]] <- combined_data
    }
    return(combined_list)
    
  } else if (format == "wide") {
    # Return single combined dataset
    if (is.null(dataset)) {
      # Return first dataset by default
      dataset <- 1
    }
    
    combined_data <- data.frame()
    for (group_name in object$group_names) {
      group_data <- object$group_results[[group_name]]$imputed_datasets[[dataset]]
      combined_data <- rbind(combined_data, group_data)
    }
    combined_data <- format_completed_dataset_output(
      combined_data,
      time_col = time_col
    )
    combined_data <- standardize_complete_column_order(
      combined_data,
      time_col = "imputed_time",
      status_col = status_col
    )
    return(combined_data)
    
  } else if (format == "long") {
    # Return long format with all datasets combined
    combined_data <- data.frame()
    n_datasets <- length(object$group_results[[1]]$imputed_datasets)
    
    for (i in seq_len(n_datasets)) {
      for (group_name in object$group_names) {
        group_data <- object$group_results[[group_name]]$imputed_datasets[[i]]
        group_data$.imp <- i
        combined_data <- rbind(combined_data, group_data)
      }
    }
    
    # Add row IDs
    combined_data$.id <- seq_len(nrow(combined_data))
    combined_data <- format_completed_dataset_output(
      combined_data,
      time_col = time_col
    )
    combined_data <- standardize_complete_column_order(
      combined_data,
      time_col = "imputed_time",
      status_col = status_col
    )
    return(combined_data)
  }
}
