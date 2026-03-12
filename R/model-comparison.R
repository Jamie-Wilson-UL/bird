#' Compare Multiple Imputation Models
#'
#' Fit multiple imputation models on the same right-censored dataset and collect
#' outputs in a single comparison object.
#'
#' Supported model names are:
#' - Parametric: `"weibull"`, `"exponential"`, `"lognormal"`
#' - Nonparametric aliases: `"nonparametric"`, `"np"`, `"lddp"`, `"nonparametric_lddp"`
#'
#' @param data Data frame or `survival_data` object.
#' @param models Character vector of model names to compare.
#' @param time Optional time column name (auto-detected if missing).
#' @param status Optional status column name (auto-detected if missing).
#' @param n_imputations Number of completed datasets per model.
#' @param time_unit Time unit label (e.g., `"days"`).
#' @param priors Parametric prior list (passed to parametric models).
#' @param prior Nonparametric prior list (passed to LDDP model).
#' @param mcmc_options Parametric MCMC options.
#' @param mcmc Nonparametric MCMC options.
#' @param verbose Logical; print progress.
#' @param ... Additional arguments passed to [impute()].
#'
#' @return An object of class `bird_model_comparison`.
#' @export
compare_models <- function(data,
                           models = c("weibull", "nonparametric"),
                           time = NULL,
                           status = NULL,
                           n_imputations = 10,
                           time_unit = "days",
                           priors = NULL,
                           prior = NULL,
                           mcmc_options = NULL,
                           mcmc = NULL,
                           verbose = TRUE,
                           ...) {
  if (!is.character(models) || length(models) < 2) {
    stop("'models' must be a character vector with at least two entries")
  }

  normalize_model <- function(x) {
    token <- tolower(trimws(x))
    if (token %in% c("weibull", "exponential", "lognormal")) return(token)
    if (token %in% c("nonparametric", "np", "lddp", "nonparametric_lddp")) return("nonparametric")
    stop(
      "Unknown model '", x, "'. Allowed values are: ",
      "'weibull', 'exponential', 'lognormal', 'nonparametric' (or aliases np/lddp/nonparametric_lddp)."
    )
  }

  model_tokens <- unique(vapply(models, normalize_model, character(1)))
  if (length(model_tokens) < 2) {
    stop("Need at least two distinct models to compare")
  }

  prepared <- data
  if (!inherits(prepared, "survival_data")) {
    prepared <- prepare_survival_data(
      data = prepared,
      time = time,
      status = status,
      time_unit = time_unit,
      verbose = FALSE
    )
  }

  fit_start <- Sys.time()
  model_results <- list()
  model_errors <- list()

  for (model_name in model_tokens) {
    if (verbose) {
      cat("Fitting model:", model_name, "\n")
    }

    fit <- tryCatch(
      {
        if (identical(model_name, "nonparametric")) {
          impute(
            data = prepared,
            model = "nonparametric",
            distribution = "nonparametric",
            n_imputations = n_imputations,
            time_unit = time_unit,
            prior = prior,
            mcmc = mcmc,
            verbose = verbose,
            ...
          )
        } else {
          impute(
            data = prepared,
            model = "parametric",
            distribution = model_name,
            n_imputations = n_imputations,
            time_unit = time_unit,
            priors = priors,
            mcmc_options = mcmc_options,
            verbose = verbose,
            ...
          )
        }
      },
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      model_errors[[model_name]] <- conditionMessage(fit)
      if (verbose) {
        cat("  - Failed:", conditionMessage(fit), "\n")
      }
    } else {
      model_results[[model_name]] <- fit
    }
  }

  if (length(model_results) < 2) {
    stop(
      "Need at least two successful model fits for comparison. ",
      "Successful: ", length(model_results), "; failed: ", length(model_errors), "."
    )
  }

  first <- model_results[[1]]
  fit_runtime <- as.numeric(difftime(Sys.time(), fit_start, units = "secs"))

  out <- list(
    requested_models = model_tokens,
    model_names = names(model_results),
    model_results = model_results,
    model_errors = model_errors,
    n_models = length(model_results),
    n_imputations = n_imputations,
    original_data = first$original_data,
    time_col = first$time_col,
    status_col = first$status_col,
    model_info = list(
      runtime = fit_runtime,
      time_unit = first$model_info$time_unit,
      data_summary = first$model_info$data_summary
    )
  )
  class(out) <- "bird_model_comparison"
  out
}

#' Print method for model comparison objects
#' @param x A `bird_model_comparison` object.
#' @param ... Unused.
#' @export
print.bird_model_comparison <- function(x, ...) {
  cat("\n")
  cat("==============================================================\n")
  cat("   BAYESIAN IMPUTATION MODEL COMPARISON\n")
  cat("==============================================================\n")
  cat("\n")

  cat("Models compared:", paste(x$model_names, collapse = ", "), "\n")
  cat("Total runtime:", format_time(x$model_info$runtime), "\n")
  cat("Imputations per model:", x$n_imputations, "\n")
  cat("\n")

  if (length(x$model_errors) > 0) {
    cat("FAILED MODELS\n")
    cat("-------------\n")
    for (nm in names(x$model_errors)) {
      cat("  -", nm, ":", x$model_errors[[nm]], "\n")
    }
    cat("\n")
  }

  cat("MODEL SUMMARY\n")
  cat("-------------\n")
  for (nm in x$model_names) {
    fit <- x$model_results[[nm]]
    waic_text <- "WAIC: NA"
    if (!is.null(fit$fit_metrics$waic) && !is.null(fit$fit_metrics$waic$waic)) {
      waic_text <- sprintf("WAIC: %.2f", as.numeric(fit$fit_metrics$waic$waic))
    }
    conv_text <- if (isTRUE(fit$diagnostics$convergence_ok)) "convergence: OK" else "convergence: check"
    cat(
      "  -", nm, "|", waic_text, "|", conv_text,
      "| runtime:", format_time(fit$model_info$runtime), "\n"
    )
  }
  cat("\n")

  cat("CLINICAL COMPARISON\n")
  cat("-------------------\n")
  time_unit <- x$model_info$time_unit %||% "days"
  horizon <- default_one_year_equivalent(time_unit)
  horizon_label <- format_time_horizon_label(horizon, time_unit)

  for (nm in x$model_names) {
    fit <- x$model_results[[nm]]
    clinical <- tryCatch(calculate_clinical_metrics(fit), error = function(e) NULL)
    ep <- tryCatch(event_probability(fit, times = horizon), error = function(e) NULL)

    if (is.null(clinical)) {
      cat("  -", nm, "| clinical summary: unavailable\n")
      next
    }

    med_mean <- as.numeric(clinical$imputed$median_summary["mean"])
    med_low <- as.numeric(clinical$imputed$median_summary["q025"])
    med_high <- as.numeric(clinical$imputed$median_summary["q975"])
    med_change <- as.numeric(clinical$change$median_change)

    cat(
      "  -", nm,
      "| median:", format_survival_time(med_mean, x = fit),
      "[", format_survival_time(med_low, x = fit), ",", format_survival_time(med_high, x = fit), "]"
    )
    if (is.finite(med_change)) {
      sign_chr <- if (med_change >= 0) "+" else "-"
      cat("| change vs observed:", paste0(sign_chr, " ", format_survival_time(abs(med_change), x = fit)))
    }
    cat("\n")

    if (!is.null(ep) && nrow(ep) >= 1 && is.finite(ep$imputed_mean[1])) {
      cat(sprintf("      event probability by %s: %.1f%%\n", horizon_label, 100 * ep$imputed_mean[1]))
    } else {
      cat(sprintf("      event probability by %s: unavailable\n", horizon_label))
    }
  }
  cat("\n")

  cat("NEXT STEPS\n")
  cat("----------\n")
  cat("- plot(result)  # Survival comparison across models\n")
  cat("- plot(result, type = 'completed_dataset_summary')\n")
  cat("- complete(result, models = 'combined', dataset = 1)\n")
  cat("\n")

  invisible(x)
}

#' Plot method for model comparison objects
#' @param x A `bird_model_comparison` object.
#' @param type Plot type: `"survival"` or `"completed_dataset_summary"`.
#' @param ... Additional arguments passed to plotting helpers.
#' @export
plot.bird_model_comparison <- function(x, type = c("survival", "completed_dataset_summary"), ...) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  switch(
    type,
    "survival" = plot_survival_curves_models(x, ...),
    "completed_dataset_summary" = plot_completed_dataset_summary_models(x, ...),
    stop("Unknown plot type: ", type)
  )
}

#' Complete method for model comparison objects
#' @param object A `bird_model_comparison` object.
#' @param dataset Dataset number to extract. `NULL` means all (except for wide, where first is used).
#' @param format Output format (`"wide"`, `"long"`, `"list"`).
#' @param models Either `"combined"`, `"separate"`, or a specific model name.
#' @param ... Additional arguments passed to [complete()].
#' @export
complete.bird_model_comparison <- function(object,
                                           dataset = NULL,
                                           format = c("wide", "long", "list"),
                                           models = "combined",
                                           ...) {
  format <- match.arg(format)

  if (identical(models, "combined")) {
    return(combine_model_datasets_safe(object, dataset = dataset, format = format))
  }

  if (identical(models, "separate")) {
    out <- list()
    for (nm in object$model_names) {
      out[[nm]] <- complete(
        object$model_results[[nm]],
        dataset = if (is.null(dataset)) "all" else dataset,
        format = format,
        ...
      )
    }
    return(out)
  }

  if (models %in% object$model_names) {
    return(complete(
      object$model_results[[models]],
      dataset = if (is.null(dataset)) "all" else dataset,
      format = format,
      ...
    ))
  }

  stop(
    "Unknown 'models' value: ", models,
    ". Use 'combined', 'separate', or one of: ",
    paste(object$model_names, collapse = ", ")
  )
}

#' Combine completed datasets across compared models
#' @param object `bird_model_comparison` object.
#' @param dataset Dataset index (optional).
#' @param format One of `"wide"`, `"long"`, `"list"`.
#' @keywords internal
combine_model_datasets_safe <- function(object, dataset = NULL, format = c("wide", "long", "list")) {
  format <- match.arg(format)
  time_col <- object$time_col
  status_col <- object$status_col

  n_per_model <- vapply(object$model_results, function(res) {
    length(res$imputed_datasets)
  }, integer(1))
  n_common <- min(n_per_model)

  if (n_common < 1) {
    stop("No completed datasets available for model comparison output")
  }

  if (is.null(dataset)) {
    dataset_idx <- if (identical(format, "wide")) 1L else seq_len(n_common)
  } else {
    if (!is.numeric(dataset) || any(dataset < 1) || any(dataset > n_common)) {
      stop("Dataset index must be between 1 and ", n_common)
    }
    dataset_idx <- as.integer(dataset)
  }

  combine_one <- function(imp_id) {
    per_model <- list()
    for (nm in object$model_names) {
      d <- object$model_results[[nm]]$imputed_datasets[[imp_id]]
      d$.model <- nm
      d$.imp <- imp_id
      if (!(".id" %in% names(d))) {
        d$.id <- seq_len(nrow(d))
      }
      per_model[[nm]] <- d
    }
    all_cols <- unique(unlist(lapply(per_model, names)))
    aligned <- lapply(per_model, function(d) {
      missing_cols <- setdiff(all_cols, names(d))
      if (length(missing_cols) > 0) {
        for (col in missing_cols) {
          d[[col]] <- NA
        }
      }
      d[, all_cols, drop = FALSE]
    })
    out <- do.call(rbind, aligned)
    standardize_complete_column_order(
      out,
      time_col = time_col,
      status_col = status_col
    )
  }

  if (identical(format, "long")) {
    long_data <- do.call(rbind, lapply(dataset_idx, combine_one))
    return(standardize_complete_column_order(
      long_data,
      time_col = time_col,
      status_col = status_col
    ))
  }

  out <- lapply(dataset_idx, combine_one)
  names(out) <- paste0("dataset_", dataset_idx)
  if (identical(format, "wide") && length(out) == 1L) {
    return(out[[1]])
  }
  out
}
