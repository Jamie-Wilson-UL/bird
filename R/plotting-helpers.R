#' @title Plotting Helper Functions
#' @description Internal functions for creating various plots for bayesian_imputation objects
#' @keywords internal
#' @name plotting-helpers
NULL

# Shared visual helpers to keep plot styling consistent across all plot types.
theme_bird <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      plot.subtitle = ggplot2::element_text(size = base_size, color = "gray35"),
      plot.caption = ggplot2::element_text(size = base_size - 2, color = "gray45"),
      axis.title = ggplot2::element_text(size = base_size + 1),
      axis.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size),
      legend.text = ggplot2::element_text(size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}

style_bird_plot <- function(p, legend_position = "right", legend_justification = NULL, base_size = 12) {
  p <- p + theme_bird(base_size = base_size) + ggplot2::theme(legend.position = legend_position)
  if (!is.null(legend_justification)) {
    p <- p + ggplot2::theme(legend.justification = legend_justification)
  }
  p
}

panel_title_bird <- function(text, base_size = 13) {
  grid::textGrob(
    label = text,
    x = 0,
    hjust = 0,
    gp = grid::gpar(fontsize = base_size, fontface = "bold", col = "gray20")
  )
}

resolve_group_colors <- function(group_names, group_colors = NULL) {
  n <- length(group_names)
  default_cols <- grDevices::hcl.colors(n, "Dark 3")
  names(default_cols) <- group_names

  if (is.null(group_colors)) {
    return(default_cols)
  }

  if (!is.null(names(group_colors)) && any(names(group_colors) %in% group_names)) {
    out <- default_cols
    nm <- intersect(names(group_colors), group_names)
    out[nm] <- group_colors[nm]
    return(out)
  }

  group_colors <- as.character(group_colors)
  if (length(group_colors) < n) {
    group_colors <- rep(group_colors, length.out = n)
  }
  out <- group_colors[seq_len(n)]
  names(out) <- group_names
  out
}

#' Plot survival curves 
#' @param x bayesian_imputation object
#' @param n_curves Number of imputed curves to show (default: 10)
#' @param alpha Transparency for imputed curves (default: 0.3)
#' @param show_original Whether to show original Kaplan-Meier curves (default: TRUE)
#' @return ggplot object
#' @keywords internal
plot_survival_curves <- function(x, n_curves = 10, alpha = 0.3, show_original = TRUE, ...) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for Kaplan-Meier curves")
  }
  
  # Column names
  time_var <- x$time_col
  status_var <- x$status_col
  
  # Original KM (for optional overlay)
  km_orig <- survival::survfit(
    as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
    data = x$original_data
  )
  orig_df <- data.frame(time = km_orig$time, survival = km_orig$surv, type = "KM estimation")
  
  # Sample imputed datasets for plotting (representative sample)
  n_show <- min(n_curves, length(x$imputed_datasets))
  selected_datasets <- sample(x$imputed_datasets, n_show)
  
  # Compute survival curves for imputed datasets
  imputed_curves <- vector("list", n_show)
  
  for (i in 1:n_show) {
    km_imp <- survival::survfit(
      as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
      data = selected_datasets[[i]]
    )
    
    imputed_curves[[i]] <- data.frame(
      time = km_imp$time,
      survival = km_imp$surv,
      type = "Bayesian imputation",
      dataset_id = i
    )
  }
  
  imputed_df <- do.call(rbind, imputed_curves)
  
  # Combine data (keep dataset_id for grouping)
  plot_data <- imputed_df[, c("time", "survival", "type", "dataset_id")]
  
  # Create plot 
  p <- ggplot2::ggplot() +
    # Imputed curves - thin gray lines with transparency
    ggplot2::geom_step(
      data = plot_data,
      ggplot2::aes(x = time, y = survival, color = type, group = interaction(type, dataset_id)),
      alpha = alpha, linewidth = 0.6
    ) +
    # Original curve - black line (if requested)
    {if (show_original) ggplot2::geom_step(
      data = orig_df,
      ggplot2::aes(x = time, y = survival, color = type),
      linewidth = 0.7
    )} +
    ggplot2::scale_color_manual(values = c("KM estimation" = "black", "Bayesian imputation" = "#7a7a7a")) +
    ggplot2::labs(
      title = paste0("Bayesian Imputation (n = ", n_show, ") vs Kaplan-Meier"),
      subtitle = NULL,
      x = paste0("Time (", x$model_info$time_unit %||% "days", ")"), 
      y = "Survival probability", color = NULL
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))

  p <- style_bird_plot(
    p,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1),
    base_size = 13
  )
  
  # Add median survival lines from imputed datasets (mean median across imputations)
  all_medians <- sapply(x$imputed_datasets, function(ds) {
    fit <- survival::survfit(as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")), data = ds)
    idx <- which(fit$surv <= 0.5)[1]
    if (!is.na(idx)) fit$time[idx] else max(fit$time)
  })
  median_mean <- mean(all_medians, na.rm = TRUE)
  median_q025 <- stats::quantile(all_medians, 0.025, na.rm = TRUE)
  median_q975 <- stats::quantile(all_medians, 0.975, na.rm = TRUE)
  if (is.finite(median_mean)) {
    hseg <- data.frame(x_start = 0, y = 0.5, x_end = median_mean, y_end = 0.5)
    vseg <- data.frame(x_start = median_mean, y = 0.5, x_end = median_mean, y_end = 0)
    p <- p +
      ggplot2::geom_segment(
        data = hseg,
        ggplot2::aes(x = x_start, y = y, xend = x_end, yend = y_end),
        inherit.aes = FALSE, linetype = "dotted", linewidth = 0.8, alpha = 0.7, color = "red"
      ) +
      ggplot2::geom_segment(
        data = vseg,
        ggplot2::aes(x = x_start, y = y, xend = x_end, yend = y_end),
        inherit.aes = FALSE, linetype = "dotted", linewidth = 0.8, alpha = 0.7, color = "red"
      )
  }
  
  # Caption with median summary (match groups style)
  cap_txt <- if (is.finite(median_mean)) {
    unit <- x$model_info$time_unit %||% "days"
    paste0("Median survival time: ", round(median_mean, 1), " ", unit, " (95% CI: ",
           round(median_q025, 1), " - ", round(median_q975, 1), ")")
  } else {
    NULL
  }
  
  if (!is.null(cap_txt)) {
    p <- p + ggplot2::labs(caption = cap_txt) +
      ggplot2::theme(plot.caption = ggplot2::element_text(size = 10, color = "gray60", hjust = 0))
  }
  
  return(p)
}

#' Plot boxplots comparison of all completed datasets
#' @param x bayesian_imputation object
#' @param n_max Number of datasets to sample for boxplots (default: 10)
#' @param dataset_indices Optional indices of datasets to include (random if NULL)
#' @return ggplot object
#' @keywords internal
plot_boxplots_comparison <- function(x, n_max = 10, dataset_indices = NULL, ...) {
  
  # Extract time variable name
  time_var <- x$time_col
  
  # Select at most n_max random datasets
  total <- length(x$imputed_datasets)
  if (is.null(dataset_indices)) {
    k <- min(n_max, total)
    dataset_indices <- sort(sample(seq_len(total), k))
  }
  
  # Prepare data for selected datasets
  all_data <- vector("list", length(dataset_indices))
  for (j in seq_along(dataset_indices)) {
    i <- dataset_indices[j]
    ds <- x$imputed_datasets[[i]]
    all_data[[j]] <- data.frame(
      time = ds[[time_var]],
      dataset = factor(paste("Dataset", i), levels = paste("Dataset", dataset_indices)),
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all datasets
  plot_data <- do.call(rbind, all_data)
  
  # Create boxplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = dataset, x = time)) +
    ggplot2::geom_boxplot(fill = "lightblue", alpha = 0.7, outlier.alpha = 0.3) +
    ggplot2::labs(
      title = "Distribution of Survival Times Across Imputed Datasets",
      subtitle = paste("Comparing", length(dataset_indices), "sampled datasets (of", total, ")"),
      x = "Survival Time",
      y = "Imputed Dataset"
    )

  p <- style_bird_plot(p, legend_position = "none", base_size = 12)
  
  return(p)
}

plot_boxplots_comparison_groups <- function(x, n_max = 10, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for arranging plots")
  }
  # Use the same dataset indices across groups when possible
  num_per_group <- sapply(x$group_results, function(gr) length(gr$imputed_datasets))
  k <- min(n_max, min(num_per_group))
  dataset_indices <- sort(sample(seq_len(min(num_per_group)), k))
  
  # Compute global x-range across both groups and selected datasets
  # Use each group's time column for extraction
  all_times <- c()
  for (gname in x$group_names) {
    gr <- x$group_results[[gname]]
    tcol <- gr$time_col
    for (i in dataset_indices) {
      all_times <- c(all_times, gr$imputed_datasets[[i]][[tcol]])
    }
  }
  all_times <- all_times[is.finite(all_times)]
  xr <- range(all_times)
  if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) {
    # Fallback/expand if needed
    xr <- c(0, max(all_times, na.rm = TRUE) * 1.05)
    if (!is.finite(xr[2])) xr <- c(0, 1)
  }
  
  grobs <- list()
  for (gname in x$group_names) {
    gr <- x$group_results[[gname]]
    p <- plot_boxplots_comparison(gr, n_max = n_max, dataset_indices = dataset_indices, ...)
    p <- p + ggplot2::labs(title = paste0("Group: ", gname, " -- Boxplots (", k, ")"), subtitle = NULL)
    p <- p + ggplot2::coord_cartesian(xlim = xr)
    grobs[[length(grobs) + 1]] <- p
  }
  ncol <- min(2, length(grobs))
  arranged <- gridExtra::grid.arrange(grobs = grobs, ncol = ncol,
                                      top = panel_title_bird("Distribution of Survival Times Across Imputed Datasets -- Groups"))
  return(arranged)
}

#' Plot MCMC trace plots for model parameters
#' @param x bayesian_imputation object
#' @return ggplot object or NULL
#' @keywords internal
plot_trace_plots <- function(x, ...) {
  
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    message("Trace plots require the 'bayesplot' package. Install with: install.packages('bayesplot')")
    return(invisible(NULL))
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    message("Trace plots require the 'posterior' package. Install with: install.packages('posterior')")
    return(invisible(NULL))
  }
  
  # Get posterior draws of model parameters
  draws <- x$posterior_samples
  
  if (is.null(draws)) {
    message("No posterior draws available for trace plots")
    return(invisible(NULL))
  }
  
  # Convert to a format bayesplot understands
  draws_array <- tryCatch(posterior::as_draws_array(draws), error = function(e) NULL)
  if (is.null(draws_array)) {
    message("Could not convert posterior samples to draws array for plotting")
    return(invisible(NULL))
  }
  
  # bayesplot::mcmc_trace returns a ggplot that can be modified
  p <- mcmc_trace_local(draws_array)
  
  return(p)
}

#' Plot parameter correlation pairs
#' @param x bayesian_imputation object
#' @return ggplot object or NULL
#' @keywords internal
plot_pairs <- function(x, ...) {
  
  if (!requireNamespace("bayesplot", quietly = TRUE)) {
    message("Pairs plots require the 'bayesplot' package. Install with: install.packages('bayesplot')")
    return(invisible(NULL))
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    message("Pairs plots require the 'posterior' package. Install with: install.packages('posterior')")
    return(invisible(NULL))
  }
  
  # Get posterior draws
  draws <- x$posterior_samples
  if (is.null(draws)) {
    message("No posterior draws available for pairs plots")
    return(invisible(NULL))
  }
  draws_array <- tryCatch(posterior::as_draws_array(draws), error = function(e) NULL)
  if (is.null(draws_array)) {
    message("Could not convert posterior samples to draws array for plotting")
    return(invisible(NULL))
  }
  
  # bayesplot::mcmc_pairs returns a grid object, not a single ggplot
  p <- mcmc_pairs_local(draws_array)
  
  return(p)
}

#' Plot posterior distribution for censored observations
#' @param result bayesian_imputation object
#' @param obs_id Specific observation ID to plot (NULL for all)
#' @return ggplot object
#' @keywords internal
plot_posterior_censored <- function(result, obs_id = NULL, dataset_id = NULL, show_title = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }

  # Validate posterior draws availability
  if (!("posterior_imputations" %in% names(result)) || !is.matrix(result$posterior_imputations)) {
    stop("posterior_imputations not available in result; cannot plot posterior distribution")
  }

  # All censored observation indices in original data (status == 0)
  all_censored_indices <- which(result$original_data[[result$status_col]] == 0)
  if (length(all_censored_indices) == 0) {
    stop("No censored observations found in the data.")
  }

  # Select observation: user-provided, else random censored obs
  if (is.null(obs_id)) {
    obs_id <- sample(all_censored_indices, 1)
  } else if (!(obs_id %in% all_censored_indices)) {
    stop("Observation ", obs_id, " is not censored.")
  }

  # Column position in posterior_imputations corresponding to this obs_id
  censored_position <- which(all_censored_indices == obs_id)
  if (length(censored_position) != 1) {
    stop("Could not locate censored observation column in posterior_imputations")
  }

  # Extract posterior draws for this censored observation
  posterior_draws <- result$posterior_imputations[, censored_position]
  n_draws <- length(posterior_draws)

  # Density info for layout helpers
  dens_est <- stats::density(posterior_draws, adjust = 1.5)
  y_max <- max(dens_est$y, na.rm = TRUE)
  x_min <- min(posterior_draws, na.rm = TRUE)
  x_max <- max(posterior_draws, na.rm = TRUE)
  x_rng <- x_max - x_min
  x_off <- if (is.finite(x_rng) && x_rng > 0) 0.02 * x_rng else 0.02

  # Censoring time from original data
  censor_time <- result$original_data[[result$time_col]][obs_id]

  # Choose which completed dataset to use for indicating the imputed time
  num_completed <- length(result$imputed_datasets)
  if (num_completed < 1) {
    stop("No completed (imputed) datasets available in result")
  }
  if (is.null(dataset_id)) {
    dataset_id <- sample(seq_len(num_completed), 1)
  } else if (dataset_id < 1 || dataset_id > num_completed) {
    stop("Invalid dataset_id; must be between 1 and ", num_completed)
  }
  imputed_time <- result$imputed_datasets[[dataset_id]][[result$time_col]][obs_id]

  # Choose label side (right by default, left if near boundary)
  place_right <- is.finite(imputed_time) && (imputed_time + x_off <= x_max)
  label_x <- if (place_right) imputed_time + x_off else imputed_time - x_off
  label_hjust <- if (place_right) 0 else 1

  # Data frame for plotting posterior distribution
  plot_data <- data.frame(imputed_time_draw = posterior_draws)

  # Build plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = imputed_time_draw)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)),
                            bins = 30, alpha = 0.7,
                            fill = "lightblue", color = "darkblue") +
    ggplot2::geom_density(color = "red", linewidth = 1.2, adjust = 1.5)

  # Add vertical reference lines with legend (censoring vs imputed)
  lines_df <- data.frame(
    xintercept = c(censor_time, imputed_time),
    line_type = c("Censoring time", "Imputed time"),
    stringsAsFactors = FALSE
  )

  p <- p +
    ggplot2::geom_vline(
      data = lines_df,
      ggplot2::aes(xintercept = xintercept, linetype = line_type, color = line_type),
      linewidth = 1.2,
      key_glyph = "path"
    ) +
    ggplot2::scale_color_manual(
      values = c("Censoring time" = "darkred", "Imputed time" = "darkgreen"),
      name = NULL,
      labels = c(
        "Censoring time" = paste0("Censoring time = ", round(censor_time, 2)),
        "Imputed time"   = paste0("Imputed time = ", round(imputed_time, 2))
      )
    ) +
    ggplot2::scale_linetype_manual(
      values = c("Censoring time" = "dashed", "Imputed time" = "solid"),
      name = NULL,
      labels = c(
        "Censoring time" = paste0("Censoring time = ", round(censor_time, 2)),
        "Imputed time"   = paste0("Imputed time = ", round(imputed_time, 2))
      )
    )

  # Label the imputed time near the top of the line
  p <- p + ggplot2::annotate(
    "text", x = label_x, y = y_max * 0.94,
    label = paste0(round(imputed_time, 2)),
    color = "darkgreen", angle = 0, vjust = 0.5, hjust = label_hjust, size = 3.2
  )

  if (show_title) {
    p <- p + ggplot2::labs(
      title = paste("Posterior Distribution for Censored Observation", obs_id),
      subtitle = paste("Full posterior distribution (n =", n_draws, "MCMC draws)",
                       "| Imputed time from dataset", dataset_id),
      x = "Imputed Event Time",
      y = "Density",
      caption = paste0(
        "Dashed red: censoring time (", round(censor_time, 2), ") | ",
        "Solid green: imputed time (", round(imputed_time, 2), ")"
      )
    )
  } else {
    # Add a small in-plot note indicating which observation was chosen
    p <- p + ggplot2::annotate("text", x = x_min + x_off, y = y_max * 0.98,
                               label = paste0("Obs ", obs_id),
                               color = "gray30", size = 3.2, hjust = 0, vjust = 1)

    p <- p + ggplot2::labs(
      x = "Imputed Event Time",
      y = "Density",
      caption = paste0(
        "Dashed red: censoring time (", round(censor_time, 2), ") | ",
        "Solid green: imputed time (", round(imputed_time, 2), ")"
      )
    )
  }

  p <- style_bird_plot(
    p,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1),
    base_size = 11
  ) +
    ggplot2::theme(
      legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", alpha.f = 0.75), color = "grey85"),
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 8)
    )

  return(p)
}

#' Normalise requested panels for completed dataset summary plots
#' @param panels Panel selection supplied by user.
#' @param default_panels Character vector used when `panels = "auto"`.
#' @param context_label Label used in validation messages.
#' @keywords internal
normalize_completed_summary_panels <- function(panels = "auto",
                                               default_panels = c("hist", "density", "survival", "boxplot"),
                                               context_label = "completed_dataset_summary") {
  valid_panels <- c("hist", "density", "survival", "boxplot")

  if (is.null(panels) || (length(panels) == 1 && identical(tolower(as.character(panels)), "auto"))) {
    return(default_panels)
  }

  panels <- as.character(panels)
  if (length(panels) < 1) {
    stop("`panels` must contain at least one panel name")
  }

  alias <- c(
    "histogram" = "hist",
    "histograms" = "hist",
    "boxplots" = "boxplot",
    "box" = "boxplot",
    "surv" = "survival",
    "km" = "survival"
  )

  mapped <- tolower(panels)
  mapped <- ifelse(mapped %in% names(alias), alias[mapped], mapped)
  mapped <- unique(mapped)

  invalid <- setdiff(mapped, valid_panels)
  if (length(invalid) > 0) {
    stop(
      "Invalid `panels` for ", context_label, ": ",
      paste(invalid, collapse = ", "),
      ". Valid options are: ",
      paste(valid_panels, collapse = ", "),
      ", or 'auto'."
    )
  }

  mapped
}

#' Arrange selected completed-summary panels
#' @param panel_grobs List of ggplot grobs.
#' @param title_text Title shown above combined panel.
#' @keywords internal
arrange_completed_summary_panels <- function(panel_grobs, title_text) {
  n <- length(panel_grobs)
  if (n < 1) {
    stop("At least one panel must be selected")
  }

  ncol <- if (n == 1) 1 else if (n == 2) 2 else 2
  do.call(
    gridExtra::grid.arrange,
    c(panel_grobs, list(ncol = ncol, top = panel_title_bird(title_text)))
  )
}

#' Plot completed dataset summary panels
#' @param x bayesian_imputation object
#' @param dataset_id Specific dataset to plot (NULL for random)
#' @param panels Which panels to display. Use `"auto"` for default 4-panel layout,
#'   or any combination of `c("hist", "density", "survival", "boxplot")`.
#' @param ... Additional arguments
#' @return ggplot object
#' @keywords internal
plot_completed_dataset_summary <- function(x, dataset_id = NULL, panels = "auto", ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  # Select dataset
  if (is.null(dataset_id)) {
    dataset_id <- sample(seq_along(x$imputed_datasets), 1)
  }
  
  if (dataset_id < 1 || dataset_id > length(x$imputed_datasets)) {
    stop("Invalid dataset_id")
  }
  
  # Get the selected dataset
  dataset <- x$imputed_datasets[[dataset_id]]
  time_var <- x$time_col
  status_var <- x$status_col
  
  selected_panels <- normalize_completed_summary_panels(
    panels = panels,
    default_panels = c("hist", "density", "survival", "boxplot"),
    context_label = "plot_completed_dataset_summary"
  )

  panel_grobs <- list()

  if ("hist" %in% selected_panels) {
    p_hist <- ggplot2::ggplot(dataset, ggplot2::aes(x = .data[[time_var]])) +
      ggplot2::geom_histogram(bins = 30, fill = "lightblue", alpha = 0.7) +
      ggplot2::labs(title = "Histogram", x = "Time", y = "Count")
    p_hist <- style_bird_plot(p_hist, legend_position = "none", base_size = 11)
    panel_grobs$hist <- p_hist
  }

  if ("density" %in% selected_panels) {
    if (requireNamespace("logspline", quietly = TRUE)) {
      x_vals <- dataset[[time_var]]
      x_vals <- x_vals[is.finite(x_vals) & x_vals >= 0]
      x_vals <- pmax(x_vals, .Machine$double.eps)

      ls_fit <- tryCatch(
        logspline::logspline(x_vals, lbound = 0, maxknots = 6),
        error = function(e) tryCatch(
          logspline::logspline(x_vals, lbound = 0),
          error = function(e2) NULL
        )
      )

      orig_times <- x$original_data[[x$time_col]]
      orig_status <- x$original_data[[x$status_col]]
      orig_events <- orig_times[is.finite(orig_times) & orig_times >= 0 & orig_status == 1]
      all_times <- c(orig_events, dataset[[time_var]])
      all_times <- all_times[is.finite(all_times) & all_times >= 0]
      all_times <- pmax(all_times, .Machine$double.eps)

      if (length(all_times) >= 2) {
        qlo <- max(0, stats::quantile(all_times, 0.005, na.rm = TRUE))
        qhi <- stats::quantile(all_times, 0.995, na.rm = TRUE)
      } else {
        qlo <- min(all_times, na.rm = TRUE)
        qhi <- max(all_times, na.rm = TRUE)
      }
      if (!is.finite(qlo) || !is.finite(qhi) || qhi <= qlo) {
        qlo <- min(all_times, na.rm = TRUE)
        qhi <- max(all_times, na.rm = TRUE)
      }
      x_grid <- seq(qlo, qhi, length.out = 400)

      if (!is.null(ls_fit)) {
        y_density <- logspline::dlogspline(x_grid, ls_fit)
      } else {
        dens <- stats::density(x_vals, from = min(x_grid), to = max(x_grid), n = length(x_grid))
        y_density <- stats::approx(dens$x, dens$y, xout = x_grid, rule = 2)$y
      }

      p_density <- ggplot2::ggplot(data.frame(x = x_grid, y = y_density), ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line(color = "darkblue", linewidth = 1) +
        ggplot2::labs(title = "Density", x = "Time", y = "Density")
      p_density <- style_bird_plot(p_density, legend_position = "none", base_size = 11)
    } else {
      p_density <- ggplot2::ggplot(dataset, ggplot2::aes(x = .data[[time_var]])) +
        ggplot2::geom_density(fill = "lightgreen", alpha = 0.7) +
        ggplot2::labs(title = "Density", x = "Time", y = "Density")
      p_density <- style_bird_plot(p_density, legend_position = "none", base_size = 11)
    }
    panel_grobs$density <- p_density
  }

  if ("survival" %in% selected_panels) {
    km_fit <- survival::survfit(
      as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
      data = dataset
    )
    p_survival <- ggplot2::ggplot(data.frame(time = km_fit$time, survival = km_fit$surv),
                                  ggplot2::aes(x = time, y = survival)) +
      ggplot2::geom_step(color = "red", linewidth = 1) +
      ggplot2::labs(title = "Survival Curve", x = "Time", y = "Survival Probability")
    p_survival <- style_bird_plot(p_survival, legend_position = "none", base_size = 11)
    panel_grobs$survival <- p_survival
  }

  if ("boxplot" %in% selected_panels) {
    p_boxplot <- ggplot2::ggplot(dataset, ggplot2::aes(x = "", y = .data[[time_var]])) +
      ggplot2::geom_boxplot(fill = "lightcoral", alpha = 0.7) +
      ggplot2::labs(title = "Boxplot", x = "", y = "Survival Time")
    p_boxplot <- style_bird_plot(p_boxplot, legend_position = "none", base_size = 11)
    panel_grobs$boxplot <- p_boxplot
  }

  ordered_panels <- panel_grobs[selected_panels]
  arrange_completed_summary_panels(
    ordered_panels,
    paste("Completed Dataset", dataset_id, "Summary")
  )
}

#' Plot completed dataset summary for groups (overlaid comparison)
#' @param x bayesian_imputation_groups object
#' @param dataset_id Dataset ID to plot (NULL for random)
#' @param panels Which panels to display. Use `"auto"` for default behaviour
#'   (full 4-panel for <=2 groups; survival+boxplot for >2 groups),
#'   or any combination of `c("hist", "density", "survival", "boxplot")`.
#' @param alpha Transparency for overlapping plots
#' @param group_colors Optional color palette for groups (auto-generated if NULL)
#' @param ... Additional arguments
#' @keywords internal
plot_completed_dataset_summary_groups <- function(x, dataset_id = NULL, panels = "auto", alpha = 0.6,
                                                  group_colors = NULL, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  # Select dataset ID if not provided
  if (is.null(dataset_id)) {
    # Use the same dataset ID across all groups for fair comparison
    max_datasets <- min(sapply(x$group_results, function(g) length(g$imputed_datasets)))
    dataset_id <- sample(1:max_datasets, 1)
  }
  
  # Check if dataset_id is valid for all groups
  for (group_name in x$group_names) {
    if (dataset_id < 1 || dataset_id > length(x$group_results[[group_name]]$imputed_datasets)) {
      stop("Invalid dataset_id for group: ", group_name)
    }
  }

  group_color_map <- resolve_group_colors(x$group_names, group_colors)
  
  # Get datasets for each group
  group_datasets <- list()
  for (group_name in x$group_names) {
    group_datasets[[group_name]] <- x$group_results[[group_name]]$imputed_datasets[[dataset_id]]
    # Add group column for plotting
    group_datasets[[group_name]]$group <- group_name
  }
  
  # Combine all datasets for plotting
  all_data <- do.call(rbind, group_datasets)
  time_var <- x$group_results[[1]]$time_col
  status_var <- x$group_results[[1]]$status_col
  
  default_panels <- if (length(x$group_names) > 2) {
    c("survival", "boxplot")
  } else {
    c("hist", "density", "survival", "boxplot")
  }
  selected_panels <- normalize_completed_summary_panels(
    panels = panels,
    default_panels = default_panels,
    context_label = "plot_completed_dataset_summary_groups"
  )
  is_auto_panels <- is.null(panels) ||
    (length(panels) == 1 && identical(tolower(as.character(panels)), "auto"))

  if (length(x$group_names) > 2 && !is_auto_panels &&
      any(selected_panels %in% c("hist", "density"))) {
    warning("Using histogram/density with >2 groups may be visually crowded")
  }

  panel_grobs <- list()

  if ("hist" %in% selected_panels) {
    p_hist <- ggplot2::ggplot(all_data, ggplot2::aes(x = .data[[time_var]], fill = group)) +
      ggplot2::geom_histogram(bins = 30, alpha = alpha, position = "identity") +
      ggplot2::scale_fill_manual(values = group_color_map) +
      ggplot2::labs(title = "Histogram", x = "Time", y = "Count", fill = "Group")
    p_hist <- style_bird_plot(p_hist, base_size = 11)
    panel_grobs$hist <- p_hist
  }

  if ("density" %in% selected_panels) {
    if (requireNamespace("logspline", quietly = TRUE)) {
      time_col <- time_var
      status_col <- status_var
      orig_times <- x$original_data[[time_col]]
      orig_status <- x$original_data[[status_col]]
      orig_events <- orig_times[is.finite(orig_times) & orig_times >= 0 & orig_status == 1]
      pooled_times <- c(orig_events, all_data[[time_var]])
      pooled_times <- pooled_times[is.finite(pooled_times) & pooled_times >= 0]
      pooled_times <- pmax(pooled_times, .Machine$double.eps)

      if (length(pooled_times) >= 2) {
        qlo <- max(0, stats::quantile(pooled_times, 0.005, na.rm = TRUE))
        qhi <- stats::quantile(pooled_times, 0.995, na.rm = TRUE)
      } else {
        qlo <- min(pooled_times, na.rm = TRUE)
        qhi <- max(pooled_times, na.rm = TRUE)
      }
      if (!is.finite(qlo) || !is.finite(qhi) || qhi <= qlo) {
        qlo <- min(pooled_times, na.rm = TRUE)
        qhi <- max(pooled_times, na.rm = TRUE)
      }
      x_grid <- seq(qlo, qhi, length.out = 400)

      density_data <- data.frame()
      for (i in seq_along(x$group_names)) {
        group_name <- x$group_names[i]
        group_data <- group_datasets[[group_name]]
        x_vals <- group_data[[time_var]]
        x_vals <- x_vals[is.finite(x_vals) & x_vals >= 0]
        x_vals <- pmax(x_vals, .Machine$double.eps)
        ls_fit <- tryCatch(
          logspline::logspline(x_vals, lbound = 0, maxknots = 6),
          error = function(e) tryCatch(
            logspline::logspline(x_vals, lbound = 0),
            error = function(e2) NULL
          )
        )

        if (!is.null(ls_fit)) {
          y_density <- logspline::dlogspline(x_grid, ls_fit)
        } else {
          dens <- stats::density(x_vals, from = min(x_grid), to = max(x_grid), n = length(x_grid))
          y_density <- stats::approx(dens$x, dens$y, xout = x_grid, rule = 2)$y
        }

        density_data <- rbind(density_data, data.frame(
          x = x_grid,
          y = y_density,
          group = group_name,
          stringsAsFactors = FALSE
        ))
      }

      p_density <- ggplot2::ggplot(density_data, ggplot2::aes(x = x, y = y, color = group)) +
        ggplot2::geom_line(linewidth = 1, alpha = alpha) +
        ggplot2::scale_color_manual(values = group_color_map) +
        ggplot2::labs(title = "Density", x = "Time", y = "Density", color = "Group")
      p_density <- style_bird_plot(p_density, base_size = 11)
    } else {
      p_density <- ggplot2::ggplot(all_data, ggplot2::aes(x = .data[[time_var]], fill = group)) +
        ggplot2::geom_density(alpha = alpha) +
        ggplot2::scale_fill_manual(values = group_color_map) +
        ggplot2::labs(title = "Density", x = "Time", y = "Density", fill = "Group")
      p_density <- style_bird_plot(p_density, base_size = 11)
    }
    panel_grobs$density <- p_density
  }

  if ("survival" %in% selected_panels) {
    survival_data <- data.frame()
    for (i in seq_along(x$group_names)) {
      group_name <- x$group_names[i]
      group_data <- group_datasets[[group_name]]
      km_fit <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = group_data
      )
      survival_data <- rbind(survival_data, data.frame(
        time = km_fit$time,
        survival = km_fit$surv,
        group = group_name,
        stringsAsFactors = FALSE
      ))
    }

    p_survival <- ggplot2::ggplot(survival_data, ggplot2::aes(x = time, y = survival, color = group)) +
      ggplot2::geom_step(linewidth = 1, alpha = alpha) +
      ggplot2::scale_color_manual(values = group_color_map) +
      ggplot2::labs(title = "Survival Curve", x = "Time", y = "Survival Probability", color = "Group")
    p_survival <- style_bird_plot(p_survival, base_size = 11)
    panel_grobs$survival <- p_survival
  }

  if ("boxplot" %in% selected_panels) {
    p_boxplot <- ggplot2::ggplot(all_data, ggplot2::aes(x = group, y = .data[[time_var]], fill = group)) +
      ggplot2::geom_boxplot(alpha = alpha) +
      ggplot2::scale_fill_manual(values = group_color_map) +
      ggplot2::labs(title = "Boxplots", x = "Group", y = "Survival Time", fill = "Group")
    p_boxplot <- style_bird_plot(p_boxplot, base_size = 11)
    panel_grobs$boxplot <- p_boxplot
  }

  ordered_panels <- panel_grobs[selected_panels]
  arrange_completed_summary_panels(
    ordered_panels,
    paste("Completed Dataset", dataset_id, "Summary - Group Comparison")
  )
}

#' Plot survival curves across compared models
#' @param x bird_model_comparison object
#' @param show_original Whether to overlay Kaplan-Meier from original observed data
#' @param n_grid Number of grid points used to summarise model curves
#' @param alpha Ribbon transparency for model uncertainty bands
#' @param model_colors Optional named/un-named vector of colors for models
#' @param ... Additional arguments (unused)
#' @keywords internal
plot_survival_curves_models <- function(x, show_original = TRUE, n_grid = 200, alpha = 0.18, model_colors = NULL, ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for Kaplan-Meier curves")
  }

  model_names <- x$model_names
  color_map <- resolve_model_colors(model_names, model_colors)

  all_times <- c()
  for (nm in model_names) {
    fit <- x$model_results[[nm]]
    for (j in seq_along(fit$imputed_datasets)) {
      all_times <- c(all_times, fit$imputed_datasets[[j]][[fit$time_col]])
    }
  }
  all_times <- all_times[is.finite(all_times) & all_times >= 0]
  if (length(all_times) < 2) {
    stop("Not enough finite times found to build survival comparison")
  }

  t_min <- min(all_times)
  t_max <- max(all_times)
  if (!is.finite(t_min) || !is.finite(t_max)) {
    stop("Could not determine plotting range for model comparison")
  }
  if (identical(t_min, t_max)) {
    t_max <- t_max + 1
  }
  time_grid <- seq(t_min, t_max, length.out = max(50, as.integer(n_grid)))

  qfun <- function(v, prob) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(NA_real_)
    as.numeric(stats::quantile(v, prob, na.rm = TRUE))
  }

  summary_data <- data.frame()
  median_stats <- list()
  time_var <- x$time_col
  status_var <- x$status_col

  for (nm in model_names) {
    fit <- x$model_results[[nm]]
    n_sets <- length(fit$imputed_datasets)
    surv_mat <- matrix(NA_real_, nrow = length(time_grid), ncol = n_sets)
    medians <- rep(NA_real_, n_sets)

    for (j in seq_len(n_sets)) {
      ds <- fit$imputed_datasets[[j]]
      km <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = ds
      )

      surv_mat[, j] <- survival_step_at_grid(km, time_grid)
      idx <- which(km$surv <= 0.5)[1]
      medians[j] <- if (!is.na(idx)) km$time[idx] else max(km$time, na.rm = TRUE)
    }

    model_df <- data.frame(
      time = time_grid,
      surv_mean = rowMeans(surv_mat, na.rm = TRUE),
      surv_q025 = apply(surv_mat, 1, qfun, prob = 0.025),
      surv_q975 = apply(surv_mat, 1, qfun, prob = 0.975),
      model = nm,
      stringsAsFactors = FALSE
    )
    summary_data <- rbind(summary_data, model_df)

    median_stats[[nm]] <- list(
      mean = mean(medians, na.rm = TRUE),
      q025 = qfun(medians, 0.025),
      q975 = qfun(medians, 0.975)
    )
  }

  p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = time, y = surv_mean, color = model, fill = model)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = surv_q025, ymax = surv_q975), alpha = alpha, color = NA) +
    ggplot2::geom_step(linewidth = 1) +
    ggplot2::scale_color_manual(values = color_map) +
    ggplot2::scale_fill_manual(values = color_map) +
    ggplot2::labs(
      title = "Survival Curve Comparison Across Models",
      subtitle = "Line: mean survival across completed datasets | Ribbon: 95% interval",
      x = paste0("Time (", x$model_info$time_unit, ")"),
      y = "Survival probability",
      color = "Model",
      fill = "Model"
    )

  if (show_original) {
    km_orig <- survival::survfit(
      as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
      data = x$original_data
    )
    orig_df <- data.frame(time = km_orig$time, survival = km_orig$surv)
    p <- p + ggplot2::geom_step(
      data = orig_df,
      ggplot2::aes(x = time, y = survival),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.9,
      linetype = "22"
    )
  }

  med_caption <- paste(
    sapply(model_names, function(nm) {
      st <- median_stats[[nm]]
      paste0(
        nm, ": median ", round(st$mean, 1), " ",
        x$model_info$time_unit, " [", round(st$q025, 1), ", ", round(st$q975, 1), "]"
      )
    }),
    collapse = " | "
  )

  p <- style_bird_plot(p, base_size = 12)

  p + ggplot2::labs(caption = med_caption) +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 9, color = "gray40", hjust = 0))
}

#' Plot boxplots comparison across models
#' @param x bird_model_comparison object
#' @param n_max Number of datasets to sample (default: 10)
#' @param dataset_indices Optional dataset indices to include (shared across models)
#' @param model_colors Optional named/un-named vector of colors for models
#' @param ... Additional arguments (unused)
#' @keywords internal
plot_boxplots_comparison_models <- function(x, n_max = 10, dataset_indices = NULL, model_colors = NULL, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for arranging plots")
  }

  model_names <- x$model_names
  color_map <- resolve_model_colors(model_names, model_colors)
  time_var <- x$time_col

  n_per_model <- vapply(x$model_results, function(res) length(res$imputed_datasets), integer(1))
  n_common <- min(n_per_model)
  if (n_common < 1) {
    stop("No completed datasets available for model comparison")
  }

  if (is.null(dataset_indices)) {
    k <- min(n_max, n_common)
    dataset_indices <- sort(sample(seq_len(n_common), k))
  } else {
    if (!is.numeric(dataset_indices) || length(dataset_indices) < 1) {
      stop("dataset_indices must be a non-empty numeric vector")
    }
    dataset_indices <- sort(unique(as.integer(dataset_indices)))
    if (any(dataset_indices < 1) || any(dataset_indices > n_common)) {
      stop("dataset_indices must be between 1 and ", n_common)
    }
  }

  all_times <- c()
  for (nm in model_names) {
    for (i in dataset_indices) {
      all_times <- c(all_times, x$model_results[[nm]]$imputed_datasets[[i]][[time_var]])
    }
  }
  all_times <- all_times[is.finite(all_times)]
  xr <- range(all_times)
  if (!is.finite(xr[1]) || !is.finite(xr[2]) || identical(xr[1], xr[2])) {
    xr <- c(0, max(all_times, na.rm = TRUE) * 1.05)
    if (!is.finite(xr[2])) xr <- c(0, 1)
  }

  grobs <- list()
  for (nm in model_names) {
    fit <- x$model_results[[nm]]
    panel_data <- do.call(rbind, lapply(dataset_indices, function(i) {
      ds <- fit$imputed_datasets[[i]]
      data.frame(
        time = ds[[time_var]],
        dataset = factor(paste("Dataset", i), levels = paste("Dataset", dataset_indices)),
        stringsAsFactors = FALSE
      )
    }))

    p <- ggplot2::ggplot(panel_data, ggplot2::aes(y = dataset, x = time)) +
      ggplot2::geom_boxplot(fill = grDevices::adjustcolor(color_map[[nm]], alpha.f = 0.45),
                            outlier.alpha = 0.3) +
      ggplot2::coord_cartesian(xlim = xr) +
      ggplot2::labs(
        title = paste0("Model: ", nm),
        subtitle = paste("Boxplots for", length(dataset_indices), "datasets"),
        x = "Survival Time",
        y = "Imputed Dataset"
      )
    p <- style_bird_plot(p, legend_position = "none", base_size = 11)
    grobs[[length(grobs) + 1]] <- p
  }

  gridExtra::grid.arrange(
    grobs = grobs,
    ncol = min(2, length(grobs)),
    top = panel_title_bird("Distribution of Survival Times Across Imputed Datasets -- Models")
  )
}

#' Plot density comparison across models (single selected dataset)
#' @param x bird_model_comparison object
#' @param dataset_id Dataset index used across all models (NULL for random)
#' @param model_colors Optional named/un-named vector of colors for models
#' @param ... Additional arguments (unused)
#' @keywords internal
plot_density_comparison_models <- function(x, dataset_id = NULL, model_colors = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("logspline", quietly = TRUE)) {
    stop("logspline package required for density comparison")
  }

  model_names <- x$model_names
  color_map <- resolve_model_colors(model_names, model_colors)
  time_var <- x$time_col
  status_var <- x$status_col

  n_per_model <- vapply(x$model_results, function(res) length(res$imputed_datasets), integer(1))
  n_common <- min(n_per_model)
  if (n_common < 1) {
    stop("No completed datasets available for model comparison")
  }
  if (is.null(dataset_id)) {
    dataset_id <- sample(seq_len(n_common), 1)
  }
  if (!is.numeric(dataset_id) || length(dataset_id) != 1 || dataset_id < 1 || dataset_id > n_common) {
    stop("dataset_id must be between 1 and ", n_common)
  }
  dataset_id <- as.integer(dataset_id)

  selected_times <- c()
  for (nm in model_names) {
    selected_times <- c(selected_times, x$model_results[[nm]]$imputed_datasets[[dataset_id]][[time_var]])
  }
  selected_times <- selected_times[is.finite(selected_times) & selected_times >= 0]
  observed_events <- x$original_data[[time_var]][x$original_data[[status_var]] == 1]
  observed_events <- observed_events[is.finite(observed_events) & observed_events >= 0]
  if (length(observed_events) < 2) {
    observed_events <- x$original_data[[time_var]]
    observed_events <- observed_events[is.finite(observed_events) & observed_events >= 0]
  }
  if (length(observed_events) < 2) {
    stop("Need at least two observed times to compute model density comparison")
  }

  x_range <- range(c(observed_events, selected_times), finite = TRUE)
  x_range[1] <- max(0, x_range[1])
  x_grid <- seq(x_range[1], x_range[2], length.out = 300)

  density_on_grid <- function(values, grid) {
    values <- values[is.finite(values)]
    if (length(values) < 2) {
      return(rep(NA_real_, length(grid)))
    }
    ls_fit <- tryCatch(
      logspline::logspline(values, lbound = 0, maxknots = 6),
      error = function(e) tryCatch(
        logspline::logspline(values, lbound = 0),
        error = function(e2) NULL
      )
    )
    if (!is.null(ls_fit)) {
      return(logspline::dlogspline(grid, ls_fit))
    }
    dens <- stats::density(values, from = min(grid), to = max(grid), n = length(grid))
    stats::approx(dens$x, dens$y, xout = grid, rule = 2)$y
  }

  plot_data <- data.frame(
    x = x_grid,
    y = density_on_grid(observed_events, x_grid),
    curve = "Original",
    stringsAsFactors = FALSE
  )
  for (nm in model_names) {
    vals <- x$model_results[[nm]]$imputed_datasets[[dataset_id]][[time_var]]
    vals <- vals[is.finite(vals) & vals >= 0]
    vals <- pmax(vals, .Machine$double.eps)
    if (length(vals) < 2) next
    plot_data <- rbind(plot_data, data.frame(
      x = x_grid,
      y = density_on_grid(vals, x_grid),
      curve = nm,
      stringsAsFactors = FALSE
    ))
  }

  manual_cols <- c("Original" = "black", unname(color_map))
  names(manual_cols) <- c("Original", names(color_map))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = curve)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_color_manual(values = manual_cols) +
    ggplot2::labs(
      title = "Density Comparison Across Models",
      subtitle = paste0("Completed dataset ", dataset_id),
      x = paste0("Time (", x$model_info$time_unit %||% "days", ")"),
      y = "Density",
      color = NULL
    )
  style_bird_plot(
    p,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1),
    base_size = 12
  )
}

#' Plot completed dataset summary for model comparison
#' @param x bird_model_comparison object
#' @param dataset_id Dataset index to compare (same index used for all models)
#' @param panels Which panels to display. Use `"auto"` for default behaviour
#'   (full 4-panel for <=2 models; survival+boxplot for >2 models),
#'   or any combination of `c("hist", "density", "survival", "boxplot")`.
#' @param alpha Transparency for overlays
#' @param model_colors Optional named/un-named vector of colors for models
#' @param ... Additional arguments (unused)
#' @keywords internal
plot_completed_dataset_summary_models <- function(x, dataset_id = NULL, panels = "auto", alpha = 0.55, model_colors = NULL, ...) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for arranging plots")
  }
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for survival curves")
  }

  model_names <- x$model_names
  n_models <- length(model_names)
  color_map <- resolve_model_colors(model_names, model_colors)
  time_var <- x$time_col
  status_var <- x$status_col

  n_per_model <- vapply(x$model_results, function(res) {
    length(res$imputed_datasets)
  }, integer(1))
  n_common <- min(n_per_model)
  if (n_common < 1) {
    stop("No completed datasets available for model comparison")
  }

  if (is.null(dataset_id)) {
    dataset_id <- sample(seq_len(n_common), 1)
  }
  if (!is.numeric(dataset_id) || dataset_id < 1 || dataset_id > n_common) {
    stop("dataset_id must be between 1 and ", n_common)
  }
  dataset_id <- as.integer(dataset_id)

  model_datasets <- list()
  for (nm in model_names) {
    d <- x$model_results[[nm]]$imputed_datasets[[dataset_id]]
    d$.model <- nm
    model_datasets[[nm]] <- d
  }
  all_data <- do.call(rbind, lapply(model_datasets, function(d) {
    data.frame(
      .time = d[[time_var]],
      .status = d[[status_var]],
      .model = d$.model,
      stringsAsFactors = FALSE
    )
  }))

  default_panels <- if (n_models <= 2) {
    c("hist", "density", "survival", "boxplot")
  } else {
    c("survival", "boxplot")
  }
  selected_panels <- normalize_completed_summary_panels(
    panels = panels,
    default_panels = default_panels,
    context_label = "plot_completed_dataset_summary_models"
  )
  is_auto_panels <- is.null(panels) ||
    (length(panels) == 1 && identical(tolower(as.character(panels)), "auto"))
  if (n_models > 2 && !is_auto_panels && any(selected_panels %in% c("hist", "density"))) {
    warning("Using histogram/density with >2 models may be visually crowded")
  }

  panel_grobs <- list()

  if ("hist" %in% selected_panels) {
    p_hist <- ggplot2::ggplot(all_data, ggplot2::aes(x = .time, fill = .model)) +
      ggplot2::geom_histogram(bins = 30, alpha = alpha, position = "identity") +
      ggplot2::scale_fill_manual(values = color_map) +
      ggplot2::labs(title = "Histogram", x = "Time", y = "Count", fill = "Model")
    p_hist <- style_bird_plot(p_hist, base_size = 11)
    panel_grobs$hist <- p_hist
  }

  if ("density" %in% selected_panels) {
    # Use logspline curves (line-only) for consistency with group comparison style.
    pooled_times <- c(all_data$.time, x$original_data[[time_var]][x$original_data[[status_var]] == 1])
    pooled_times <- pooled_times[is.finite(pooled_times) & pooled_times >= 0]
    pooled_times <- pmax(pooled_times, .Machine$double.eps)

    if (length(pooled_times) >= 2) {
      qlo <- max(0, stats::quantile(pooled_times, 0.005, na.rm = TRUE))
      qhi <- stats::quantile(pooled_times, 0.995, na.rm = TRUE)
    } else {
      qlo <- min(pooled_times, na.rm = TRUE)
      qhi <- max(pooled_times, na.rm = TRUE)
    }
    if (!is.finite(qlo) || !is.finite(qhi) || qhi <= qlo) {
      qlo <- min(pooled_times, na.rm = TRUE)
      qhi <- max(pooled_times, na.rm = TRUE)
    }
    x_grid <- seq(qlo, qhi, length.out = 400)

    density_data <- data.frame()
    for (nm in model_names) {
      model_vals <- model_datasets[[nm]][[time_var]]
      model_vals <- model_vals[is.finite(model_vals) & model_vals >= 0]
      model_vals <- pmax(model_vals, .Machine$double.eps)

      y_density <- NULL
      if (requireNamespace("logspline", quietly = TRUE)) {
        ls_fit <- tryCatch(
          logspline::logspline(model_vals, lbound = 0, maxknots = 6),
          error = function(e) tryCatch(
            logspline::logspline(model_vals, lbound = 0),
            error = function(e2) NULL
          )
        )
        if (!is.null(ls_fit)) {
          y_density <- logspline::dlogspline(x_grid, ls_fit)
        }
      }
      if (is.null(y_density)) {
        dens <- stats::density(model_vals, from = min(x_grid), to = max(x_grid), n = length(x_grid))
        y_density <- stats::approx(dens$x, dens$y, xout = x_grid, rule = 2)$y
      }

      density_data <- rbind(density_data, data.frame(
        x = x_grid,
        y = y_density,
        model = nm,
        stringsAsFactors = FALSE
      ))
    }

    p_density <- ggplot2::ggplot(density_data, ggplot2::aes(x = x, y = y, color = model)) +
      ggplot2::geom_line(linewidth = 1, alpha = alpha) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::labs(title = "Density", x = "Time", y = "Density", color = "Model")
    p_density <- style_bird_plot(p_density, base_size = 11)
    panel_grobs$density <- p_density
  }

  if ("survival" %in% selected_panels) {
    surv_data <- data.frame()
    for (nm in model_names) {
      d <- model_datasets[[nm]]
      km <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = d
      )
      surv_data <- rbind(surv_data, data.frame(
        time = km$time,
        survival = km$surv,
        model = nm,
        stringsAsFactors = FALSE
      ))
    }

    p_survival <- ggplot2::ggplot(surv_data, ggplot2::aes(x = time, y = survival, color = model)) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::scale_color_manual(values = color_map) +
      ggplot2::labs(title = "Survival Curve", x = "Time", y = "Survival Probability", color = "Model")
    p_survival <- style_bird_plot(p_survival, base_size = 11)
    panel_grobs$survival <- p_survival
  }

  if ("boxplot" %in% selected_panels) {
    p_boxplot <- ggplot2::ggplot(all_data, ggplot2::aes(x = .model, y = .time, fill = .model)) +
      ggplot2::geom_boxplot(alpha = alpha) +
      ggplot2::scale_fill_manual(values = color_map) +
      ggplot2::labs(title = "Boxplots", x = "Model", y = "Survival Time", fill = "Model")
    p_boxplot <- style_bird_plot(p_boxplot, base_size = 11)
    panel_grobs$boxplot <- p_boxplot
  }

  ordered_panels <- panel_grobs[selected_panels]
  arrange_completed_summary_panels(
    ordered_panels,
    paste("Completed Dataset", dataset_id, "Summary - Model Comparison")
  )
}

#' Resolve colors for model comparisons
#' @param model_names Character vector of model names.
#' @param model_colors Optional color vector.
#' @keywords internal
resolve_model_colors <- function(model_names, model_colors = NULL) {
  n <- length(model_names)
  default_cols <- grDevices::hcl.colors(n, "Dark 3")
  names(default_cols) <- model_names

  if (is.null(model_colors)) {
    return(default_cols)
  }

  if (!is.null(names(model_colors)) && any(names(model_colors) %in% model_names)) {
    out <- default_cols
    nm <- intersect(names(model_colors), model_names)
    out[nm] <- model_colors[nm]
    return(out)
  }

  model_colors <- as.character(model_colors)
  if (length(model_colors) < n) {
    model_colors <- rep(model_colors, length.out = n)
  }
  out <- model_colors[seq_len(n)]
  names(out) <- model_names
  out
}

#' Evaluate a Kaplan-Meier survival step function on a grid
#' @param km `survfit` object.
#' @param grid Numeric grid.
#' @keywords internal
survival_step_at_grid <- function(km, grid) {
  idx <- findInterval(grid, km$time)
  out <- rep(1, length(grid))
  positive <- idx > 0
  out[positive] <- km$surv[idx[positive]]
  out
}

#' Plot density comparison using logspline
#' @param x bayesian_imputation object
#' @param dataset_id Dataset index to use for the imputed density curve (NULL for random)
#' @param n_curves Deprecated and ignored (kept for backward compatibility)
#' @param alpha Deprecated and ignored (kept for backward compatibility)
#' @param ... Additional arguments
#' @return ggplot object
#' @keywords internal
plot_density_comparison <- function(x, dataset_id = NULL, n_curves = NULL, alpha = 0.3, ...) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  if (!requireNamespace("logspline", quietly = TRUE)) {
    stop("logspline package required for density comparison")
  }
  
  time_var <- x$time_col
  n_datasets <- length(x$imputed_datasets)
  if (n_datasets < 1) {
    stop("No imputed datasets available for density plot")
  }
  
  if (is.null(dataset_id)) {
    dataset_id <- sample(seq_len(n_datasets), 1)
  }
  if (!is.numeric(dataset_id) || length(dataset_id) != 1 || !is.finite(dataset_id)) {
    stop("dataset_id must be a single numeric value between 1 and ", n_datasets)
  }
  dataset_id <- as.integer(dataset_id)
  if (dataset_id < 1 || dataset_id > n_datasets) {
    stop("dataset_id must be between 1 and ", n_datasets)
  }
  
  # Original data density
  original_data <- x$original_data[x$original_data[[x$status_col]] == 1, ]  # Observed events only
  if (nrow(original_data) < 2) {
    stop("Need at least two observed events to compute density plot")
  }
  selected_dataset <- x$imputed_datasets[[dataset_id]]
  
  # Create grid for plotting (stabilised bounds)
  x_range <- range(
    c(original_data[[time_var]], selected_dataset[[time_var]]),
    finite = TRUE
  )
  x_range[1] <- max(0, x_range[1])
  x_grid <- seq(x_range[1], x_range[2], length.out = 200)
  
  density_on_grid <- function(values, grid) {
    ls_fit <- tryCatch(
      logspline::logspline(values, lbound = 0, maxknots = 6),
      error = function(e) tryCatch(
        logspline::logspline(values, lbound = 0),
        error = function(e2) NULL
      )
    )
    
    if (!is.null(ls_fit)) {
      return(logspline::dlogspline(grid, ls_fit))
    }
    
    dens <- stats::density(values, from = min(grid), to = max(grid), n = length(grid))
    stats::approx(dens$x, dens$y, xout = grid, rule = 2)$y
  }
  
  # Original density
  y_orig <- density_on_grid(original_data[[time_var]], x_grid)
  
  # Imputed density (single selected dataset)
  y_imp <- density_on_grid(selected_dataset[[time_var]], x_grid)
  
  # Create plot data
  imputed_label <- paste0("Imputed (dataset ", dataset_id, ")")
  plot_data <- data.frame(
    x = rep(x_grid, 2),
    y = c(y_orig, y_imp),
    type = rep(c("Original", imputed_label), each = length(x_grid)),
    stringsAsFactors = FALSE
  )
  
  # Create plot 
  time_unit <- x$model_info$time_unit %||% "days"
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = type)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_color_manual(values = stats::setNames(c("black", "#7a7a7a"), c("Original", imputed_label)), name = NULL) +
    ggplot2::labs(
      title = NULL,
      subtitle = NULL,
      x = paste0("Time (", time_unit, ")"),
      y = "Density"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))
  p <- style_bird_plot(
    p,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1),
    base_size = 13
  ) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 13),
      axis.title = ggplot2::element_text(size = 16)
    )
  
  return(p)
} 

# ============================================================================
# GROUP PLOTTING FUNCTIONS
# ============================================================================

#' Plot survival curves by group
#' @param x bayesian_imputation_groups object
#' @param n_curves Number of imputed curves to show per group (default: 10)
#' @param alpha Transparency for imputed curves (default: 0.3)
#' @param show_original Whether to show original KM curves (default: TRUE)
#' @param combine_groups Whether to show all groups on same plot (TRUE) or separate plots (FALSE)
#' @param ... Additional arguments
#' @keywords internal
plot_survival_curves_by_group_safe <- function(x, n_curves = 10, alpha = 0.3, show_original = TRUE,
                                               combine_groups = TRUE, ...) {
  
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("survival package required for Kaplan-Meier curves")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  
  # Check if we have groups to plot
  if (length(x$group_names) < 1) {
    stop("No successful groups to plot")
  }
  
  # Define color palette for groups
  group_colors <- resolve_group_colors(
    x$group_names,
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  )
  
  if (combine_groups) {
    # Combined plot - all groups on same plot
    return(plot_survival_curves_combined(x, n_curves, alpha, show_original, group_colors, ...))
  } else {
    # Separate plots - one plot per group
    return(plot_survival_curves_separate(x, n_curves, alpha, show_original, group_colors, ...))
  }
}

#' Plot survival curves for all groups combined
#' @param x bayesian_imputation_groups object
#' @param n_curves Number of curves per group
#' @param alpha Transparency
#' @param show_original Whether to show original curves
#' @param group_colors Color palette
#' @param ... Additional arguments
#' @keywords internal
plot_survival_curves_combined <- function(x, n_curves, alpha, show_original, group_colors, ...) {
  group_color_map <- resolve_group_colors(x$group_names, group_colors)
  
  # Collect all plot data
  all_plot_data <- data.frame()
  
  for (i in seq_along(x$group_names)) {
    group_name <- x$group_names[i]
    group_result <- x$group_results[[group_name]]
    group_color <- group_colors[i]
    
    # Get column names from this group
    time_var <- group_result$time_col
    status_var <- group_result$status_col
    
    # Calculate median survival times and survival curves from ALL imputed datasets
    all_median_times <- numeric(length(group_result$imputed_datasets))
    all_km_curves <- list()
    
    for (i in seq_along(group_result$imputed_datasets)) {
      dataset <- group_result$imputed_datasets[[i]]
      km_temp <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = dataset
      )
      
      # Store KM curve
      all_km_curves[[i]] <- data.frame(
        time = km_temp$time,
        survival = km_temp$surv,
        stringsAsFactors = FALSE
      )
      
      # Calculate median survival time
      median_idx <- which(km_temp$surv <= 0.5)[1]
      if (!is.na(median_idx)) {
        all_median_times[i] <- km_temp$time[median_idx]
      } else {
        all_median_times[i] <- max(km_temp$time)  # If survival never drops below 0.5
      }
    }
    

    
    # Store median survival statistics for annotation
    if (!exists("median_stats")) median_stats <- list()
    median_stats[[group_name]] <- list(
      median_times = all_median_times,
      mean_median = mean(all_median_times, na.rm = TRUE),
      sd_median = sd(all_median_times, na.rm = TRUE),
      n_imputations = length(all_median_times)
    )
    
    # Original Kaplan-Meier for this group (optional)
    if (show_original) {
      km_orig <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = group_result$original_data
      )
      
      orig_df <- data.frame(
        time = km_orig$time,
        survival = km_orig$surv,
        group = group_name,
        type = "KM estimation",
        dataset_id = 0,
        legend_label = "KM estimation",
        stringsAsFactors = FALSE
      )
      all_plot_data <- rbind(all_plot_data, orig_df)
    }
    
    # Sample imputed datasets for plotting (representative sample)
    n_show <- min(n_curves, length(group_result$imputed_datasets))
    selected_datasets <- sample(group_result$imputed_datasets, n_show)
    
    # Compute survival curves for imputed datasets
    for (j in 1:n_show) {
      km_imp <- survival::survfit(
        as.formula(paste("survival::Surv(", time_var, ",", status_var, ") ~ 1")),
        data = selected_datasets[[j]]
      )
      
      imp_df <- data.frame(
        time = km_imp$time,
        survival = km_imp$surv,
        group = group_name,
        type = "Bayesian imputation",
        dataset_id = j,
        legend_label = paste0("Bayesian imputation (group = ", group_name, ")"),
        stringsAsFactors = FALSE
      )
      all_plot_data <- rbind(all_plot_data, imp_df)
    }
  }
  
  # Calculate actual number of curves shown per group (count within each group)
  actual_per_group <- sapply(x$group_names, function(g) {
    length(unique(all_plot_data$dataset_id[all_plot_data$type == "Bayesian imputation" & all_plot_data$group == g]))
  })
  # Use the minimum across groups for a conservative consistent subtitle
  actual_curves_per_group <- min(actual_per_group)
  
  # Create combined plot 
  p <- ggplot2::ggplot() +
    # Imputed curves - colored lines with transparency (sample)
    ggplot2::geom_step(data = subset(all_plot_data, type == "Bayesian imputation"), 
                      ggplot2::aes(x = time, y = survival, color = legend_label, group = interaction(group, dataset_id)), 
                      alpha = alpha, linewidth = 0.6)
  
  # Add KM curves if requested
  if (show_original) {
    p <- p + ggplot2::geom_step(data = subset(all_plot_data, type == "KM estimation"), 
                      ggplot2::aes(x = time, y = survival, color = legend_label, group = group), 
                      linewidth = 0.7)
  }
  
  p <- p +
    # Add median survival lines (horizontal from y=0.5 to curve, then vertical to x-axis)
    ggplot2::geom_segment(data = data.frame(
      legend_label = paste0("Bayesian imputation (group = ", names(median_stats), ")"),
      median_time = sapply(median_stats, function(x) x$mean_median),
      x_start = 0,
      y_start = 0.5,
      x_end = sapply(median_stats, function(x) x$mean_median),
      y_end = 0.5
    ), ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = legend_label), 
    linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    ggplot2::geom_segment(data = data.frame(
      legend_label = paste0("Bayesian imputation (group = ", names(median_stats), ")"),
      median_time = sapply(median_stats, function(x) x$mean_median),
      x_start = sapply(median_stats, function(x) x$mean_median),
      y_start = 0.5,
      x_end = sapply(median_stats, function(x) x$mean_median),
      y_end = 0
    ), ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = legend_label), 
    linetype = "dotted", linewidth = 0.8, alpha = 0.7) +
    # Color scheme: group colors for imputed, black for KM
    ggplot2::scale_color_manual(
      values = c(
        setNames(unname(group_color_map[x$group_names]), paste0("Bayesian imputation (group = ", x$group_names, ")")),
        "KM estimation" = "black"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0("Bayesian Imputation (n = ", actual_curves_per_group, ") vs Kaplan-Meier by Group"),
      subtitle = NULL,
      x = paste0("Time (", x$group_results[[1]]$model_info$time_unit %||% "days", ")"),
      y = "Survival probability",
      caption = paste("Median survival times: ", 
                     paste(sapply(names(median_stats), function(g) {
                       stats <- median_stats[[g]]
                       unit <- x$group_results[[g]]$model_info$time_unit %||% "days"
                       if (stats$n_imputations <= 2) {
                         paste(g, "=", round(stats$mean_median, 1), unit)
                       } else {
                         paste(g, "=", round(stats$mean_median, 1), "+/-", round(stats$sd_median, 1), unit)
                       }
                     }), collapse = " | "))
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, alpha = 1)))
  p <- style_bird_plot(
    p,
    legend_position = c(0.98, 0.98),
    legend_justification = c(1, 1),
    base_size = 13
  ) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 13),
      axis.title = ggplot2::element_text(size = 16),
      plot.caption = ggplot2::element_text(size = 10, color = "gray60", hjust = 0)
    )
  
  return(p)
}

#' Plot survival curves for groups separately
#' @param x bayesian_imputation_groups object
#' @param n_curves Number of curves per group
#' @param alpha Transparency
#' @param show_original Whether to show original curves
#' @param group_colors Color palette
#' @param ... Additional arguments
#' @keywords internal
plot_survival_curves_separate <- function(x, n_curves, alpha, show_original, group_colors, ...) {
  
  # Create list to store plots
  plot_list <- list()
  
  for (i in seq_along(x$group_names)) {
    group_name <- x$group_names[i]
    group_result <- x$group_results[[group_name]]
    
    # Calculate actual number of curves for this group
    actual_curves <- min(n_curves, length(group_result$imputed_datasets))
    
    # Use the original single-group plotting function for each group
    plot_list[[group_name]] <- plot_survival_curves(group_result, n_curves, alpha, show_original = show_original, ...) +
      ggplot2::labs(
        title = paste("Group:", group_name),
        subtitle = paste("Solid black = Original KM | Faint red =", actual_curves, "imputed datasets")
      )
  }
  
  # Return list of plots
  return(plot_list)
} 

#' Plot individual groups
#' @param x bayesian_imputation_groups object
#' @param type Plot type
#' @param ... Additional arguments
#' @keywords internal
plot_individual_groups_safe <- function(x, type, ...) {
  # Plot each group separately for other plot types
  rendered_plots <- list()
  for (group_name in x$group_names) {
    cat("Plotting group:", group_name, "with type:", type, "\n")
    p <- plot(x$group_results[[group_name]], type = type, ...)
    if (!is.null(p)) {
      print(p)
      rendered_plots[[group_name]] <- p
    }
  }
  
  invisible(rendered_plots)
} 

#' Plot posterior distributions for censored observations by group (side-by-side)
#' @param x bayesian_imputation_groups object
#' @param dataset_id Optional completed dataset index used to mark the imputed time in each panel (random if NULL)
#' @param ... Additional arguments passed through to plot_posterior_censored
#' @keywords internal
plot_posterior_censored_groups <- function(x, dataset_id = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for arranging plots")
  }

  # Build a list of per-group posterior plots
  plot_list <- list()

  for (group_name in x$group_names) {
    group_result <- x$group_results[[group_name]]

    # Identify censored observations in this group
    censored_idx <- which(group_result$original_data[[group_result$status_col]] == 0)

    if (length(censored_idx) == 0) {
      # Create a placeholder plot indicating no censored observations
      p_empty <- ggplot2::ggplot() +
        ggplot2::geom_text(ggplot2::aes(x = 0.5, y = 0.5, label = paste("No censored observations in", group_name))) +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
        ggplot2::theme_void()
      plot_list[[group_name]] <- p_empty
      next
    }

    # Pick a random censored observation for this group
    obs_id_group <- sample(censored_idx, 1)

    # Generate the plot using the single-group helper (use title, then override it to be concise)
    p <- plot_posterior_censored(group_result, obs_id = obs_id_group, dataset_id = dataset_id, show_title = TRUE, ...)

    # Set a concise per-panel title and remove subtitle/caption to avoid clutter/overlap
    p <- p + ggplot2::labs(
      title = paste0("Group: ", group_name, " -- Obs ", obs_id_group),
      subtitle = NULL,
      caption = NULL
    )

    plot_list[[group_name]] <- p
  }

  # Arrange side by side (2 columns) or in a grid if more groups, with a common title
  n_groups <- length(plot_list)
  ncol <- min(2, n_groups)
  common_title <- paste0("Posterior distributions by group (random censored obs)")
  arranged <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol, top = panel_title_bird(common_title))
  return(arranged)
} 

#' Plot density comparison by group (side-by-side panels)
#' @param x bayesian_imputation_groups object
#' @param n_curves Number of imputed datasets to include (NULL for all)
#' @param alpha Transparency for CI ribbon in each panel
#' @param ... Additional arguments passed to plot_density_comparison
#' @keywords internal
plot_density_comparison_groups <- function(x, n_curves = NULL, alpha = 0.3, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("gridExtra package required for arranging plots")
  }

  plot_list <- list()

  for (group_name in x$group_names) {
    grp <- x$group_results[[group_name]]
    p <- plot_density_comparison(grp, n_curves = n_curves, alpha = alpha, ...)
    # Make the panel title concise with group name
    # Minimal title (group only) and dynamic unit handled in child
    p <- p + ggplot2::labs(title = paste0("Group: ", group_name), subtitle = NULL)
    plot_list[[group_name]] <- p
  }

  n_groups <- length(plot_list)
  ncol <- min(2, n_groups)
  common_title <- NULL
  arranged <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol, top = common_title)
  return(arranged)
} 

mcmc_trace_local <- function(draws_array) {
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    bayesplot::mcmc_trace(draws_array)
  } else {
    stop("Package 'bayesplot' is required for trace plots. Install and load it to use this feature.")
  }
}

mcmc_pairs_local <- function(draws_array) {
  if (requireNamespace("bayesplot", quietly = TRUE)) {
    bayesplot::mcmc_pairs(draws_array)
  } else {
    stop("Package 'bayesplot' is required for pairs plots. Install and load it to use this feature.")
  }
} 
