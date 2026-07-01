#' Unified impute() wrapper
#'
#' A single entry point that routes to either parametric (`bayesian_impute`)
#' or nonparametric LDDP (`bayes_np_impute`) engines based on arguments.
#'
#' - If `model = "nonparametric"` or `distribution` is one of
#'   `c("nonparametric","np","lddp","nonparametric_lddp")`, uses the
#'   nonparametric engine.
#' - Otherwise, uses the parametric engine with the chosen `distribution`.
#'
#' @param data Data frame or `survival_data` prepared by `prepare_survival_data()`
#' @param time Time column name. Required unless `data` is a `survival_data` object.
#' @param status Status column name (1=event, 0=censored). Required unless `data` is a `survival_data` object.
#' @param groups Optional group variable name for group analyses
#' @param model One of c("auto","parametric","nonparametric"). Default "auto".
#' @param distribution For parametric: one of c("weibull","exponential","lognormal").
#'   Nonparametric can be requested via distribution tokens in
#'   `c("nonparametric","np","lddp","nonparametric_lddp")`.
#' @param n_imputations Number of completed datasets to generate
#' @param time_unit Display unit for times (e.g., "days","months","years")
#' @param priors Parametric priors list (see `get_default_priors()` / `get_adaptive_priors()`)
#' @param prior Nonparametric prior list (LDDP hyperparameters)
#' @param mcmc_options Parametric MCMC options list
#' @param mcmc Nonparametric MCMC options list
#' @param verbose Logical; print progress
#' @param ... Passed through to the underlying engine
#'
#' @return A `bayesian_imputation` or `bayesian_imputation_groups` object
#' @seealso [bayesian_impute()], [bayes_np_impute()], [complete()]
#' @export
impute <- function(data,
                   time = NULL,
                   status = NULL,
                   groups = NULL,
                   model = c("auto","parametric","nonparametric"),
                   distribution = c("weibull","exponential","lognormal",
                                    "nonparametric","np","lddp","nonparametric_lddp"),
                   n_imputations = 10,
                   time_unit = "days",
                   priors = NULL,          # parametric
                   prior = NULL,           # nonparametric
                   mcmc_options = NULL,    # parametric
                   mcmc = NULL,            # nonparametric
                   verbose = TRUE,
                   ...) {

  model <- match.arg(model)
  distribution <- match.arg(distribution)

  # Capture additional args
  dots <- list(...)

  # Prepare/normalise data
  prepared_here <- FALSE
  if (!inherits(data, "survival_data")) {
    if (is.null(time) || is.null(status)) {
      stop("Please specify both 'time' and 'status'; auto-detection is not yet supported.")
    }
    data <- prepare_survival_data(data, time = time, status = status,
                                  time_unit = time_unit, verbose = FALSE)
    prepared_here <- TRUE
  }

  time_col <- attr(data, "time_col")
  status_col <- attr(data, "status_col")

  if (is.null(time_col) || is.null(status_col)) {
    stop("Please specify both 'time' and 'status'; auto-detection is not yet supported.")
  }

  if (verbose && prepared_here) {
    cat(sprintf("Using columns: time='%s', status='%s'\n", time_col, status_col))
  }

  # Decide engine
  np_tokens <- c("nonparametric","np","lddp","nonparametric_lddp")
  is_np <- (model == "nonparametric") || (distribution %in% np_tokens)

  # Warn on mismatched options
  if (is_np && (!is.null(priors) || !is.null(mcmc_options))) {
    warning("Ignoring parametric options (priors/mcmc_options) for nonparametric model")
  }
  if (!is_np && (!is.null(prior) || !is.null(mcmc))) {
    warning("Ignoring nonparametric options (prior/mcmc) for parametric model")
  }

  # Internal helper
  invoke_matching <- function(fun, arglist) {
    formals_names <- tryCatch(names(formals(fun)), error = function(e) NULL)
    if (!is.null(formals_names)) {
      arglist <- arglist[intersect(names(arglist), formals_names)]
    }
    do.call(fun, arglist)
  }

  if (is_np) {
    base_args <- list(
      data = data,
      time_col = time_col,
      status_col = status_col,
      groups = groups,
      n_imputations = n_imputations,
      prior = prior,
      verbose = verbose,
      time_unit = time_unit
    )
    if (!is.null(mcmc)) {
      base_args$mcmc <- mcmc
    }
    arglist <- c(base_args, dots)
    return(invoke_matching(bayes_np_impute, arglist))
  } else {
    # For parametric models, use the specified distribution
    arglist <- c(list(
      data = data,
      time_col = time_col,
      status_col = status_col,
      groups = groups,
      n_imputations = n_imputations,
      distribution = distribution,
      priors = priors,
      mcmc_options = mcmc_options,
      verbose = verbose,
      time_unit = time_unit
    ), dots)
    return(invoke_matching(bayesian_impute, arglist))
  }
}

