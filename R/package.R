#' @useDynLib bird, .registration = TRUE
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang .data
#' @importFrom stats median quantile as.formula density
#' @importFrom utils head
#' @importFrom graphics plot
#' @importFrom survival survreg Surv
#' @importFrom stats aggregate coefficients model.matrix model.response na.fail optimize plnorm rexp rlnorm runif rweibull sd setNames ts uniroot vcov
#' @importFrom graphics axis hist layout lines par polygon segments
#' @importFrom methods is
#' @import methods
#' @importFrom gridExtra grid.arrange
#' @importFrom haven write_dta write_sav
#' @importFrom scales percent
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
## usethis namespace: end

# Global variables for ggplot2 NSE
utils::globalVariables(
  c(
    "time", "survival", "type", "dataset_id", "density", "imputed_time",
    "dataset", "group", "y", "ymin", "ymax", "imputed_time_draw",
    "xintercept", "line_type", "lower", "upper", "x_start", "y_start",
    "x_end", "y_end", "legend_label"
  )
)

NULL

# Package environment for storing compiled Stan models
.bird_env <- new.env(parent = emptyenv())

#' bird: Bayesian Imputation for Right-Censored Data
#'
#' @description
#' The bird package implements the Bayesian imputation methodology 
#' from Moghaddam et al. (2022) for handling right-censored survival data. 
#' The package treats censored observations as missing data and uses Bayesian 
#' methods to generate posterior distributions for each censored observation,
#' enabling complete dataset generation.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{impute}}}{Unified function for Bayesian imputation (parametric and nonparametric)}
#'   \item{\code{\link{bayesian_impute}}}{Parametric Bayesian imputation of censored survival data}
#'   \item{\code{\link{bayes_np_impute}}}{Nonparametric Bayesian imputation using Dirichlet Process}
#'   \item{\code{\link{complete}}}{Extract complete datasets from bayesian_imputation objects}
#'   \item{\code{\link{export}}}{Export completed datasets in multiple formats (CSV, RDS)}
#'   \item{\code{\link{generate_complete_datasets}}}{Generate additional complete datasets from posterior distributions}
#' }
#'
#'
#' @section S3 methods:
#' \describe{
#'   \item{\code{\link{print.bayesian_imputation}}}{Print method for bayesian_imputation objects}
#'   \item{\code{\link{plot.bayesian_imputation}}}{Plot method with multiple visualisation types}
#'   \item{\code{\link{complete.bayesian_imputation}}}{Extract completed datasets}
#' }
#'
#' @section Key features:
#' \itemize{
#'   \item Bayesian imputation for right-censored survival data using Weibull, exponential, or lognormal distributions
#'   \item Posterior-based imputation: one imputation per censored observation per MCMC iteration
#'   \item Unlimited dataset generation by sampling from stored posterior distributions
#'   \item Comprehensive visualisations: survival curves, density plots, posterior plots, trace plots, pairs plots
#'   \item MCMC convergence diagnostics and validation
#'   \item Simulation framework with realistic censoring patterns (exponential and uniform)
#'   \item Posterior diagnostics for censored observations and posterior draws
#' }
#' 
#' @section Supported distributions:
#' \itemize{
#'   \item \strong{Weibull}: Flexible hazard shapes (increasing, decreasing, constant)
#'   \item \strong{Exponential}: Constant hazard (special case of Weibull)
#'   \item \strong{Lognormal}: Alternative flexible distribution for survival data
#' }
#'
#' 
#' @section Methodology:
#' The implementation generates one imputation per censored observation per MCMC iteration, 
#' creating a full posterior distribution for each censored observation. This approach:
#' \enumerate{
#'   \item Fits a Bayesian survival model using MCMC
#'   \item Generates one imputation per censored observation per MCMC iteration
#'   \item Stores the full posterior distribution for each censored observation
#'   \item Generates complete datasets by sampling from these posteriors
#'   \item Provides posterior diagnostics for imputed datasets and model fits
#' }
#'
#' @references
#' Moghaddam, S., Newell, J., & Hinde, J. (2022). A Bayesian Approach for 
#' Imputation of Censored Survival Data. Stats, 5(1), 89-107. 
#' doi:10.3390/stats5010006
#'
#' @keywords internal
"_PACKAGE" 
