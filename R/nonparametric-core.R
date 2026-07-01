#' Non-parametric DP Survival (Wrapper)
#'
#' Internal wrapper for the vendored `LDDPsurvival()` function.
#' All arguments, return value and behaviour are identical to the original
#' function from **DPpackage**; we only forward the call.
#'
#' @param ... All arguments passed straight to `LDDPsurvival()`.
#'
#' @return The object returned by the internal `LDDPsurvival()` call (class
#'   `LDDPsurvival`).
#'
#' @seealso Internally calls [LDDPsurvival.default()]. A higher-level helper
#'   [bayes_np_impute()] provides user-friendly data handling.
#' @keywords internal
#' @examples
#' # Not run: see bayes_np_impute() for a full example.
#' # dp_np_LDDPsurvival(ymat ~ 1, zpred = matrix(1,1,1), ...)
#'
dp_np_LDDPsurvival <- function(...) {
  LDDPsurvival(...)
} 
