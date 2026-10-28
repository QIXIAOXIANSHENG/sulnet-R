#' Extract coefficients from a sulnet model
#'
#' Delegate coefficient extraction to the next path class in the fitted object.
#'
#' @param object Object inheriting from \code{sulnet}.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type Coefficient output requested from the next method:
#'   \code{"coefficients"} or \code{"nonzero"}.
#' @param ... Additional arguments passed to the next coefficient method.
#'
#' @return The next method's result: a coefficient matrix or a list of
#'   matrices, or nonzero-index lists. Intercepts, when present, occupy the
#'   first row. Coefficient scale and grid selection follow the path method.
#' @importFrom stats coef
#' @method coef sulnet
#' @export
coef.sulnet <- function(object, s = NULL,
                        type = c("coefficients", "nonzero"),
                        ...) {
  NextMethod("coef")
}
