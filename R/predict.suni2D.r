#' Predict on a two-dimensional penalty grid
#'
#' Delegate prediction to the underlying grid path method.
#'
#' @details Coefficients are interpolated linearly in each grid coordinate, with
#'   out-of-range requests clamped to the fitted endpoints. The requested
#'   coordinates form a Cartesian grid. Original-predictor coefficients are
#'   \eqn{\gamma_j = \beta_j\theta_j}; their intercept is \eqn{\gamma_0 =
#'   \theta_0 + \sum_j \beta_{0j}\theta_j}.
#'
#' @param object Object inheriting from \code{suni2D}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param penalty1 Values of the first grid coordinate: zeta, lambda, or sum
#'   according to \code{object$method}. \code{NULL} uses the fitted grid.
#' @param penalty2 Values of the second grid coordinate: tau, alpha, or rate
#'   according to \code{object$method}. \code{NULL} uses the fitted grid.
#' @param type Type forwarded to the next method. Currently only \code{links}
#' is supported.
#' @param ... Additional arguments, currently unused by the path method.
#'
#' @return A named list of numeric matrices, one per selected second-coordinate
#'   value, with new observations in rows and first-coordinate values in
#'   columns.
#'
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' fit <- suni2D(x, y)
#' newx <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' predict(fit, newx)
#'
#' @method predict suni2D
#' @export
predict.suni2D <- function(
  object, newx, penalty1 = NULL,
  penalty2 = NULL, type = c("link"),
  ...
) {
  NextMethod("predict")
}
