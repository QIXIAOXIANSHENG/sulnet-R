#' Predict on a two-dimensional penalty grid
#'
#' Evaluate reconstructed original-predictor coefficients at new observations.
#'
#' @details Coefficients are interpolated linearly in each grid coordinate, with
#'   out-of-range requests clamped to the fitted endpoints. The requested
#'   coordinates form a Cartesian grid. Original-predictor coefficients are
#'   \eqn{\gamma_j = \beta_j\theta_j}; their intercept is \eqn{\gamma_0 =
#'   \theta_0 + \sum_j \beta_{0j}\theta_j}.
#'
#' @param object Object inheriting from \code{sunipath}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param penalty1 Values of the first grid coordinate: zeta, lambda, or sum
#'   according to \code{object$use}. \code{NULL} uses the fitted grid.
#' @param penalty2 Values of the second grid coordinate: tau, alpha, or rate
#'   according to \code{object$use}. \code{NULL} uses the fitted grid.
#' @param type Only \code{"link"} is accepted.
#' @param ... Additional arguments, currently unused by the path method.
#'
#' @return A named list of numeric matrices, one per selected second-coordinate
#'   value, with new observations in rows and first-coordinate values in
#'   columns.
#' @method predict sunipath
#' @export
predict.sunipath <- function(
  object, newx,
  penalty1 = NULL, penalty2 = NULL, type = c("link"), ...
) {
  type <- match.arg(type)
  gammas <- coef(object, penalty1, penalty2)
  nnewx <- as.matrix(cbind2(1, newx))
  nfit <- lapply(seq_along(gammas), function(x) {
    as.matrix(nnewx %*% gammas[[x]])
  })
  names(nfit) <- names(gammas)
  nfit
}
