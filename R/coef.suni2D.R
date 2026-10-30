#' Extract coefficients from a suni2D model
#'
#' Extract coefficients from a fitted two-dimensional univariate-guided model.
#' Gaussian fits delegate coefficient extraction to \code{coef.sunipath} and
#' binomial fits to \code{coef.sunilogpath}; both work the same way.
#'
#' @param object A fitted object inheriting from \code{suni2D}.
#' @param penalty1 Values of the first grid coordinate at which coefficients
#'   are requested. \code{NULL} uses the fitted grid.
#' @param penalty2 Values of the second grid coordinate at which coefficients
#'   are requested. \code{NULL} uses the fitted grid.
#' @param type Coefficient output requested: \code{"coefficients"},
#'   \code{"nonzero"}, or \code{"inner"}.
#' @param ... Additional arguments passed to \code{coef.sunipath} or
#'   \code{coef.sunilogpath}.
#'
#' @return The value returned by \code{coef.sunipath} (Gaussian) or
#'   \code{coef.sunilogpath} (binomial): for \code{type = "coefficients"}, a
#'   list with one sparse matrix per second-coordinate value, with the
#'   intercept in the first row and the requested first-coordinate values in
#'   columns.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' fit <- suni2D(x, y)
#' penalty1 <- fit$zeta[c(5, 90)] * 0.5
#' penalty2 <- fit$tau[1:5] * 0.5
#' coef(fit, penalty1, penalty2)
#' @method coef suni2D
#' @export
coef.suni2D <- function(
  object, penalty1 = NULL, penalty2 = NULL,
  type = c("coefficients", "nonzero", "inner"), ...
) {
  NextMethod("coef")
}
