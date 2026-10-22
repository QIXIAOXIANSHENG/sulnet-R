#' Extract coefficients at cross-validated penalties for a
#'   \code{cv.sulnet2D} object.
#'
#' Use the penalty selection stored by cross-validation to query the full-data
#'   model.
#'
#' @details Both coordinates of the selected penalty pair are evaluated using
#'   the full-data fit.
#'
#' @param object Cross-validation object of class \code{cv.sulnet2D}.
#' @param s \code{"cv.1se"} (default) or \code{"cv.min"}, selecting the stored
#'   penalty pair. Numeric input is not accepted.
#' @param ... Additional arguments passed to the fitted model's coef
#'   method.
#'
#' @return A list with the fitted model' coefficients. Names of the list is the
#'   corresponding \code{alpha} or \code{lamPos} value of \code{s}. Coefficient
#'   matrices include an intercept first when the model has one; nonzero output
#'   excludes intercepts.
#' @method coef cv.sulnet2D
#' @export
coef.cv.sulnet2D <- function(object, s = c("cv.1se", "cv.min"), ...) {
  if (is.character(s)) {
    s <- match.arg(s)
    la <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  coef(object$sulnet2D.fit, s = la[[1]], alpha = la[[2]], ...)
}
