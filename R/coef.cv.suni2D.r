#' Extract coefficients at cross-validated penalties for a
#' \code{cv.suni2D} object.
#'
#' Use the penalty selection stored by cross-validation to query the full-data
#'   model.
#'
#' @details Both coordinates of the selected penalty pair are evaluated using
#'   the full-data fit.
#'
#' @param object Cross-validation object of class \code{cv.suni2D}.
#' @param s \code{"cv.1se"} (default) or \code{"cv.min"}, selecting the stored
#'   penalty pair. Numeric input is not accepted.
#' @param ... Additional arguments passed to the fitted model's coef
#'   method.
#'
#' @return A list with the fitted model' coefficients. Names of the list is the
#'   corresponding second penalty of the \code{cv.suni2D} object. Coefficient
#'   matrices include an intercept first when the model has one; nonzero output
#'   excludes intercepts.
#' @method coef cv.suni2D
#' @export
coef.cv.suni2D <- function(object, s = c("cv.1se", "cv.min"), ...) {
  if (is.character(s)) {
    s <- match.arg(s)
    la <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  coef(object$suni2D.fit, penalty1 = la[[1]], penalty2 = la[[2]], ...)
}
