#' Extract coefficients at cross-validated penalties for a
#'   \code{cv.sulnet} object
#'
#' Use the penalty selection stored by cross-validation to query the full-data
#'   model.
#'
#' @details The selected lambda values are evaluated using the full-data fit.
#'
#' @param object Cross-validation object of class \code{cv.sulnet}.
#' @param s \code{"lambda.1se"} (default), \code{"lambda.min"}, or numeric
#'   lambda values.
#' @param ... Additional arguments passed to the fitted model's coef
#'   method.
#'
#' @return A matrix of the fitted model' coefficients.Coefficient matrices
#'   include an intercept first when the model has one; nonzero output
#'   excludes intercepts.
#' @method coef cv.sulnet
#' @export
coef.cv.sulnet <- function(object, s = c("lambda.1se", "lambda.min"), ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  coef(object$sulnet.fit, s = lambda, ...)
}
