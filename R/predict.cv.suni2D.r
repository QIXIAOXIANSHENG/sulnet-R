#' Predict at cross-validated penalties for a \code{cv.suni2D} fit
#'
#' Use the penalty selection stored by cross-validation to predict from the
#' full-data model.
#'
#' @details Both coordinates of the selected penalty pair are evaluated using
#'   the full-data fit.
#'
#' @param object Cross-validation object of class \code{cv.suni2D}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s \code{"cv.1se"} (default) or \code{"cv.min"}, selecting the stored
#'   penalty pair. Numeric input is not accepted.
#' @param ... Additional arguments passed to the fitted model's predict method.
#'
#' @return A named list containing a numeric prediction matrix for the selected
#'   second-coordinate value. Rows correspond to new observations and the
#'   column corresponds to the selected first-coordinate value.
#'
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' cvfit <- cv.suni2D(x, y, nfolds = 5)
#' newx <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' predict(cvfit, newx)
#'
#' @method predict cv.suni2D
#' @export
predict.cv.suni2D <- function(
  object, newx,
  s = c("cv.1se", "cv.min"),
  ...
) {
  if (is.character(s)) {
    s <- match.arg(s)
    pen <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  predict(object$suni2D.fit, newx,
    penalty1 = pen[[1]], penalty2 = pen[[2]], ...
  )
}
