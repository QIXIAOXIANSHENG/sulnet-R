#' Predict at cross-validated penalties
#'
#' Use the penalty selection stored by cross-validation to query the full-data
#'   model.
#'
#' @details Both coordinates of the selected penalty pair are evaluated using
#'   the full-data fit. This method currently looks for \code{lambda.1se} or
#'   \code{lambda.min} inside the selected pair, whereas \code{cv.sulnet2D}
#'   stores its first coordinate as \code{lambda}.
#'
#' @param object Cross-validation object of class \code{cv.sulnet2D}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s \code{"cv.1se"} (default) or \code{"cv.min"}, selecting the stored
#'   penalty pair. Numeric input is not accepted.
#' @param ... Additional arguments passed to the fitted model's predict
#'   method.
#'
#' @return The fitted model method's result. Prediction matrices have new
#'   observations in rows and selected lambda values in columns.
#'
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' cvfit <- cv.sulnet2D(x, y, nfolds = 5)
#' newx <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' predict(cvfit, newx)
#'
#' @method predict cv.sulnet2D
#' @export
predict.cv.sulnet2D <- function(object, newx, s = c("cv.1se", "cv.min"),
                                ...) {
  if (is.character(s)) {
    s <- match.arg(s)
    if (s == "cv.1se") {
      lambda <- object$cv.1se$lambda.1se
      alpha <- object$cv.1se[[2]]
    } else {
      lambda <- object$cv.min$lambda.min
      alpha <- object$cv.min[[2]]
    }
  } else {
    stop("Invalid form for s")
  }
  predict(object$sulnet2D.fit, newx, s = lambda, alpha = alpha, ...)
}
