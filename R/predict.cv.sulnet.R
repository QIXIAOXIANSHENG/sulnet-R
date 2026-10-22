#' Predict at cross-validated penalties
#'
#' Use the penalty selection stored by cross-validation to query the full-data
#'   model.
#'
#' @details The selected lambda values are evaluated using the full-data fit.
#'
#' @param object Cross-validation object of class \code{cv.sulnet}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s \code{"lambda.1se"} (default), \code{"lambda.min"}, or numeric
#'   lambda values.
#' @param ... Additional arguments passed to the fitted model's predict
#'   method.
#'
#' @return The fitted model method's result. Prediction matrices have new
#'   observations in rows and selected lambda values in columns.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' cvfit <- cv.sulnet(x, y, nfolds = 5)
#' newx <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' predict(cvfit, newx)
#' @method predict cv.sulnet
#' @export
predict.cv.sulnet <- function(object, newx, s = c("lambda.1se", "lambda.min"),
                              ...) {
  if (is.numeric(s)) {
    lambda <- s
  } else if (is.character(s)) {
    s <- match.arg(s)
    lambda <- object[[s]]
  } else {
    stop("Invalid form for s")
  }
  predict(object$sulnet.fit, newx, s = lambda, ...)
}
