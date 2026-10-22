#' Predict from a sulnet model
#'
#' Delegate prediction to the next path class in the fitted object.
#'
#' @param object Object inheriting from \code{sulnet}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type Prediction type forwarded to the next method. Available choices
#'   depend on the path class: Gaussian paths accept \code{"link"}, logistic
#'   paths also accept \code{"class"}, and Cox paths accept \code{"link"} or
#'   \code{"response"}.
#' @param ... Additional arguments passed to the next prediction method.
#'
#' @return The next prediction method's result, usually a matrix or a
#'   list of matrices with observations in rows and selected path values in
#'   columns.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' fit <- sulnet(x, y)
#' newx <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' predict(fit, newx)
#'
#' @method predict sulnet
#' @export
predict.sulnet <- function(object, newx, s = NULL, type = c("class", "link"),
                           ...) {
  NextMethod("predict")
}

#' Predict from a fitted model
#'
#' Dispatch to the method associated with the class of the fitted object.
#'
#' @param object Fitted model object.
#' @param ... Additional arguments passed to the selected method.
#'
#' @return The value returned by the selected method; see the corresponding
#'   class method for coefficient scales and dimensions.
#' @export
predict <- function(object, ...) UseMethod("predict")
