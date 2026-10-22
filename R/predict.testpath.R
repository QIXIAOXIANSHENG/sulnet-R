#' Predict from a testpath fit
#'
#' Evaluate fitted coefficients at new predictor values.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used. The
#'   stored intercept is added to each prediction.
#'
#' @param object Fitted object inheriting from \code{testpath}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type Prediction type; only \code{"link"} is accepted.
#' @param ... Additional arguments, currently unused.
#'
#' @return A numeric matrix. Rows correspond to new observations and columns to
#'   requested lambda values.
#' @method predict testpath
#' @export
predict.testpath <- function(object, newx, s = NULL, type = c("link"), ...) {
  type <- match.arg(type)
  b0 <- t(as.matrix(object$b0))
  rownames(b0) <- "(Intercept)"
  nbeta <- rbind2(b0, object$beta)
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)
    nbeta <- nbeta[, lamlist$left, drop = FALSE] %*%
      Diagonal(x = lamlist$frac) +
      nbeta[, lamlist$right, drop = FALSE] %*%
      Diagonal(x = 1-lamlist$frac)
    dimnames(nbeta) <- list(vnames, paste(seq(along = s)))
  }
  nfit <- as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta)
  nfit
}
