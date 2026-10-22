#' Extract coefficients from a lspath fit
#'
#' Extract stored coefficients on the input-predictor scale, optionally
#'   interpolating along the lambda path.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used.
#'
#' @param object Fitted object inheriting from \code{lspath}.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type \code{"coefficients"} returns coefficients; \code{"nonzero"}
#'   returns one-based indices of nonzero predictor coefficients, excluding the
#'   intercept.
#' @param ... Additional arguments, currently unused.
#'
#' @return A coefficient matrix. Each matrix has p + 1 rows, with the intercept
#'   first, and one column per requested lambda. With \code{type = "nonzero"},
#'   returns lists of predictor indices by path step.
#' @importFrom methods rbind2
#' @note Nonzero-index results are computed with apply and may simplify to vectors
#'   or arrays for small or equally sized selections.
#' @method coef lspath
#' @export
coef.lspath <- function(object, s = NULL, type = c("coefficients", "nonzero"),
                        ...) {
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
  if (type == "coefficients")
    return(nbeta)
  if (type == "nonzero")
    return(nonzero(nbeta[-1, , drop = FALSE], bystep = TRUE))
}
