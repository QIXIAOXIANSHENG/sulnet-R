#' Extract coefficients from a coxsunipath fit
#'
#' Extract stored coefficients on the input-predictor scale, optionally
#'   interpolating along the lambda path.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used.
#'   Interpolation uses the separate \code{lambdamat} vector for each path.
#'
#' @param object Fitted object inheriting from \code{coxsunipath}.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints. Default is \code{NULL}.
#' @param type \code{"coefficients"} returns coefficients; \code{"nonzero"}
#'   returns one-based indices of nonzero predictor coefficients, excluding the
#'   intercept.
#' @param alpha Optional values selecting fitted \code{alpha} paths, or fitted
#'   \code{lamPos} paths when \code{object$use_alpha} is false. Only exact
#'   matches are selected; if none match, all paths are retained.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list of coefficient matrices, one per selected second-coordinate
#'   value (\code{alpha}). Each matrix has p rows and one column per requested lambda;
#'   Cox models have no intercept. With \code{type = "nonzero"}, returns lists of
#'   predictor indices by path step nested within second-coordinate lists.
#' @note Nonzero-index results are computed with apply and may simplify to vectors
#'   or arrays for small or equally sized selections.
#' @method coef coxsunipath
#' @export
coef.coxsunipath <- function(object, s = NULL,
                             type = c("coefficients", "nonzero"),
                             alpha, ...) {
    type <- match.arg(type)
    nbeta <- object$beta
    name <- names(nbeta)
    if (!is.null(s)) {
        nbeta <- lapply(seq_along(nbeta), function(x) {
            lambda <- object$lambdamat[[x]]
            lamlist <- lambda.interp(lambda, s)
            xuse <- nbeta[[x]]
            vnames <- dimnames(xuse)[[1]]
            dimnames(xuse) <- list(NULL, NULL)

            xuse <- xuse[, lamlist$left, drop = FALSE] %*%
                Diagonal(x = lamlist$frac) +
                xuse[, lamlist$right, drop = FALSE] %*%
                Diagonal(x = 1 - lamlist$frac)
            dimnames(xuse) <- list(vnames, paste(seq(along = s)))
            xuse
        })
        names(nbeta) <- name
    }
    if (!missing(alpha)) {
        nbeta_2 <- list()
        if (object$use_alpha) {
            for (a in alpha) {
                if (a %in% object$alpha) {
                    nbeta_2[[paste0("alpha_", a)]] <- nbeta[[paste0("alpha_", a)]]
                }
            }
        } else {
            for (a in alpha) {
                if (a %in% object$lamPos) {
                    nbeta_2[[paste0("lamPos_", a)]] <- nbeta[[paste0("lamPos_", a)]]
                }
            }
        }
        if (length(nbeta_2) != 0) {
            nbeta <- nbeta_2
        }
    }
    if (type == "coefficients") {
        return(nbeta)
    }
    if (type == "nonzero") {
        new_nbeta <- lapply(nbeta, function(x) {
            nonzero(x[, , drop = FALSE], bystep = TRUE)
        })
        return(new_nbeta)
    }
}
