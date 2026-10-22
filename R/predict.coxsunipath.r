#' Predict from a coxsunipath fit
#'
#' Evaluate fitted coefficients at new predictor values.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used.
#'   Interpolation uses the separate \code{lambdamat} vector for each path. No
#'   intercept column is added for Cox predictions.
#'
#' @param object Fitted object inheriting from \code{coxsunipath}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type \code{"link"} returns the Cox linear predictor; \code{"response"}
#'   returns its exponential (relative risk), without a baseline hazard.
#' @param alpha Optional values selecting fitted \code{alpha} paths, or fitted
#'   \code{lamPos} paths when \code{object$use_alpha} is false. Only exact
#'   matches are selected; if none match, all paths are retained.
#' @param ... Additional arguments, currently unused.
#'
#' @return A named list of numeric matrices, one per selected second-coordinate
#'   value. Rows correspond to new observations and columns to requested lambda
#'   values.
#' @method predict coxsunipath
#' @export
predict.coxsunipath <- function(
    object, newx, s = NULL,
    type = c("link", "response"),
    alpha,
    ...) {
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


    nfit <- lapply(seq_along(nbeta), function(x) {
        as.matrix(newx %*% nbeta[[x]])
    })

    nfit <- lapply(nfit, function(x) {
        switch(type,
            link = x,
            response = exp(x)
        )
    })
    names(nfit) <- name
    nfit
}
