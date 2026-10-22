#' Predict from a logsunipath fit
#'
#' Evaluate fitted coefficients at new predictor values.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used. The
#'   stored intercept is added to each prediction.
#'
#' @param object Fitted object inheriting from \code{logsunipath}.
#' @param newx Predictor matrix for new observations, with the same columns and
#'   order as the training predictors and no intercept column.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type \code{"link"} returns the linear predictor; \code{"class"}
#'   returns 1 for positive values and -1 otherwise. Probabilities are not
#'   returned.
#' @param alpha Optional values selecting fitted \code{alpha} paths, or fitted
#'   \code{lamPos} paths when \code{object$use_alpha} is false. Only exact
#'   matches are selected; if none match, all paths are retained.
#' @param ... Additional arguments, currently unused.
#'
#' @return A named list of numeric matrices, one per selected second-coordinate
#'   value. Rows correspond to new observations and columns to requested lambda
#'   values.
#' @method predict logsunipath
#' @export
predict.logsunipath <- function(
    object, newx, s = NULL,
    type = c("class", "link"),
    alpha,
    ...) {
    type <- match.arg(type)
    b0 <- lapply(object$b0, function(x) {
        b <- t(as.matrix(x))
        rownames(b) <- "(Intercept)"
        b
    })
    nbeta <- lapply(seq_along(b0), function(x) {
        rbind2(b0[[x]], object$beta[[x]])
    })
    names(nbeta) <- names(b0)
    if (!is.null(s)) {
        lambda <- object$lambda
        lamlist <- lambda.interp(lambda, s)

        nbeta <- lapply(nbeta, function(x) {
            vnames <- dimnames(x)[[1]]
            dimnames(x) <- list(NULL, NULL)

            x <- x[, lamlist$left, drop = FALSE] %*%
                Diagonal(x = lamlist$frac) +
                x[, lamlist$right, drop = FALSE] %*%
                Diagonal(x = 1 - lamlist$frac)
            dimnames(x) <- list(vnames, paste(seq(along = s)))
            x
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


    nnewx <- as.matrix(cbind2(1, newx))
    nfit <- lapply(seq_along(nbeta), function(x) {
        as.matrix(nnewx %*% nbeta[[x]])
    })

    nfit <- lapply(nfit, function(x) {
        switch(type,
            link = x,
            class = ifelse(x > 0, 1, -1)
        )
    })
    names(nfit) <- names(nbeta)
    nfit
}
