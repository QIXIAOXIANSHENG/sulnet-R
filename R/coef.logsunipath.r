#' Extract coefficients from a logsunipath fit
#'
#' Extract stored coefficients on the input-predictor scale, optionally
#'   interpolating along the lambda path.
#'
#' @details The stored \code{beta} contains the final coefficients, including
#'   reconstruction from the univariate stage when that stage was used.
#'
#' @param object Fitted object inheriting from \code{logsunipath}.
#' @param s Lambda values at which to evaluate the path. \code{NULL} uses all
#'   fitted values. Intermediate values are linearly interpolated and values
#'   outside the fitted range are clamped to its endpoints.
#' @param type \code{"coefficients"} returns coefficients; \code{"nonzero"}
#'   returns one-based indices of nonzero predictor coefficients, excluding the
#'   intercept.
#' @param alpha Optional values selecting fitted \code{alpha} paths, or fitted
#'   \code{lamPos} paths when \code{object$use_alpha} is false. Only exact
#'   matches are selected; if none match, all paths are retained.
#' @param ... Additional arguments, currently unused.
#'
#' @return A list of coefficient matrices, one per selected second-coordinate
#'   value. Each matrix has p + 1 rows, with the intercept first, and one column
#'   per requested lambda. With \code{type = "nonzero"}, returns lists of
#'   predictor indices by path step nested within second-coordinate lists.
#' @note Nonzero-index results are computed with apply and may simplify to vectors
#'   or arrays for small or equally sized selections.
#' @method coef logsunipath
#' @export
coef.logsunipath <- function(object, s = NULL,
                           type = c("coefficients", "nonzero"),
                           alpha, ...) {
  type <- match.arg(type)
  b0 <- lapply(object$b0, function(x){
    b <- t(as.matrix(x))
    rownames(b) <- "(Intercept)"
    b
  })
  nbeta <- lapply(seq_along(b0), function(x){
    rbind2(b0[[x]], object$beta[[x]])
  })
  names(nbeta) <- names(b0)
  if (!is.null(s)) {
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda, s)

    nbeta <- lapply(nbeta, function(x){
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
  if(!missing(alpha)){
    nbeta_2 <- list()
    if(object$use_alpha){
      for(a in alpha){
        if(a %in% object$alpha){
          nbeta_2[[paste0("alpha_", a)]] <- nbeta[[paste0("alpha_", a)]]
        }
      }
    }else{
      for(a in alpha){
        if(a %in% object$lamPos){
          nbeta_2[[paste0("lamPos_", a)]] <- nbeta[[paste0("lamPos_", a)]]
        }
      }
    }
    if(length(nbeta_2) != 0){
      nbeta <- nbeta_2
    }
  }
  if (type == "coefficients")
    return(nbeta)
  if (type == "nonzero"){
    new_nbeta <- lapply(nbeta, function(x){
      nonzero(x[-1, , drop = FALSE], bystep = TRUE)
    })
    return(new_nbeta)
  }
}
