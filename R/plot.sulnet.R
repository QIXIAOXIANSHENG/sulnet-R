#####################################################################
# This function is adapted/modified based on the plot function
# from the `glmnet` package:
#
# Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
# Regularization Paths for Generalized Linear Models via Coordinate Descent.
#   Journal of Statistical Software, 33(1), 1-22.
#   URL: https://www.jstatsoft.org/v33/i01/.

#' Plot coefficient paths for a sulnet fit
#'
#' Draw slope trajectories, excluding the intercept and predictors that are zero
#'   over the entire plotted path.
#'
#' @details The top axis shows the stored nonzero-coefficient counts along the
#'   path.
#'
#' @param x A fitted \code{sulnet} object.
#' @param xvar Horizontal coordinate: \code{"norm"} for the L1 norm of plotted
#'   slopes, or \code{"lambda"} for log lambda. Default is \code{"lambda"}.
#' @param color When FALSE, explicitly use the package's rainbow palette;
#'   when TRUE, leave colors to \code{matplot} and arguments in \code{...}.
#' @param label Whether to label curves with one-based predictor indices.
#' @param ... Additional graphical arguments passed to \code{matplot}.
#'
#' @return Called for its graphical side effect; no fitted model is returned.
#' @importFrom graphics matplot axis text
#' @importFrom grDevices gray.colors rainbow
#' @importFrom stats approx
#' @method plot sulnet
#' @export
plot.sulnet <- function(x, xvar = c("lambda", "norm"), color = FALSE,
                        label = FALSE, ...) {
  beta <- x$beta
  lambda <- x$lambda
  df <- x$df
  xvar <- match.arg(xvar)
  ## beta should be in 'dgCMatrix' format
  which <- nonzero(beta)
  beta <- as.matrix(beta[which, ])
  xvar <- match.arg(xvar)
  switch(xvar,
    norm = {
      index <- apply(abs(beta), 2, sum)
      iname <- "L1 Norm"
    },
    lambda = {
      index <- log(lambda)
      iname <- "Log Lambda"
    }
  )
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE) {
      matplot(index, t(beta),
        lty = 1, xlab = xlab,
        ylab = ylab, type = "l", pch = 500,
        col = rainbow(12, start = 0.7, end = 0.95),
        ...
      )
    } else {
      matplot(index, t(beta),
        lty = 1,
        xlab = xlab, ylab = ylab,
        type = "l", pch = 500, ...
      )
    }
  } else {
    matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  }
  atdf <- pretty(index)
  prettydf <- trunc(approx(
    x = index, y = df, xout = atdf,
    rule = 2
  )$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}
