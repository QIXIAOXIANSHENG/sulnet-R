#####################################################################
# This function is adapted/modified based on the plot function
# from the `glmnet` package:
#
# Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
# Regularization Paths for Generalized Linear Models via Coordinate Descent.
#   Journal of Statistical Software, 33(1), 1-22.
#   URL: https://www.jstatsoft.org/v33/i01/.

#' Plot coefficient paths for a sulnet2D fit
#'
#' Draw slope trajectories, excluding the intercept and predictors that are zero
#'   over the entire plotted path.
#'
#' @details The top axis shows the stored nonzero-coefficient counts along the
#'   path.
#'
#' @param x A fitted \code{sulnet2D} object.
#' @param xvar Horizontal coordinate: \code{"norm"} for the L1 norm of plotted
#'   slopes, or \code{"lambda"} for log lambda.
#' @param yvar \code{"final_coef"} plots stored final \code{beta};
#'   \code{"middle_coef"} plots \code{fbeta} on the transformed design, when
#'   available.
#' @param color When FALSE, explicitly use the package's rainbow palette;
#'   when TRUE, leave colors to \code{matplot} and arguments in \code{...}.
#' @param label Whether to label curves with one-based predictor indices.
#' @param ... Additional graphical arguments passed to \code{matplot}.
#'
#' @return Creates one panel per second-coordinate value, restores the previous
#'   graphics settings, and invisibly returns the settings replaced by the final
#'   \code{par} call.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' cvfit <- cv.sulnet2D(x, y, nfolds = 5)
#' plot(cvfit$sulnet2D.fit, xvar = "lambda", yvar = "final_coef")
#'
#' @importFrom graphics matplot axis text
#' @importFrom grDevices gray.colors rainbow
#' @importFrom stats approx
#' @method plot sulnet2D
#' @export
plot.sulnet2D <- function(
  x, xvar = c("norm", "lambda"),
  yvar = c("final_coef", "middle_coef"),
  color = FALSE,
  label = FALSE, ...
) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  yvar <- match.arg(yvar)
  if (yvar == "final_coef") {
    betaall <- x$beta
    ylab <- "Final Coefficients"
  } else {
    betaall <- x$fbeta
    ylab <- "Middle Coefficients"
  }
  lambda <- x$lambda
  dfall <- x$df
  xvar <- match.arg(xvar)
  if (x$use_alpha) {
    plotname <- paste0("alpha: ", x$alpha)
  } else {
    plotname <- paste0("lamPos: ", x$lamPos)
  }
  ncanvas <- length(plotname)

  ## beta should be in 'dgCMatrix' format
  par(mfrow = canvas_structure(ncanvas), mar = c(3, 3, 2, 1))
  for (i in seq_along(betaall)) {
    beta <- betaall[[i]]
    df <- dfall[[i]]
    which <- nonzero(beta)
    beta <- as.matrix(beta[which, ])
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
    ylab <- ylab
    dotlist <- list(...)
    type <- dotlist$type


    if (is.null(type)) {
      if (color == FALSE) {
        matplot(index, t(beta),
          lty = 1, xlab = xlab,
          ylab = ylab, type = "l", pch = 500,
          col = rainbow(12, start = 0.7, end = 0.95),
          main = plotname[i],
          ...
        )
      } else {
        matplot(index, t(beta),
          lty = 1,
          xlab = xlab, ylab = ylab,
          type = "l", pch = 500,
          main = plotname[i], ...
        )
      }
    } else {
      matplot(index, t(beta),
        lty = 1, xlab = xlab, ylab = ylab,
        main = plotname[i], ...
      )
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
  par(op)
}
