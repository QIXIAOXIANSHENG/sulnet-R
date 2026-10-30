#' Plot coefficient paths for a suni2D fit
#'
#' Draw slope trajectories, excluding the intercept and predictors that are zero
#'   over the entire plotted path.
#'
#' @details The top axis shows the stored nonzero-coefficient counts along the
#'   path.
#'
#' @param x A fitted \code{suni2D} object.
#' @param xvar Horizontal coordinate: \code{"norm"} for the L1 norm of plotted
#'   slopes, or \code{"penalty"} for the log first grid coordinate.
#'   Default is \code{"penalty"}.
#' @param yvar \code{"final_coef"} plots reconstructed \code{gamma};
#'   \code{"middle_coef"} plots \code{inner$theta} on the transformed design.
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
#' cvfit <- cv.suni2D(x, y, nfolds = 5)
#' plot(cvfit$suni2D.fit, xvar = "penalty", yvar = "final_coef")
#'
#' @method plot suni2D
#' @export
plot.suni2D <- function(x, xvar = c("penalty", "norm"), yvar = c("final_coef", "middle_coef"), color = FALSE, label = FALSE, ...) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  yvar <- match.arg(yvar)
  xvar <- match.arg(xvar)
  if (yvar == "final_coef") {
    gammaall <- x$gamma
    ylab <- "Final Coefficients"
  } else {
    gammaall <- x$inner$theta
    ylab <- "Middle Coefficients"
  }
  method <- x$method
  penname <- strsplit(method, "_")[[1]]
  pen1 <- x[[penname[1]]]
  pen2 <- x[[penname[2]]]
  dfall <- x$df
  plotname <- paste0(penname[2], pen2)
  ncanvas <- length(plotname)
  par(mfrow = canvas_structure(ncanvas), mar = c(3, 3, 2, 1))
  for (i in seq_along(gammaall)) {
    gamma <- gammaall[[i]]
    df <- dfall[, i]
    which <- nonzero(gamma)
    if (length(which) == 0L) {
      ## no active predictor anywhere on this path: draw a flat zero line
      gamma <- matrix(0, 1L, ncol(gamma))
    } else {
      gamma <- as.matrix(gamma[which, , drop = FALSE])
    }
    switch(xvar,
      norm = {
        index <- apply(abs(gamma), 2, sum)
        iname <- "L1 Norm"
      },
      penalty = {
        index <- log(pen1)
        iname <- paste0("Log ", penname[1])
      }
    )
    xlab <- iname
    ylab <- ylab
    dotlist <- list(...)
    type <- dotlist$type


    if (is.null(type)) {
      if (color == FALSE) {
        matplot(index, t(gamma),
          lty = 1, xlab = xlab,
          ylab = ylab, type = "l", pch = 500,
          col = rainbow(12, start = 0.7, end = 0.95),
          main = plotname[i],
          ...
        )
      } else {
        matplot(index, t(gamma),
          lty = 1,
          xlab = xlab, ylab = ylab,
          type = "l", pch = 500,
          main = plotname[i], ...
        )
      }
    } else {
      matplot(index, t(gamma),
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
    if (label && length(which) > 0L) {
      nnz <- length(which)
      xpos <- max(index)
      pos <- 4
      if (xvar == "penalty") {
        xpos <- min(index)
        pos <- 2
      }
      xpos <- rep(xpos, nnz)
      ypos <- gamma[, ncol(gamma)]
      text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
    }
  }
  par(op)
}
