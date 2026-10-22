#######################################################################
# This function is adapted/modified based on the plot.cv function
# from the `gcdnet` package:
#
#' Plot a cross-validation curve for a \code{cv.sulnet} object.
#'
#' Display mean validation loss, plus and minus one standard error, and the
#'   selected lambda values.
#'
#' @param x A \code{cv.sulnet} object.
#' @param sign.lambda Multiplier for log lambda on the horizontal axis; default
#'   1. Negative values reverse its direction.
#' @param ... Graphical arguments overriding the initial plot settings.
#'
#' @return Invisibly returns \code{NULL} after drawing the plot.
#' @importFrom graphics points abline
#' @method plot cv.sulnet
#' @export
plot.cv.sulnet <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0) {
    xlab <- paste("-", xlab, sep = "")
  }
  plot.args <- list(
    x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
    ylim = range(cvobj$cvupper, cvobj$cvlo), xlab = xlab,
    ylab = cvobj$name, type = "n"
  )
  new.args <- list(...)
  if (length(new.args)) {
    plot.args[names(new.args)] <- new.args
  }
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper,
    cvobj$cvlo,
    width = 0.01, col = "darkgrey"
  )
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm,
    pch = 20,
    col = "red"
  )
  axis(
    side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
    tick = FALSE, line = 0
  )
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}
