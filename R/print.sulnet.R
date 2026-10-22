#' Print a fitted coefficient-path summary
#'
#' Display the stored call, number of nonzero slopes, and lambda for each path
#'   step.
#'
#' @param x A \code{sulnet} fit with vector components \code{df} and
#'   \code{lambda}.
#' @param digits Number of significant digits used for lambda.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns the printed matrix with columns Df and Lambda.
#' @method print sulnet
#' @export
print.sulnet <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}
