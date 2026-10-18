##' @export
predict.suni2D <- function(
  object, newx, penalty1 = NULL,
  penalty2 = NULL, type = c("class", "link"),
  ...
) {
  NextMethod("predict")
}

##' Model predictions
##'
##' \code{predict} is a generic function for predictions from the results of
##' various model fitting functions. The function invokes particular
##' \emph{methods} which depend on the \code{\link{class}} of the first
##' argument.
##'
##' @param object a model object for which prediction is desired.
##' @param \dots additional arguments affecting the predictions produced.
##'
##' @return The form of the value returned by \code{predict} depends on the
##'   class of its argument. See the documentation of the particular methods for
##'   details of what is produced by that method.
##'
##' @seealso \code{\link{predict.sulnet}}, \code{\link{predict.erpath}},
##'   \code{\link{predict.lspath}}, \code{\link{predict.hsvmpath}},
##'   \code{\link{predict.logitpath}}, \code{\link{predict.sqsvmpath}}.
##'
##' @export
predict <- function(object, ...) UseMethod("predict")
