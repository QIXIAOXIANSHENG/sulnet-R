##' Make predictions from a "cv.sulnet" object.
##'
##' This function makes predictions from a cross-validated sulnet model, using
##' the stored \code{"sulnet2D.fit"} object, and the optimal value chosen for
##' \code{lambda} and \code{alpha}/\code{lamPos}.
##'
##' This function makes it easier to use the results of cross-validation to
##' make a prediction.
##'
##' @param object fitted \code{\link{cv.sulnet}} object.
##' @param newx matrix of new values for \code{x} at which predictions are to be
##'   made. Must be a matrix. See documentation for \code{predict.sulnet}.
##' @param s value(s) of the penalty parameter \code{lambda} at which
##'   predictions are required. Default is the value \code{s="lambda.1se"}
##'   stored on the CV object. Alternatively \code{s="lambda.min"} can be used.
##'   If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
##'   used.
##' @param \dots not used. Other arguments to predict.
##' @return The object returned depends the \dots{} argument which is passed on
##'   to the \code{\link{predict}} method for \code{\link{sulnet}} objects.
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##' @seealso \code{\link{cv.sulnet}}, and \code{\link{coef.cv.sulnet}} methods.
##' @references Yang, Y. and Zou, H. (2012).
##'   "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
##'   \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
##'   BugReport: \url{https://github.com/emeryyi/sulnet}\cr
##'
##'   Gu, Y., and Zou, H. (2016).
##'   "High-dimensional generalizations of asymmetric least squares regression and their applications."
##'   \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr
##'
##'   Friedman, J., Hastie, T., and Tibshirani, R. (2010).
##'   "Regularization paths for generalized linear models via coordinate descent."
##'   \emph{Journal of Statistical Software, 33, 1.}\cr
##'   \url{https://www.jstatsoft.org/v33/i01/}
##' @keywords models regression
##' @examples
##'
##' data(FHT)
##' set.seed(2011)
##' cv=cv.sulnet2D(FHT$x, FHT$y, lambda2 = 1, pred.loss="misclass",
##'              lambda.factor=0.05, nfolds=5)
##' pre = predict(cv$sulnet.fit, newx = FHT$x, s = "cv.min")
##'
##' @export
predict.cv.sulnet2D <- function(object, newx, s = c("cv.1se", "cv.min"),
                              ...) {
  if (is.character(s)) {
      s <- match.arg(s)
      if(s == "cv.1se"){
        lambda <- object$cv.1se$lambda.1se
        alpha <- object$cv.1se[[2]]
      }else{
        lambda <- object$cv.min$lambda.min
        alpha <- object$cv.min[[2]]
      }

    } else stop("Invalid form for s")
  predict(object$sulnet2D.fit, newx, s = lambda, alpha = alpha, ...)
}
