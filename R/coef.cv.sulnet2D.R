##' Get coefficients or make coefficient predictions from a "cv.sulnet" object.
##'
##' This function gets coefficients or makes coefficient predictions from a
##' cross-validated sulnet model, using the stored \code{"sulnet.fit"} object,
##' and the optimal value chosen for \code{lambda}.
##'
##' This function makes it easier to use the results of cross-validation to get
##' coefficients or make coefficient predictions.
##'
##' @param object fitted \code{\link{cv.sulnet}} object.
##' @param s value(s) of the penalty parameter \code{lambda} at which
##'   predictions are required. Default is the value \code{s="lambda.1se"}
##'   stored on the CV \code{object}, it is the largest value of \code{lambda}
##'   such that error is within 1 standard error of the minimum. Alternatively
##'   \code{s="lambda.min"} can be used, it is the optimal value of
##'   \code{lambda} that gives minimum cross validation error \code{cvm}. If
##'   \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
##'   used.
##' @param \dots not used. Other arguments to predict.
##'
##' @return The object returned depends the \dots{} argument which is passed on
##'   to the \code{\link{predict}} method for \code{\link{sulnet}} objects.
##'
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##'   Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##'
##' @seealso \code{\link{cv.sulnet}}, and \code{\link{predict.cv.sulnet}}
##'   methods.
##'
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
##'
##' @keywords models regression
##'
##' @examples
##'
##' data(FHT)
##' set.seed(2011)
##' cv <- cv.sulnet2D(FHT$x, FHT$y, nfolds = 5)
##' coef(cv, s = "lambda.min")
##'
##' @export
coef.cv.sulnet2D <- function(object, s = c("cv.1se", "cv.min"), ...) {
  if (is.character(s)) {
      s <- match.arg(s)
      la <- object[[s]]
    } else stop("Invalid form for s")
  coef(object$sulnet2D.fit, s = la[[1]], alpha = la[[2]], ...)
}
