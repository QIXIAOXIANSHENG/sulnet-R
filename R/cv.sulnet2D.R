##' Cross-validation for sulnet2D
##'
##' Does k-fold cross-validation for sulnet, produces a plot, and returns a
##' value for \code{lambda}. This function is modified based on the \code{cv}
##' function from the \code{glmnet} package.
##'
##' The function runs \code{\link{sulnet2D}} \code{nfolds}+1 times; the first to
##' get the \code{lambda} sequence, and then the remainder to compute the fit
##' with each of the folds omitted. The average error and standard deviation
##' over the folds are computed.
##'
##' @aliases cv.sulnet cv.hsvmpath cv.sqsvmpath cv.logitpath cv.lspath cv.erpath
##'
##' @param x \code{x} matrix as in \code{\link{sulnet2D}}.
##' @param y response variable or class label \code{y} as in
##'   \code{\link{sulnet}}.
##' @param lambda optional user-supplied lambda sequence; default is
##'   \code{NULL}, and \code{\link{sulnet}} chooses its own sequence.
##' @param alpha alpha.
##' @param nfolds number of folds - default is 5. Although \code{nfolds} can be
##'   as large as the sample size (leave-one-out CV), it is not recommended for
##'   large datasets. Smallest value allowable is \code{nfolds=3}.
##' @param foldid an optional vector of values between 1 and \code{nfold}
##'   identifying what fold each observation is in. If supplied, \code{nfold}
##'   can be missing.
##' @param pred.loss loss function to use for cross-validation error. Valid
##'   options are: \itemize{ \item \code{"loss"} Margin based loss function.
##'   When use least square loss \code{"ls"}, it gives mean square error (MSE).
##'   When use expectile regression loss \code{"er"}, it gives asymmetric mean
##'   square error (AMSE). \item \code{"misclass"} only available for
##'   classification: it gives misclassification error. } Default is
##'   \code{"loss"}.
##' @param \dots other arguments that can be passed to sulnet.
##' @return an object of class \code{\link{cv.sulnet}} is returned, which is a
##'   list with the ingredients of the cross-validation fit. \item{lambda}{the
##'   values of \code{lambda} used in the fits.} \item{cvm}{the mean
##'   cross-validated error - a vector of length \code{length(lambda)}.}
##'   \item{cvsd}{estimate of standard error of \code{cvm}.}
##'   \item{cvupper}{upper curve = \code{cvm+cvsd}.} \item{cvlower}{lower curve
##'   = \code{cvm-cvsd}.} \item{nzero}{number of non-zero coefficients at each
##'   \code{lambda}.} \item{name}{a text string indicating type of measure (for
##'   plotting purposes).} \item{sulnet.fit}{a fitted \code{\link{sulnet}}
##'   object for the full data.} \item{lambda.min}{The optimal value of
##'   \code{lambda} that gives minimum cross validation error \code{cvm}.}
##'   \item{lambda.1se}{The largest value of \code{lambda} such that error is
##'   within 1 standard error of the minimum.}
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr Maintainer: Yi Yang
##'   <yi.yang6@mcgill.ca>
##' @seealso \code{\link{sulnet}}, \code{\link{plot.cv.sulnet}},
##'   \code{\link{predict.cv.sulnet}}, and \code{\link{coef.cv.sulnet}} methods.
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
##' # fit an elastic net penalized HHSVM with lambda2 = 0.1 for the L2 penalty.
##' # Use the misclassification rate as the cross validation prediction loss.
##' # Use five-fold CV to choose the optimal lambda for the L1 penalty.
##'
##' data(FHT)
##'
##' # fit an elastic net penalized least squares
##' # with lambda2 = 0.1 for the L2 penalty. Use the
##' # least square loss as the cross validation
##' # prediction loss. Use five-fold CV to choose
##' # the optimal lambda for the L1 penalty.
##'
##' set.seed(2011)
##' cv1 <- cv.sulnet2D(FHT$x, FHT$y_reg, method ="ls",
##'                  lambda2 = 0.1, pred.loss = "loss",
##'                  nfolds = 5)
##' plot(cv1)
##'

##'
##' @export
cv.sulnet2D <- function(x, y, lambda = NULL, alpha = seq(0,0.5,length.out = 11),
                      pred.loss = c("misclass", "loss"),
                      nfolds = 5, foldid, ...) {
  if (missing(pred.loss))
    pred.loss <- "default" else pred.loss <- match.arg(pred.loss)

  N <- nrow(x)
  ## Fit the model once to get dimensions etc of output
  y <- drop(y)
  sulnet2D.object <- sulnet2D(x, y, lambda = lambda,alpha = alpha,
                          ...)
  lambda <- sulnet2D.object$lambda

  if(sulnet2D.object$use_alpha){
    alpha = sulnet2D.object$alpha
  }else{
    alpha = sulnet2D.object$lamPos
  }

  ## predict -> coef
  nz <- sapply(coef(sulnet2D.object, type = "nonzero"),
               function(i){lapply(i, function(x){sapply(x,length)})})
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    if(sulnet2D.object$use_alpha){
      outlist[[i]] <- sulnet2D(x = x[!which, , drop = FALSE],
                               y = y_sub, lambda = lambda,alpha = alpha,
                               ...)
    }else{
      outlist[[i]] <- sulnet2D(x = x[!which, , drop = FALSE],
                               y = y_sub, lambda = lambda,alpha = NULL,
                               ...)
    }

  }
  ## What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(sulnet2D.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, alpha, x, y, foldid,
                               pred.loss))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  if(sulnet2D.object$use_alpha){
    out <- list(lambda = lambda, alpha = alpha, cvm = cvm, cvsd = cvsd,
                cvupper = cvm + cvsd, cvlower = cvm - cvsd,
                nzero = nz,
                name = cvname, sulnet2D.fit = sulnet2D.object)
  }else{
    out <- list(lambda = lambda, lamPos = alpha, cvm = cvm, cvsd = cvsd,
                cvupper = cvm + cvsd, cvlower = cvm - cvsd,
                nzero = nz,
                name = cvname, sulnet2D.fit = sulnet2D.object)
  }
  lamin <- getmin2D(lambda, alpha, cvm, cvsd, sulnet2D.object$use_alpha )
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sulnet2D"
  obj
}
