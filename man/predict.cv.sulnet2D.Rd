% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.cv.sulnet2D.R
\name{predict.cv.sulnet2D}
\alias{predict.cv.sulnet2D}
\title{Make predictions from a "cv.sulnet" object.}
\usage{
\method{predict}{cv.sulnet2D}(object, newx, s = c("cv.1se", "cv.min"), ...)
}
\arguments{
\item{object}{fitted \code{\link{cv.sulnet}} object.}

\item{newx}{matrix of new values for \code{x} at which predictions are to be
made. Must be a matrix. See documentation for \code{predict.sulnet}.}

\item{s}{value(s) of the penalty parameter \code{lambda} at which
predictions are required. Default is the value \code{s="lambda.1se"}
stored on the CV object. Alternatively \code{s="lambda.min"} can be used.
If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
used.}

\item{\dots}{not used. Other arguments to predict.}
}
\value{
The object returned depends the \dots{} argument which is passed on
  to the \code{\link{predict}} method for \code{\link{sulnet}} objects.
}
\description{
This function makes predictions from a cross-validated sulnet model, using
the stored \code{"sulnet2D.fit"} object, and the optimal value chosen for
\code{lambda} and \code{alpha}/\code{lamPos}.
}
\details{
This function makes it easier to use the results of cross-validation to
make a prediction.
}
\examples{

data(FHT)
set.seed(2011)
cv=cv.sulnet2D(FHT$x, FHT$y, lambda2 = 1, pred.loss="misclass",
             lambda.factor=0.05, nfolds=5)
pre = predict(cv$sulnet.fit, newx = FHT$x, s = "cv.min")

}
\references{
Yang, Y. and Zou, H. (2012).
  "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
  \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
  BugReport: \url{https://github.com/emeryyi/sulnet}\cr

  Gu, Y., and Zou, H. (2016).
  "High-dimensional generalizations of asymmetric least squares regression and their applications."
  \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr

  Friedman, J., Hastie, T., and Tibshirani, R. (2010).
  "Regularization paths for generalized linear models via coordinate descent."
  \emph{Journal of Statistical Software, 33, 1.}\cr
  \url{https://www.jstatsoft.org/v33/i01/}
}
\seealso{
\code{\link{cv.sulnet}}, and \code{\link{coef.cv.sulnet}} methods.
}
\author{
Yi Yang, Yuwen Gu and Hui Zou\cr

Maintainer: Yi Yang <yi.yang6@mcgill.ca>
}
\keyword{models}
\keyword{regression}
