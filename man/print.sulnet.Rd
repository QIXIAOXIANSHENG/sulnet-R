% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.sulnet.R
\name{print.sulnet}
\alias{print.sulnet}
\title{Print a sulnet object}
\usage{
\method{print}{sulnet}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{x}{fitted \code{\link{sulnet}} object}

\item{digits}{significant digits in printout}

\item{\dots}{additional print arguments}
}
\value{
a two-column matrix, the first columns is the number of nonzero
  coefficients and the second column is \code{Lambda}.
}
\description{
Print a summary of the sulnet path at each step along the path.
}
\details{
The call that produced the \code{\link{sulnet}} object is printed, followed
by a two-column matrix with columns \code{Df} and \code{Lambda}. The
\code{Df} column is the number of nonzero coefficients.
}
\examples{

data(FHT)
m1 <- sulnet(x = FHT$x, y = FHT$y, delta = 1, lambda2 = 0.1)
print(m1)

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
\author{
Yi Yang, Yuwen Gu and Hui Zou\cr

Maintainer: Yi Yang <yi.yang6@mcgill.ca>
}
\keyword{models}
\keyword{regression}
