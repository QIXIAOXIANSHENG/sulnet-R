#' Fitting univariate response
#'
#' @param x Input X matrix of dimension nxp
#' @param y Input y vector of dimension n
#' @param loo Logical, whether apply leave one out or not. Default is \code{TRUE}
#'
#' @returns A list with \begin{\itemize}
#'  \item{\code{f}} Univariate fitted matrix with same dimension as \code{x}
#'  \item{\code{beta}} Univariate estimate of each variable
#'  \item{\code{b0}} Univariate intercept of each variable
#'  \item{\code{loo}} Logical, whether apply leave one out or not.
#'
#' @examples
#' data(FHT)
#' unifit <- uniFit(FHT$x, FHT$y_reg)
#' @export
uniFit <- function(x,y,loo = TRUE){
  y <- as.double(drop(y))
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)

  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])

  unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
                     beta0 = double(nvars),
                     beta  = double(nvars),
                     fit   = double(nobs * nvars),
                     PACKAGE = "sulnet")
  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
  return(list(fitted_values = f,
              beta = unifit$beta,
              b0 = unifit$beta0,
              loo = loo))
}
