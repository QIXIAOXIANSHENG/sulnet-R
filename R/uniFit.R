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
uniFit <- function(x,y,loo = TRUE,family = c("gaussian", "binomial")){
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)
  family <- match.arg(family)

  np = dim(x)
  nobs = as.integer(np[1])
  nvars = as.integer(np[2])
  if(nobs != length(y)){
    stop("X and y have different number of rows in call to uniFit",call.=FALSE)
  }
  if(family == "binomial"){
    y <- factor(y)
    if(length(levels(y)) != 2){
      stop("There should only be two classes in y")
    }
    levels(y) <- c(0,1)
    y <- as.double(drop(as.vector(y)))
  }
  unifit <- switch(family,
                  gaussian = .Fortran("loofit", nobs, nvars, x, y, loo,
                     beta0 = double(nvars),
                     beta  = double(nvars),
                     fit   = double(nobs * nvars),
                     PACKAGE = "sulnet"),
                  binomial = .Fortran("loofit_binom", nobs, nvars, x, y, loo,
                     beta0 = double(nvars),
                     beta  = double(nvars),
                     fit   = double(nobs * nvars),
                     PACKAGE = "sulnet")
                  )

  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
  return(list(fitted_values = f,
              beta = unifit$beta,
              b0 = unifit$beta0,
              loo = loo))
}
