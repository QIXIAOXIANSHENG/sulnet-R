#' A testing code for conducting transfer learning
#'
#' @description
#' A short description...
#'
#'
#' @param d0 Target data, a list of two objects \code{x} and \code{y}
#' @param dA Source data, a list of all source data. Each source data should have the same
#'  structure as \code{d0}.
#' @param method Transferring method, default is \code{"suni_tsloo"}.
#' @param ... optional arguments to function \code{sulnet2D}.
#'
#' @returns A S3 class same as what \code{sulnet2D} with \code{negOnly = TRUE} would return.
#' @export
#'
#' @examples
#' n = 20
#' p = 10
#' beta = runif(p)
#' nA = 5
#' x = matrix(rnorm(n*p),nrow = n)
#' d0 = list(x = x,
#'           y = x %*% beta)
#' dA = lapply(seq(nA),function(a){
#'   x = matrix(rnorm(n*p),nrow = n)
#'   list(x = x,
#'   y = x %*% (beta + rnorm(p,0, 0.05)))
#' })
#' fit = transulnet(d0, dA)
transulnet = function(d0, dA, method = c("suni_tsloo","suni_sloo","suni_s"),
                      ...){
  ################################################################################
  ## target data setup
  method <- match.arg(method)
  this.call <- match.call()
  y0 <- drop(d0$y)
  x0 <- as.matrix(d0$x)
  np <- dim(x0)
  nobs0 <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x0)
  if (is.null(vnames))
    vnames <- paste("V", seq(nvars), sep = "")
  if (NROW(y0) != nobs0)
    stop("x and y have different number of observations")
  if (NCOL(y0) > 1L) stop("Multivariate response is not supported now")
  ################################################################################
  ## source data setup
  pA <- unique(sapply(dA, function(x){dim(x$x)[2]}))
  if(length(pA)!=1 | pA[1]!=nvars){
    stop("Source data has different number of variables from target data")
  }else{
    yA <- unlist(lapply(dA, function(item) item$y))
    xA <- do.call(rbind, lapply(dA, function(item) item$x))
  }
  ################################################################################
  ## Step 1, univariate fits
  nobsA <- nrow(xA)
  nobs <- as.integer(nobs0 + nobsA)
  new_y <- switch(method,
                  suni_tsloo = as.double(c(y0,yA)),
                  suni_sloo = as.double(y0),
                  suni_s = as.double(y0))
  unifit <- switch(method,
                   suni_tsloo = .Fortran("loofit_ts", nobs, nobs0, nvars, as.double(x0), as.double(y0),
                                         as.double(rbind(x0,xA)),
                                         as.double(c(y0,yA)), loo,
                                         beta0 = double(nvars),
                                         beta = double(nvars),
                                         fit = double(nobs * nvars),
                                         PACKAGE = "sulnet"),
                   suni_sloo  = .Fortran("loofit_s", nobs, nobs0, nvars, as.double(x0), as.double(y0),
                                         as.double(rbind(x0,xA)),
                                         as.double(c(y0,yA)), loo,
                                         beta0 = double(nvars),
                                         beta = double(nvars),
                                         fit = double(nobs0 * nvars),
                                         PACKAGE = "sulnet"),
                   suni_s  = .Fortran("loofit_st", nobs0, nvars, nobsA, as.double(x0), as.double(y0),
                                      as.double(xA),
                                      as.double(yA),
                                      beta0 = double(nvars),
                                      beta = double(nvars),
                                      fit = double(nobs0 * nvars),
                                      PACKAGE = "sulnet")
  )
  unibeta = unifit$beta
  unibeta0 = unifit$beta0
  f <- matrix(unifit$fit, ncol = nvars)
  colnames(f) = vnames


  ################################################################################
  ## parameter setup
  args <- list(...)
  args$negOnly <- TRUE
  args$x <- f
  args$y <- new_y

  ################################################################################
  ## Step 2, lasso using univariate fits
  lsfit = do.call(sulnet2D, args)

  ################################################################################
  ## Step 3, lasso using univariate fits
  forseq = ifelse(lsfit$use_alpha, length(lsfit$alpha), length(lsfit$lamPos))
  ones = rep(1,length(lsfit$lambda))
  unibeta <- outer(unibeta, ones)
  unibeta0 <- outer(unibeta0, ones)
  new_lsfit = lapply(seq(forseq), function(x){
    multuni(lsfit$beta[[x]], lsfit$b0[[x]],
            unibeta, unibeta0)
  })
  beta_list <- lapply(new_lsfit, function(x) x$beta)
  names(beta_list) = names(lsfit$beta)
  b0_list <- lapply(new_lsfit, function(x) x$b0)
  names(b0_list)  = names(lsfit$b0)

  outlist = lsfit
  outlist$beta = beta_list
  outlist$b0 = b0_list
  outlist$fbeta = lsfit$beta
  outlist$fb0 = lsfit$b0
  class(outlist) = "sulnet2D.sunipath_2"
  outlist
}

