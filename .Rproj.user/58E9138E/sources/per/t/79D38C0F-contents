##' @import Matrix

sunipath <- function(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, pf,
                   pf2, maxit, lam2, lamPos, loo, negOnly, nobs, nvars, vnames,
                   alpha, ignore_lamPos) {
  ################################################################################
  ## data setup
  y <- as.double(y)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)

  ################################################################################
  ## if only computing the negative steps

  if(negOnly){
    if(ignore_lamPos){
      fit <- .Fortran("suniwalpha", lam2, lamPos, nobs, nvars, x, as.double(y), jd, pf, pf2,
                      dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                      nalam = integer(1), b0 = double(nlam),
                      beta = double(pmax * nlam), ibeta = integer(pmax),
                      nbeta = integer(nlam), alam = double(nlam),
                      npass = integer(1), jerr = integer(1),
                      alpha = as.double(alpha), iglamPos = as.logical(ignore_lamPos),
                      PACKAGE = "sulnet")
    }else{
      fit <- .Fortran("suniwalpha", lam2,lamPos, nobs, nvars, x, as.double(y), jd, pf, pf2,
                      dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                      nalam = integer(1), b0 = double(nlam),
                      beta = double(pmax * nlam), ibeta = integer(pmax),
                      nbeta = integer(nlam), alam = double(nlam),
                      npass = integer(1), jerr = integer(1),
                      alpha = as.double(0), iglamPos = as.logical(ignore_lamPos),
                      PACKAGE = "sulnet")
    }
    outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
    if(ignore_lamPos){
      outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr,
                                 alpha = alpha, negOnly = negOnly))
    }else{
      outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr,
                                 lamPos = lamPos, negOnly = negOnly))
    }

    class(outlist) <- c("sunipath")
    return(outlist)

  }

  ################################################################################
  ## univariate fit

  unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
                     beta0 = double(nvars),
                     beta  = double(nvars),
                     fit   = double(nobs * nvars),
                     PACKAGE = "sulnet")
  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)


  ################################################################################
  ## call Fortran core
  if(ignore_lamPos){

    fit <- .Fortran("suniwalpha", lam2, lamPos, nobs, nvars, f, as.double(y), jd, pf, pf2,
                    dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                    nalam = integer(1), b0 = double(nlam),
                    beta = double(pmax * nlam), ibeta = integer(pmax),
                    nbeta = integer(nlam), alam = double(nlam),
                    npass = integer(1), jerr = integer(1),
                    alpha = as.double(alpha), iglamPos = as.logical(ignore_lamPos),
                    PACKAGE = "sulnet")
  }else{
    fit <- .Fortran("suniwalpha", lam2,lamPos, nobs, nvars, f, as.double(y), jd, pf, pf2,
                    dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                    nalam = integer(1), b0 = double(nlam),
                    beta = double(pmax * nlam), ibeta = integer(pmax),
                    nbeta = integer(nlam), alam = double(nlam),
                    npass = integer(1), jerr = integer(1),
                    alpha = as.double(0), iglamPos = as.logical(ignore_lamPos),
                    PACKAGE = "sulnet")
  }
  ################################################################################
  ## output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))

  if(ignore_lamPos){
    outlist[["alpha"]] = alpha
  }else{
    outlist[["lamPos"]] = lamPos
  }
  outlist <- c(outlist, list(LOO = loo,
                             univariate.fit = list(beta = unifit$beta,
                                                   beta0 = unifit$beta0,
                                                   fitted.values = f)))

  ones = rep(1,fit$nalam)
  unibeta <- outer(unifit$beta, ones)
  unibeta0 <- outer(unifit$beta0, ones)
  beta_temp = outlist$beta
  beta0_temp = outlist$b0

  row_idx <- beta_temp@i + 1
  col_ptrs <- beta_temp@p
  col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
  beta_result <- beta_temp
  beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]

  outlist$beta = beta_result
  outlist$b0 = beta0_temp + colSums(unibeta0 * beta_temp)

  outlist$fbeta = beta_temp; outlist$fb0 = beta0_temp
  class(outlist) <- c("sunipath")
  outlist
}


