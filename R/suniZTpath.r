##' @import Matrix

suniZTpath <- function(x, y, nzeta, ntau, flminz, flmint, uzeta, utau, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3, maxit, lam2, loo, nobs, nvars, vnames) {
  y <- as.double(y)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)

  betaFIT <- .Fortran("betaFIT", as.integer(nobs), as.integer(nvars), x, y, beta0 = double(nvars), beta = beta, beta_ju, as.double(median(beta[beta_ju == 1])), fit = double(nobs * nvars), PACKAGE = "sulnet")

  lam <- getztR(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family = "gaussian")
  uzeta <- lam$uzeta
  utau <- lam$utau
}
