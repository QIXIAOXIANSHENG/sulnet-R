##' @import Matrix

suniZTpath <- function(x, y, nzeta, ntau, flminz, flmint, uzeta, utau, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3, maxit, lam2, nobs, nvars, vnames) {
  y <- as.double(y)
  storage.mode(x) <- "double"

  medpj <- median(abs(beta[beta_ju != 0]))
  if (is.na(medpj)) medpj <- as.double(1)

  betaFIT <- .Fortran("betaFIT", as.integer(nobs), as.integer(nvars), x, y, beta0 = double(nvars), beta = beta, beta_ju, medpj, fit = double(nobs * nvars), PACKAGE = "sulnet")

  f <- matrix(betaFIT$fit, nobs, nvars)
  storage.mode(f) <- "double"

  lam <- getztR(nobs, nvars, nzeta, ntau, f, y, uzeta, utau, pf, flminz, flmint, isd, intr, family = "gaussian")
  uzeta <- as.double(lam$uzeta)
  utau <- as.double(lam$utau)
  diffz <- as.double(lam$diffz)
  difft <- as.double(lam$difft)

  fit <- .Fortran("suniZTNET", lam2, nobs, nvars, f, y, jd, pf, pf2, pf3, pmax,
    nzeta = nzeta, ntau = ntau, uzeta = uzeta, utau = utau, eps, isd, intr, maxit,
    napenalty = integer(1),
    theta0 = double(1 * nzeta * ntau), theta = double(pmax * nzeta * ntau), itheta = integer(pmax * nzeta * ntau), ntheta = integer(nzeta * ntau), apenalty = double(nzeta * ntau * 2), npass = integer(1), diffz, difft, jerr = integer(1),
    PACKAGE = "sulnet"
  )
  output <- getoutputZT(fit, maxit, pmax, nvars, vnames)
  gam <- theta2gamma(output, betaFIT)
  outlist <- list(
    gamma = gam$gamma,
    gamma0 = gam$gamma0,
    df = gam$df,
    dim = gam$dim,
    zeta = gam$zeta,
    tau = gam$tau,
    npass = fit$npass,
    jerr = fit$jerr,
    use = "zeta_tau",
    inner = list(
      theta0 = output$theta0,
      theta = output$theta,
      f = f
    ),
    penalty.path = output$penalty.path,
    beta = betaFIT$beta,
    beta0 = betaFIT$beta0
  )
  class(outlist) <- c("sunipath")
  return(outlist)
}
