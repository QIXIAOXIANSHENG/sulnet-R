#' @import Matrix
#' @importFrom stats median

suniSRlogpath <- function(
  x, y, nsum, nrate, flmins, flminr, usum,
  urate, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3, maxit,
  lam2, nobs, nvars, vnames
) {
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
  if (!all(y %in% c(-1, 1))) {
    stop("y should be a factor with two levels")
  }
  y <- as.double(y)
  storage.mode(x) <- "double"

  medpj <- median(abs(beta[beta_ju != 0]))
  if (is.na(medpj)) medpj <- as.double(1)

  betaFIT <- .Fortran("betalogFIT", as.integer(nobs),
    as.integer(nvars), x, y,
    beta0 = double(nvars), beta = beta, beta_ju, medpj,
    fit = double(nobs * nvars), jerr = integer(nvars),
    PACKAGE = "sulnet"
  )
  betalogcheck(betaFIT$jerr, vnames)

  f <- matrix(betaFIT$fit, nobs, nvars)
  storage.mode(f) <- "double"

  lam <- getsrR(
    nobs, nvars, nsum, nrate, f, y, usum,
    urate, gridpf(pf, jd), flmins, flminr, isd, intr,
    family = "binomial"
  )
  usum <- as.double(lam$usum)
  urate <- as.double(sort(lam$urate))

  fit <- .Fortran("suniSRlogNET",
    lam2, nobs, nvars, f, y, jd, pf, pf2, pf3, pmax,
    nsum = nsum, nrate = nrate, usum = usum, urate = urate,
    eps, isd, intr, maxit,
    napenalty = integer(1),
    theta0 = double(1 * nsum * nrate), theta = double(pmax * nsum * nrate),
    itheta = integer(pmax * nsum * nrate), ntheta = integer(nsum * nrate),
    apenalty = double(nsum * nrate * 2), npass = integer(1), jerr = integer(1),
    PACKAGE = "sulnet"
  )
  output <- getoutputSR(fit, maxit, pmax, nvars, vnames)
  tol <- sqrt(.Machine$double.eps)
  which <- abs(output$penalty.path[, 4] - 0) < tol
  output$penalty.path[which, 4] <- 0
  gam <- theta2gamma(output, betaFIT)
  outlist <- list(
    gamma = gam$gamma,
    gamma0 = gam$gamma0,
    df = gam$df,
    dim = gam$dim,
    sum = gam$sum,
    rate = gam$rate,
    npass = fit$npass,
    jerr = fit$jerr,
    method = "sum_rate",
    inner = list(
      theta0 = output$theta0,
      theta = output$theta,
      f = f
    ),
    penalty.path = output$penalty.path,
    beta = betaFIT$beta,
    beta0 = betaFIT$beta0
  )
  class(outlist) <- c("sunilogpath")
  return(outlist)
}
