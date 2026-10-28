#' @import Matrix
#' @importFrom stats median

suniLAlogpath <- function(
  x, y, nlambda, nalpha, flmin, flmina,
  ulambda, ualpha, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3,
  maxit, lam2, nobs, nvars, vnames
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

  betaFIT <- .Fortran("betalogFIT",
    as.integer(nobs), as.integer(nvars), x, y,
    beta0 = double(nvars), beta = beta, beta_ju, medpj,
    fit = double(nobs * nvars),
    PACKAGE = "sulnet"
  )

  f <- matrix(betaFIT$fit, nobs, nvars)
  storage.mode(f) <- "double"

  lam <- getlaR(
    nobs, nvars, nlambda, nalpha, f, y, ulambda, ualpha,
    pf, flmin, flmina, intr,
    family = "binomial"
  )
  ulambda <- as.double(lam$ulambda)
  ualpha <- as.double(sort(lam$ualpha))

  fit <- .Fortran("suniLAlogNET",
    lam2, nobs, nvars, f, y, jd, pf, pf2, pf3, pmax,
    nlambda = nlambda, nalpha = nalpha, ulambda = ulambda, ualpha = ualpha,
    eps, isd, intr, maxit,
    napenalty = integer(1),
    theta0 = double(1 * nlambda * nalpha),
    theta = double(pmax * nlambda * nalpha),
    itheta = integer(pmax * nlambda * nalpha),
    ntheta = integer(nlambda * nalpha),
    apenalty = double(nlambda * nalpha * 2), npass = integer(1),
    jerr = integer(1),
    PACKAGE = "sulnet"
  )
  output <- getoutputLA(fit, maxit, pmax, nvars, vnames)
  tol <- sqrt(.Machine$double.eps)
  which <- abs(output$penalty.path[, 4] - 0) < tol
  output$penalty.path[which, 4] <- 0
  gam <- theta2gamma(output, betaFIT)
  outlist <- list(
    gamma = gam$gamma,
    gamma0 = gam$gamma0,
    df = gam$df,
    dim = gam$dim,
    lambda = gam$lambda,
    alpha = gam$alpha,
    npass = fit$npass,
    jerr = fit$jerr,
    method = "lambda_alpha",
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
