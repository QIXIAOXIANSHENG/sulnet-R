suni2D <- function(
  x, y, beta = NULL,
  method = c("zeta-tau", "lambda-alpha", "sum-rate"), family = c("gaussian", "binomial"),
  nzeta = 100, ntau = 20, nlambda = 100, nalpha = 20, nsum = 100, nrate = 20,
  lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  zeta = NULL, tau = NULL, lambda = NULL, alpha = NULL, sum = NULL, rate = NULL,
  lambda2 = 0, pf = rep(1, nvars),
  pf2 = rep(1, nvars), exclude, dfmax = nvars + 1,
  pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
  intercept = TRUE, eps = 1e-08, maxit = 1e+05, loo = FALSE
) {
  method <- match.arg(method)
  family <- match.arg(family)
  this.call <- match.call()
  y <- drop(y)
  x <- as.matrix(x)
  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  vnames <- colnames(x)
  if (is.null(vnames)) {
    vnames <- paste("V", seq(nvars), sep = "")
  }
  if (NROW(y) != nobs) {
    stop("x and y have different number of observations")
  }
  if (NCOL(y) > 1L && family != "cox") stop("Multivariate response is not supported now")
  ## parameter setup
  if (length(pf) != nvars) {
    stop("The size of L1 penalty factor must be same as the number of input variables")
  }
  if (length(pf2) != nvars) {
    stop("The size of L2 penalty factor must be same as the number of input variables")
  }
  if (lambda2 < 0) {
    stop("lambda2 must be non-negative")
  }
  lam2 <- as.double(lambda2)
  pf <- as.double(pf)
  pf2 <- as.double(pf2)
  pf3 <- as.double(pf3)
  isd <- as.integer(standardize)
  intr <- as.integer(intercept)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  maxit <- as.integer(maxit)
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) {
      stop("Some excluded variables out of range")
    }
    jd <- as.integer(c(length(jd), jd))
  } else {
    jd <- as.integer(0)
  }
  if (!is.null(beta)) {
    if (is.list(beta)) {
      index <- as.integer(pmin(pmax(1, beta[[1]]), nvars))
      if (anyDuplicated(index) > 0) {
        stop("Duplicate indices in beta[[1]] are not allowed")
      }
      if (length(beta[[2]]) != length(index)) {
        stop("The size of beta[[2]] must be same as the number of indices in beta[[1]]")
      }
      beta_full <- double(nvars)
      beta_full[index] <- as.double(beta[[2]])
      beta <- beta_full
      beta_ju <- integer(nvars)
      beta_ju[index] <- as.integer(1)
    } else if (is.vector(beta)) {
      if (length(beta) != nvars) {
        stop("The size of beta must be same as the number of input variables")
      }
      beta_ju <- as.integer(which(!is.na(beta)))
      beta[-beta_ju] <- 0
      beta <- as.double(beta)
    } else {
      stop("beta must be either a list or a vector")
    }
    beta_ju[beta == 0] <- as.integer(0)
  } else {
    beta_ju <- integer(nvars)
    beta <- double(nvars)
  }


  ## lambda setup
  nzeta <- as.integer(nzeta)
  ntau <- as.integer(ntau)
  if (is.null(zeta)) {
    lambda.factor <- as.double(lambda.factor)
    if (lambda.factor >= 1) {
      stop("lambda.factor should be less than 1")
    }
    if (lambda.factor <= 0) {
      stop("lambda.factor should be greater than 0")
    }
    if (nzeta < 3) {
      message("nzeta should be at least 3, set to 3")
      nzeta <- 3
    }
    flminz <- as.double(lambda.factor)
    uzeta <- double(nzeta)
  } else {
    flminz <- as.double(1)
    if (any(zeta < 0)) {
      stop("zetas should be non-negative")
    }
    uzeta <- as.double(rev(sort(zeta)))
    nzeta <- as.integer(length(zeta))
  }
  if (is.null(tau)) {
    lambda.factor <- as.double(lambda.factor)
    if (lambda.factor >= 1) {
      stop("lambda.factor should be less than 1")
    }
    if (lambda.factor <= 0) {
      stop("lambda.factor should be greater than 0")
    }
    if (ntau < 3) {
      message("ntau should be at least 3, set to 3")
      ntau <- 3
    }
    flmint <- as.double(lambda.factor)
    utau <- double(ntau)
  } else {
    flmint <- as.double(1)
    if (any(tau < 0)) {
      stop("taus should be non-negative")
    }
    utau <- as.double(rev(sort(tau)))
    ntau <- as.integer(length(tau))
  }
}
