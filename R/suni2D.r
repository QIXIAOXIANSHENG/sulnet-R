suni2D <- function(
  x, y, beta = NULL,
  method = c("zeta_tau", "lambda_alpha", "sum_rate"),
  family = c("gaussian", "binomial"),
  nzeta = 100, ntau = 20, nlambda = 100, nalpha = 20, nsum = 100, nrate = 20,
  lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  zeta = NULL, tau = NULL, lambda = NULL, alpha = NULL, sum = NULL, rate = NULL,
  lambda2 = 0, pf = rep(1, nvars),
  pf2 = rep(1, nvars), pf3 = rep(1, nvars), exclude,
  pmax = nvars, standardize = TRUE,
  intercept = TRUE, eps = 1e-08, maxit = 1e+05
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
  if (length(pf3) != nvars) {
    stop("The size of L1 penalty factor must be same as the number of input variables")
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
  if (method == "zeta_tau") {
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
  if (method == "sum_rate") {
    nsum <- as.integer(nsum)
    nrate <- as.integer(nrate)
    if (is.null(sum)) {
      lambda.factor <- as.double(lambda.factor)
      if (lambda.factor >= 1) {
        stop("lambda.factor should be less than 1")
      }
      if (lambda.factor <= 0) {
        stop("lambda.factor should be greater than 0")
      }
      if (nsum < 3) {
        message("nzeta should be at least 3, set to 3")
        nsum <- 3
      }
      flmins <- as.double(lambda.factor)
      usum <- double(nsum)
    } else {
      flmins <- as.double(1)
      if (any(sum < 0)) {
        stop("zetas should be non-negative")
      }
      usum <- as.double(rev(sort(sum)))
      nsum <- as.integer(length(sum))
    }
    if (is.null(rate)) {
      lambda.factor <- as.double(lambda.factor)
      if (lambda.factor >= 1) {
        stop("lambda.factor should be less than 1")
      }
      if (lambda.factor <= 0) {
        stop("lambda.factor should be greater than 0")
      }
      flminr <- as.double(lambda.factor)
      urate <- double(nrate)
    } else {
      flminr <- as.double(1)
      if (any(rate < 1)) {
        stop("rate should be larger than 1")
      } else if (any(is.infinite(rate))) {
        if (nobs < nvars) {
          stop("rate should not be infinite in high dimension settings")
        } else {
          warning("using alpha = 1 for infinite rate cases")
          method <- "lambda_alpha"
          ualpha <- as.double(rate / (1 + rate))
          ualpha[which(is.na(ualpha))] <- as.double(1)
          ualpha <- as.double(sort(unique(ualpha)))
          alpha <- ualpha
          lambda <- zeta * 2
        }
      }
      urate <- as.double(sort(unique(rate)))
      nrate <- as.integer(length(urate))
    }
  }
  if (method == "lambda_alpha") {
    nlambda <- as.integer(nlambda)
    nalpha <- as.integer(nalpha)
    if (is.null(lambda)) {
      lambda.factor <- as.double(lambda.factor)
      if (lambda.factor >= 1) {
        stop("lambda.factor should be less than 1")
      }
      if (lambda.factor <= 0) {
        stop("lambda.factor should be greater than 0")
      }
      if (nlambda < 3) {
        message("nzeta should be at least 3, set to 3")
        nlambda <- 3
      }
      flmin <- as.double(lambda.factor)
      ulambda <- double(nlambda)
    } else {
      flmin <- as.double(1)
      if (any(lambda < 0)) {
        stop("zetas should be non-negative")
      }
      ulambda <- as.double(rev(sort(lambda)))
      nlambda <- as.integer(length(lambda))
    }
    if (is.null(alpha)) {
      lambda.factor <- as.double(lambda.factor)
      if (lambda.factor >= 1) {
        stop("lambda.factor should be less than 1")
      }
      if (lambda.factor <= 0) {
        stop("lambda.factor should be greater than 0")
      }
      flmina <- as.double(lambda.factor)
      ualpha <- double(nalpha)
    } else {
      flmina <- as.double(1)
      if (any(alpha < 0.5)) {
        stop("alpha should be larger than 0.5")
      } else if (any(alpha > 1)) {
        stop("alpha should be smaller than 1")
      } else if (any(alpha == 1) && nobs < nvars) {
        stop("alpha should not be equal to 1 in high dimension settings")
      }
      ualpha <- as.double(sort(unique(alpha)))
      nalpha <- as.integer(length(ualpha))
    }
  }

  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) {
      stop("Some excluded variables out of range")
    }
    jd <- as.integer(c(length(jd), jd))
  } else {
    jd <- as.integer(0)
  }

  fit <- switch(method,
    zeta_tau = suniZTpath(
      x, y, nzeta, ntau, flminz, flmint, uzeta,
      utau, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3,
      maxit, lam2, nobs, nvars, vnames
    ),
    lambda_alpha = suniLApath(
      x, y, nlambda, nalpha, flmin, flmina,
      ulambda, ualpha, isd, intr, beta, beta_ju, eps, pmax, jd, pf,
      pf2, pf3, maxit, lam2, nobs, nvars, vnames
    ),
    sum_rate = suniSRpath(
      x, y, nsum, nrate, flmins, flminr, usum,
      urate, isd, intr, beta, beta_ju, eps, pmax, jd, pf, pf2, pf3,
      maxit, lam2, nobs, nvars, vnames
    )
  )

  fit$call <- this.call
  ################################################################################
  class(fit) <- c("suni2D", class(fit))
  fit
}
