library(devtools)
clean_dll()
document()
load_all()
set.seed(1)
nobs <- as.integer(100)
nvars <- as.integer(10)
nzeta <- as.integer(10)
ntau <- as.integer(5)
nlambda <- as.integer(10)
nalpha <- as.integer(5)
nsum <- as.integer(10)
nrate <- as.integer(5)
x <- matrix(rnorm(nobs * nvars), nrow = nobs, ncol = nvars)
storage.mode(x) <- "double"
y <- as.double(rnorm(nobs)) + x %*% rep(c(2, 0), each = nvars / 2)
pf <- as.double(rep(1, nvars))
pf2 <- as.double(rep(1, nvars))
pf3 <- as.double(rep(1, nvars))
beta <- as.double(c(runif(5, -1, 1), rep(0, nvars - 5)))
beta_ju <- as.integer(rep(c(1, 0), c(5, nvars - 5)))

########################
beta = NULL
uzeta <- double(nzeta)
utau <- double(ntau)
ulambda <- double(nlambda)
ualpha <- double(nalpha)
usum <- double(nsum)
urate <- double(nrate)
flminz <- 0.0001
flmint <- 0.0001
flmin <- 0.0001
flmina <- 0.0001
flmins <- 0.0001
flminr <- 0.0001
family <- "gaussian"
method = "zeta_tau"
y <- drop(y)
x <- as.matrix(x)
np <- dim(x)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(x)
if (is.null(vnames)) {
  vnames <- paste("V", seq(nvars), sep = "")
}
lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04)
lam2 <- as.double(0)
isd <- as.integer(1)
intr <- as.integer(1)
eps <- as.double(1e-6)
pmax <- as.integer(nvars)
maxit <- as.integer(1e+05 * ntau)
jd <- as.integer(0)
y <- as.double(y)
storage.mode(x) <- "double"
betaFIT <- .Fortran("betaFIT", as.integer(nobs), as.integer(nvars), x, y, beta0 = double(nvars), beta = beta, beta_ju, as.double(median(abs(beta[beta_ju == 1]))), fit = double(nobs * nvars), PACKAGE = "sulnet")

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



lam <- getlaR(nobs, nvars, nlambda, nalpha, f, y, ulambda, ualpha, pf, flmin, flmina, isd, intr)
ulambda <- as.double(lam$ulambda)
ualpha <- as.double(sort(lam$ualpha))
fit <- .Fortran("suniLANET", lam2, nobs, nvars, f, y, jd, pf, pf2, pf3, pmax,
  nlambda = nlambda, nalpha = nalpha, ulambda = ulambda, ualpha = ualpha, eps, isd, intr, maxit,
  napenalty = integer(1),
  theta0 = double(1 * nlambda * nalpha), theta = double(pmax * nlambda * nalpha), itheta = integer(pmax * nlambda * nalpha), ntheta = integer(nlambda * nalpha), apenalty = double(nlambda * nalpha * 2), npass = integer(1), jerr = integer(1),
  PACKAGE = "sulnet"
)
output <- getoutputLA(fit, maxit, pmax, nvars, vnames)
gam <- theta2gamma(output, betaFIT)
# gridind <- .Fortran("gridindexla", nlambda, nalpha, ulambda, ualpha, ind = integer(nlambda * nalpha * 4), PACKAGE = "sulnet")
# ind <- matrix(gridind$ind, nlambda * nalpha, 4)
# storage.mode(ind) <- "integer"


lam <- getsrR(nobs, nvars, nsum, nrate, f, y, usum, urate, pf, flmins, flminr, isd, intr)
usum <- as.double(lam$usum)
urate <- as.double(sort(lam$urate))

###########################
fitzt <- suni2D(x, y, nzeta = nzeta, ntau = ntau)
fitzt$beta
uniFit(x, y)$beta
fitzt$gamma[[ntau]]
fitzt$gamma0[, ntau]
coef(fitzt, penalty2 = fitzt$tau[ntau])
object <- fitzt
penalty1 <- object$zeta[c(5, 90)] * 0.5
penalty2 <- object$tau[1:5] * 0.5
coef(fitzt, penalty1, penalty2)
mean(abs(as.vector(predict(fitzt, x, penalty1 = fitzt$zeta[nzeta], penalty2 = fitzt$tau[ntau])[[1]]) - y))


coef1 <- coef(fitzt)
suni2D.object <- fitzt
pred.loss <- "loss"
N <- nrow(x)
nfolds <- 5
foldid <- sample(rep(seq(nfolds), length = N))
outlist <- as.list(seq(nfolds))
for (i in seq(nfolds)) {
  which <- foldid == i
  y_sub <- y[!which]
  args <- c(
    list(x = x[!which, , drop = FALSE], y = y_sub, method = method),
    setNames(pen, strsplit(method, "_")[[1]]) # , # Dynamic arguments
    # list(...) # Additional arguments
  )
  outlist[[i]] <- do.call(suni2D, args)
}

zeta <- NULL
tau <- NULL
lambda <- NULL
alpha <- NULL
sum = NULL
rate = NULL
beta = NULL
N <- nrow(x)
y <- drop(y)
method = "zeta_tau"
lambda2 = 0
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
lam2 <- as.double(lambda2)
pf <- as.double(pf)
pf2 <- as.double(pf2)
pf3 <- as.double(pf3)
isd <- as.integer(standardize)
intr <- as.integer(intercept)
eps <- as.double(eps)
pmax <- as.integer(pmax)
maxit <- as.integer(maxit)
pmax = nvars
standardize = TRUE
intercept = TRUE
eps <- 1e-08
maxit = 1e+05
family = "gaussian"
jd <- as.integer(0)
lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04)
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


suni2D.object <- suni2D(x, y,
  zeta = zeta, tau = tau,
  lambda = lambda, alpha = alpha, sum = sum, rate = rate,
  method = method # , nzeta = nzeta, ntau = ntau
)
method <- suni2D.object$use
pen <- getpenalty(suni2D.object)
pen1 = pen$pen1
pen2 = pen$pen2
nfolds <- 5
pred.loss <- "loss"
foldid <- sample(rep(seq(nfolds), length = N))
outlist <- as.list(seq(nfolds))
for (i in seq(nfolds)) {
  which <- foldid == i
  y_sub <- y[!which]
  args <- c(
    list(x = x[!which, , drop = FALSE], y = y_sub, method = method),
    setNames(pen, strsplit(method, "_")[[1]])
  )
  outlist[[i]] <- do.call(suni2D, args)
}
fun <- paste("cv", class(suni2D.object)[[2]], sep = ".")
cvstuff <- do.call(fun, list(
  outlist, method, pen, x, y, foldid,
  pred.loss
))
cvm <- cvstuff$cvm
cvsd <- cvstuff$cvsd
cvname <- cvstuff$name
lamin <- getmin2D(pen$pen1, pen$pen2, cvm, cvsd, method)
out <- list(setNames(pen, strsplit(method, "_")[[1]]),
  cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
  cvlower = cvm - cvsd, name = cvname,
  suni2D.fit = suni2D.object
)


cvfitzt <- cv.suni2D(x, y)
cvfitzt$cv.min
cvfitzt$cv.1se














template <- sulnet2D(x, y, method = "septhresh", nlambda = 10, lamPos = c(3, 2, 1), alpha = NULL)
template$b0
template$beta
