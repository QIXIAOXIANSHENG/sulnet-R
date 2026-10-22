#' Fit a penalized least-squares coefficient path
#'
#' Prepare dense regression data and fit the selected one-dimensional penalty
#'   path.
#'
#' @details The implemented \code{ls} branch calls the native \code{lslassoNET}
#'   solver with L1 and quadratic penalties. The transformation-based legacy
#'   branches have different coefficient parameterizations; their names do not
#'   imply support by the newer two-dimensional \code{sunipath} methods.
#'   Automatically generated paths may stop early. The default
#'   \code{standardize} is TRUE and the default \code{maxit} is 1e6.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}.
#' @param nlambda Number of lambda values to generate when \code{lambda} is
#'   \code{NULL}; default 100.
#' @param method Path implementation: \code{"ls"} (default), \code{"suni"}, or
#'   \code{"test"}. The latter two currently refer to legacy helpers that are
#'   not defined in \code{R/}.
#' @param lambda.factor Ratio controlling the lower end of an automatically
#'   generated penalty path; defaults to 0.01 when observations are fewer than
#'   predictors, and 0.0001 otherwise.
#' @param lambda Optional nonnegative lambda values, sorted in decreasing order.
#'   Supplying these overrides \code{nlambda} and \code{lambda.factor}.
#' @param lambda2 Nonnegative quadratic-penalty tuning parameter passed to the
#'   solver as \code{lam2}; default zero.
#' @param pf Length-p vector of L1 penalty factors, defaulting to ones. The
#'   native least-squares solvers truncate negative entries to zero and require
#'   at least one positive entry.
#' @param pf2 Length-p vector of quadratic penalty factors, defaulting to ones.
#'   The native least-squares solvers truncate negative entries to zero and
#'   require at least one positive entry, even when \code{lambda2 = 0}.
#' @param exclude Optional vector of one-based predictor indices to exclude from
#'   the penalized fit.
#' @param dfmax Maximum number of nonzero coefficients used as a path stopping
#'   criterion; default p + 1.
#' @param pmax Maximum active-set size allocated by the solver; default
#'   \code{min(dfmax * 1.2, p)}.
#' @param standardize Whether the native solver scales its input predictors.
#'   Returned coefficients are rescaled by the solver.
#' @param intercept Whether to fit an intercept in the penalized stage.
#'   Univariate preprocessing, when used, computes its own predictor-specific
#'   intercepts.
#' @param eps Convergence tolerance passed to the native coordinate solver;
#'   default 1e-8.
#' @param maxit Coordinate-pass limit passed to the solver.
#'
#' @return A list inheriting from \code{sulnet} and the path class.
#'   \describe{
#'     \item{b0}{Intercepts along the fitted path.}
#'     \item{beta}{A p by fitted-path-length sparse coefficient matrix on the
#'       input-predictor scale.}
#'     \item{df}{Number of nonzero coefficients at each fitted lambda value.}
#'     \item{dim}{Dimensions of the coefficient matrix.}
#'     \item{lambda}{Penalty values corresponding to the fitted path.}
#'     \item{npasses}{Number of coordinate-descent passes made by the solver.}
#'     \item{jerr}{Integer error code returned by the solver.}
#'     \item{call}{The matched function call.}
#'   }
#' @export
#' @useDynLib sulnet, .registration = TRUE
sulnet <- function(x, y, nlambda = 100,
                   method = c("ls"),
                   lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
                   lambda = NULL, lambda2 = 0, pf = rep(1, nvars),
                   pf2 = rep(1, nvars), exclude, dfmax = nvars + 1,
                   pmax = min(dfmax * 1.2, nvars), standardize = TRUE,
                   intercept = TRUE, eps = 1e-08, maxit = 1e+06) {
  ################################################################################
  ## data setup
  method <- match.arg(method)
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
  if (NCOL(y) > 1L) stop("Multivariate response is not supported now")
  ################################################################################
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
  maxit <- as.integer(maxit)
  lam2 <- as.double(lambda2)
  pf <- as.double(pf)
  pf2 <- as.double(pf2)
  isd <- as.integer(standardize)
  intr <- as.integer(intercept)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  if (!missing(exclude)) {
    jd <- match(exclude, seq(nvars), 0)
    if (!all(jd > 0)) {
      stop("Some excluded variables out of range")
    }
    jd <- as.integer(c(length(jd), jd))
  } else {
    jd <- as.integer(0)
  }
  ################################################################################
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) {
      stop("lambda.factor should be less than 1")
    }
    flmin <- as.double(lambda.factor)
    ulam <- double(1)
  } else {
    ## flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) {
      stop("lambdas should be non-negative")
    }
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ################################################################################
  fit <- switch(method,
    ls = lspath(
      x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
      pmax, jd, pf, pf2, maxit, lam2, nobs, nvars, vnames
    )
  )
  if (is.null(lambda)) {
    fit$lambda <- lamfix(fit$lambda)
  }
  fit$call <- this.call
  ################################################################################
  class(fit) <- c("sulnet", class(fit))
  fit
}
