#' Fit legacy two-dimensional penalized regression paths
#'
#' Fit a collection of lambda paths indexed by alpha or a separate
#'   negative-coefficient penalty.
#'
#' @details For the Gaussian \code{suni_2} branch, univariate slopes and
#'   intercepts construct the transformed predictors. The solver coefficients
#'   are reconstructed on the original scale before storage in \code{beta} and
#'   \code{b0}; these names therefore differ from the initial weights named beta
#'   in the newer grid API. With non-NULL alpha, this branch uses an L1
#'   threshold \code{lambda * (1 - alpha) * pf} and adds \code{lambda * alpha}
#'   to the negative threshold. Other branches have their own native
#'   implementations. Not every accepted method/family pair has a fitting
#'   branch. A repeated multi-value alpha grid has zero range in the rescaling
#'   step.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}.
#' @param nlambda Number of lambda values to generate when \code{lambda} is
#'   \code{NULL}; default 100.
#' @param method Legacy path implementation. Gaussian branches defined in
#'   \code{R/} include \code{"suni_2"}, \code{"septhresh"}.
#'  Binomial fitting has a \code{"septhresh"} branch; Cox fitting
#'   has \code{"suni_2"} and \code{"septhresh"} branches.
#' @param family Response family: \code{"gaussian"}, \code{"binomial"}, or
#'   \code{"cox"}. Cox responses are right-censored survival responses.
#'   Available method branches depend on the family.
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
#' @param standardize Whether to scale the predictors used by the penalized
#'   solver; default FALSE.
#' @param intercept Whether to fit an intercept in the penalized stage.
#'   Univariate preprocessing, when used, computes its own predictor-specific
#'   intercepts.
#' @param eps Convergence tolerance passed to the native coordinate solver;
#'   default 1e-8.
#' @param maxit Base coordinate-pass limit, default 1e5. This wrapper multiplies
#'   it by the number of lambda values for its current method choices.
#' @param lamPos Negative-penalty grid, sorted increasingly and used when
#'   \code{alpha = NULL}; default 0.1.
#' @param loo Whether univariate preprocessing uses leave-one-out fitted values;
#'   ignored when that preprocessing is bypassed.
#' @param alpha Mixing grid, defaulting to 21 values from 0.01 to 0.99. Absolute
#'   values are sorted. A vector of length greater than one is rescaled to span
#'   0.01 to 0.99; a scalar is clamped to [0.5, 1]. Use \code{NULL} to select
#'   the \code{lamPos} grid.
#' @param negOnly Whether to bypass univariate preprocessing and fit directly on
#'   the supplied predictors in the univariate-guided path routines.
#'
#' @return A list inheriting from \code{sulnet2D} and the selected path class.
#'   Components may include:
#'   \describe{
#'     \item{beta}{Coefficient paths. For univariate-guided Gaussian fits, a
#'       list of p by fitted-path-length matrices on the original-predictor
#'       scale, one for each second-coordinate value.}
#'     \item{b0}{Intercept paths on the original-response scale, one for each
#'       second-coordinate value. Cox fits omit this component.}
#'     \item{lambda}{Penalty values corresponding to the fitted paths.}
#'     \item{lambdamat}{Branch-specific lambda vectors, when separate paths use
#'       different penalty values.}
#'     \item{alpha, lamPos}{The second-coordinate grid; \code{alpha} is returned
#'       for alpha-based fits and \code{lamPos} for negative-penalty fits.}
#'     \item{use_alpha}{Whether the alpha-based parameterization was used.}
#'     \item{df}{Numbers of nonzero coefficients along the fitted paths.}
#'     \item{dim}{Dimensions of the fitted coefficient paths.}
#'     \item{npasses, jerr}{Solver diagnostics giving the number of coordinate
#'       passes and the integer error code, respectively.}
#'     \item{fbeta, fb0}{Intermediate coefficient and intercept paths produced
#'       when univariate preprocessing is used.}
#'     \item{univariate.fit}{The first-stage fit used for univariate
#'       preprocessing.}
#'     \item{call}{The matched function call.}
#'   }
#' @export
#' @useDynLib sulnet, .registration = TRUE
sulnet2D <- function(x, y, nlambda = 100,
                     method = c(
                       "suni_2", "septhresh"
                     ),
                     family = c("gaussian", "binomial", "cox"),
                     lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
                     lambda = NULL, lambda2 = 0, pf = rep(1, nvars),
                     pf2 = rep(1, nvars), exclude, dfmax = nvars + 1,
                     pmax = min(dfmax * 1.2, nvars), standardize = FALSE,
                     intercept = TRUE, eps = 1e-08, maxit = 1e+05, lamPos = 0.1,
                     loo = TRUE, alpha = seq(0.01, 0.99, length.out = 21),
                     negOnly = FALSE) {
  ################################################################################
  ## data setup
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
  lam2 <- as.double(lambda2)
  pf <- as.double(pf)
  pf2 <- as.double(pf2)
  isd <- as.integer(standardize)
  intr <- as.integer(intercept)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  lamPos <- as.double(sort(lamPos))
  ignore_lamPos <- FALSE
  if (!is.null(alpha)) {
    alpha <- as.double(sort(abs(alpha)))
    ignore_lamPos <- TRUE
    if (length(alpha) > 1) {
      alpha <- (alpha - min(alpha)) / diff(range(alpha)) * (0.99 - 0.01) + 0.01
    } else {
      alpha <- min(1, max(alpha, 0.5))
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
  ################################################################################
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1) {
      stop("lambda.factor should be less than 1")
    }
    flmin <- as.double(lambda.factor)
    ulam <- double(nlam)
  } else {
    ## flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0)) {
      stop("lambdas should be non-negative")
    }
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  maxit <- ifelse(grepl("cold", method), as.integer(maxit),
    as.integer(maxit * nlam)
  )
  if (grepl("adasuni", method)) {
    if (lambda2 == 0) lambda2 <- 0.1
  }
  ################################################################################
  if (family == "binomial") {
    fit <- switch(method,
      septhresh = logsunipath(
        x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
        pmax, jd, pf, pf2, maxit, lam2, lamPos, loo, negOnly,
        nobs, nvars, vnames, alpha, ignore_lamPos
      )
    )
  } else if (family == "gaussian") {
    fit <- switch(method,
      suni_2 = sunipath_2(
        x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
        pmax, jd, pf, pf2, maxit, lam2, lamPos, loo, negOnly,
        nobs, nvars, vnames, alpha, ignore_lamPos
      ),
      septhresh = septhreshpath(
        x, y, nlam, flmin, ulam, isd, intr, eps, dfmax,
        pmax, jd, pf, pf2, maxit, lam2, lamPos, loo, negOnly,
        nobs, nvars, vnames, alpha, ignore_lamPos
      )
    )
  } else if (family == "cox") {
    cl <- rbind(-9e30, 9e30)
    storage.mode(cl) <- "double"
    if (!inherits(x, "sparseMatrix")) {
      storage.mode(x) <- "double"
    }
    fit <- switch(method,
      suni_2 = coxnet(
        x, FALSE, y, rep(1, nobs), NULL, 1, nobs, nvars,
        jd, pf2, cl, dfmax, as.integer(pmax), nlam, flmin, ulam, eps,
        isd, vnames, maxit
      ),
      septhresh = coxsuninet(
        x, FALSE, y, rep(1, nobs), NULL, 1, nobs,
        nvars, jd, pf2, dfmax, as.integer(pmax), nlam, flmin, ulam, eps,
        isd, vnames, maxit, lamPos, alpha, ignore_lamPos, negOnly, loo
      )
    )
  }

  if (is.null(lambda) && family != "cox") {
    fit$lambda <- lamfix(fit$lambda)
  }
  fit$call <- this.call
  ################################################################################
  class(fit) <- c("sulnet2D", class(fit))
  fit
}
