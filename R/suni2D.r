#' Fit a two-dimensional univariate-guided regression grid
#'
#' Build univariate transformed predictors, solve the selected penalty grid, and
#'   reconstruct coefficients on the original predictor scale.
#'
#' @details Consider a linear regression problem \eqn{y = X\gamma + \epsilon},
#' where \eqn{X} is of dimension \eqn{n\times p}, \eqn{y} is a length-p vector,
#' and \eqn{epsilon} is a length-p random noise with zero mean.
#'   The aim is that the sign of the coefficients \eqn{\gamma_j} to have higher
#'   chance of conforming the underlying truth while maintaining high
#'   prediction accuracy. As the underlying truth might not be available or the
#'   variable itself does not carry sign-specific information, the univariate
#'   slope can act as a good prior for the sign. Users are allowed to input
#'   their own prior for certain variables that may carry sign-specific
#'   meanings. In the code, the sign priors are denoted as \code{beta}.
#'
#'   For nonconstant predictors without priors, beta is the univariate
#'   slope: correlation of \eqn{X_j} and \eqn{y} times \eqn{sd(y)/sd(x)}.
#'   Supplied priors are multiplied by the median absolute correlation over
#'   nonconstant predictors divided by the median absolute supplied
#'   prior \eqn{p_s},
#'   \eqn{\frac{\operatorname{med}_j\|\operatorname{Cor}(X_j,y)\|}
#'   {\operatorname{med}_k|p_k|}p_s}.
#'   An NA prior median is replaced by 1.
#'   It then constructs \eqn{F_{ij} = X_{ij}\beta_j +
#'   \beta_{0j}}, where \eqn{\beta_{0j} = \operatorname{mean}(y - X_j\beta_j)}.
#'   The native solver optionally centers/scales F and returns theta on its
#'   input scale. Reconstruction uses \eqn{\gamma_j = \beta_j\theta_j} and
#'   \eqn{\gamma_0 = \theta_0 + \sum_j \beta_{0j}\theta_j}; disabling the
#'   penalized intercept does not remove the first-stage offsets.
#'
#' The general loss function is \eqn{\frac{1}{2n}\lVert y-F\theta\rVert_2^2
#' +\zeta\sum_{j=1}^p \operatorname{pf1}_j\|\theta_j\|
#' +\tau\sum_{j=1}^p \operatorname{pf3}_j(\theta_j)_-
#' +\lambda_2\sum_{j=1}^p \operatorname{pf2}_j\theta_j^2,}
#' where \eqn{(t)_-=\max(-t,0)} is the negative part. Three set of penalty
#' formation are provided: zeta-tau, lambda-alpha, and sum-rate.
#'
#' In zeta-tau coordinates, the positive threshold is zeta * pf1 and the negative
#'   threshold is zeta * pf1 + tau * pf3. In lambda-alpha coordinates they are
#'   lambda * (1 - alpha) * pf1 and lambda * alpha * pf3. In sum-rate coordinates
#'   they are sum/(1 + rate) * pf1 and sum * rate/(1 + rate) * pf3. The quadratic
#'   coordinate denominator adds lambda2 * pf2. These descriptions follow the
#'   native updates without imposing a loss-normalization convention. When
#'   standardized, pf1 and pf3 are also divided by the predictor scale. The
#'   reported zeta-tau coordinates in the latter two paths assume equal side
#'   factors: the extra negative penalty is lambda * (2 * alpha - 1), or sum *
#'   (rate - 1)/(1 + rate). They are not a general conversion for unequal pf1 and
#'   pf3. In particular alpha = 0.5 (rate = 1) gives symmetric thresholds with
#'   equal side factors, and alpha = 1 gives a zero positive-side threshold.
#'   \code{alpha = 1} is not recommended in high dimension cases.
#'
#' Lambda-alpha and sum-rate traversal visits every second-coordinate value
#'   within each decreasing first-coordinate level. Warm starts copy a
#'   previously computed model with componentwise stronger side penalties when
#'   available. The returned coefficient lists use rectangular grid order;
#'   penalty.path records the actual traversal. Positive native error codes
#'   stop; negative codes warn and return partial output. Uncomputed grid cells
#'   have missing intercepts and zero sparse slope placeholders. This fitter is
#'   currently internal to the package.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}. For
#'   \code{family = "binomial"}, any two-class vector (numeric, logical,
#'   character, or factor) without missing values; it is converted to a factor
#'   and the first level is coded -1 and the second +1.
#' @param beta Optional prior weights: a two-element list of predictor indices
#'   and values, or a length-p vector. List indices are clamped to [1, p] and
#'   duplicate indices rejected. Zero priors fall back to data-derived slopes.
#' @param method Penalty parameterization: \code{"zeta_tau"},
#'   \code{"lambda_alpha"}, or \code{"sum_rate"}. Default is \code{"zeta_tau"}.
#' @param family Response family, \code{"gaussian"} or \code{"binomial"}.
#'   \code{"gaussian"} uses the least-squares grid routines; \code{"binomial"}
#'   uses the logistic first-stage fit and logistic grid routines.
#' @param nzeta Requested number of automatically generated zeta values; default
#'   100. Ignored when \code{zeta} or other \code{method} is supplied.
#' @param ntau Requested number of automatically generated tau values; default
#'   20. Ignored when \code{tau} or other \code{method} is supplied.
#' @param nlambda Requested number of automatically generated lambda values;
#' default 100. Ignored when \code{lambda} or other \code{method} is supplied.
#' @param nalpha Requested number of automatically generated alpha values;
#'   default 20. Ignored when \code{alpha} or other \code{method} is supplied.
#' @param nsum Requested number of automatically generated sum values; default
#'   100. Ignored when \code{sum} or other \code{method} is supplied.
#' @param nrate Requested number of automatically generated rate values; default
#'   20. Ignored when \code{rate} or other \code{method} is supplied.
#' @param lambda.factor Ratio controlling the lower end of an automatically
#'   generated penalty path; defaults to 0.01 when observations are fewer than
#'   predictors, and 0.0001 otherwise.
#' @param zeta Optional nonnegative zeta grid, sorted in decreasing order.
#' \code{NULL} requests an automatically generated grid of length \code{nzeta}.
#' @param tau Optional nonnegative tau grid, sorted in decreasing order.
#' \code{NULL} requests an automatically generated grid of length \code{ntau}.
#' @param lambda Optional nonnegative lambda values, sorted in decreasing order.
#' \code{NULL} requests an automatically generated grid of length \code{nlambda}.
#' @param alpha Optional increasing, deduplicated alpha grid in [0.5, 1]. alpha
#'   equal to 1 is rejected when n < p.
#' \code{NULL} requests an automatically generated grid of length \code{nalpha}.
#' @param sum Optional nonnegative sum grid, sorted in decreasing order.
#'   \code{NULL} requests an automatically generated grid of length \code{nsum}.
#' @param rate Optional increasing, deduplicated rate grid with values at least
#'   1. Infinite rates are rejected when n < p; otherwise the wrapper switches
#'   to lambda-alpha and sets lambda to \code{sum}.
#' @param lambda2 Nonnegative quadratic-penalty tuning parameter passed to the
#'   solver as \code{lam2}; default zero.
#' @param pf1 Length-p positive-side/L1 penalty factors; default ones. Negative
#'   entries are truncated to zero by the solver and at least one positive entry
#'   is required.
#' @param pf2 Length-p quadratic penalty factors; default ones. Negative entries
#'   are truncated to zero and at least one positive entry is required even if
#'   lambda2 is zero.
#' @param pf3 Length-p negative-side penalty factors; default ones. Negative
#'   entries are truncated to zero. These multiply the extra negative penalty in
#'   zeta-tau and the entire negative threshold in the other two
#'   parameterizations.
#' @param exclude Optional vector of one-based predictor indices to exclude from
#'   the penalized fit.
#' @param pmax Maximum active-set size allocated for each grid fit; default p.
#' @param standardize Whether the native solver scales its input predictors.
#'   Returned coefficients are rescaled by the solver.
#' @param intercept Whether to fit an intercept in the penalized stage.
#'   Univariate preprocessing, when used, computes its own predictor-specific
#'   intercepts.
#' @param eps Convergence tolerance passed to the native coordinate solver;
#'   default 1e-8.
#' @param maxit Coordinate-pass limit across the grid; default 1e5.
#'
#' @return A list of classes \code{suni2D} and \code{sunipath} with components:
#'   \describe{
#'     \item{gamma}{A list indexed by the second grid coordinate. Each element
#'       is a row p and column by first-grid-length sparse
#'       matrix of reconstructed coefficients
#'       on the original-predictor scale.}
#'     \item{gamma0}{A first-grid-length by second-grid-length matrix of
#'       reconstructed intercepts.}
#'     \item{df}{A first-grid-length by second-grid-length matrix giving the
#'       number of nonzero coefficients at each grid point.}
#'     \item{dim}{Dimensions of the coefficient grid.}
#'     \item{zeta, tau}{Penalty-coordinate vectors for a \code{"zeta_tau"}
#'       fit.}
#'     \item{lambda, alpha}{Penalty-coordinate vectors for a
#'       \code{"lambda_alpha"} fit.}
#'     \item{sum, rate}{Penalty-coordinate vectors for a \code{"sum_rate"}
#'       fit.}
#'     \item{method}{The penalty parameterization used for the fit.}
#'     \item{family}{The response family requested for the fit.}
#'     \item{npass}{Total number of coordinate passes made by the solver.}
#'     \item{jerr}{Integer error code returned by the solver.}
#'     \item{call}{The matched function call.}
#'     \item{penalty.path}{A matrix recording the penalty coordinates in the
#'       order traversed by the solver.}
#'     \item{beta, beta0}{The first-stage predictor weights and offsets used to
#'       construct the transformed design matrix.}
#'     \item{inner}{A list containing \code{theta}, the coefficient paths on
#'       the transformed-predictor scale; \code{theta0}, their intercepts; and
#'       \code{f}, the n by p transformed design matrix.}
#'   }
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' fit <- suni2D(x, y)
#' @export
#' @useDynLib sulnet, .registration = TRUE
suni2D <- function(
  x, y, beta = NULL,
  method = c("zeta_tau", "lambda_alpha", "sum_rate"),
  family = c("gaussian", "binomial"),
  nzeta = 100, ntau = 20, nlambda = 100, nalpha = 20, nsum = 100, nrate = 20,
  lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  zeta = NULL, tau = NULL, lambda = NULL, alpha = NULL, sum = NULL, rate = NULL,
  lambda2 = 0, pf1 = rep(1, nvars),
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
  if (family == "binomial") {
    if (anyNA(y)) {
      stop("y should not contain missing values")
    }
    y <- factor(y)
    if (nlevels(y) != 2L) {
      stop("y should have exactly two classes when family = \"binomial\"")
    }
    y <- c(-1, 1)[as.integer(y)]
  }
  ## parameter setup
  if (length(pf1) != nvars) {
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
  pf1 <- as.double(pf1)
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

  ### beta setup
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
  if (any(is.infinite(beta))) stop("beta should not contain infinite values")

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
      uzeta <- as.double(rev(sort(unique(zeta))))
      nzeta <- as.integer(length(uzeta))
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
      utau <- as.double(rev(sort(unique(tau))))
      ntau <- as.integer(length(utau))
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
        message("nsum should be at least 3, set to 3")
        nsum <- 3
      }
      flmins <- as.double(lambda.factor)
      usum <- double(nsum)
    } else {
      flmins <- as.double(1)
      if (any(sum < 0)) {
        stop("sums should be non-negative")
      }
      usum <- as.double(rev(sort(unique(sum))))
      nsum <- as.integer(length(usum))
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
          lambda <- sum
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
        message("nlambda should be at least 3, set to 3")
        nlambda <- 3
      }
      flminl <- as.double(lambda.factor)
      ulambda <- double(nlambda)
    } else {
      flminl <- as.double(1)
      if (any(lambda < 0)) {
        stop("lambdas should be non-negative")
      }
      ulambda <- as.double(rev(sort(unique(lambda))))
      nlambda <- as.integer(length(ulambda))
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

  fit <- switch(family,
    gaussian = switch(method,
      zeta_tau = suniZTpath(
        x, y, nzeta, ntau, flminz, flmint, uzeta,
        utau, isd, intr, beta, beta_ju, eps, pmax, jd, pf1, pf2, pf3,
        maxit, lam2, nobs, nvars, vnames
      ),
      lambda_alpha = suniLApath(
        x, y, nlambda, nalpha, flminl, flmina,
        ulambda, ualpha, isd, intr, beta, beta_ju, eps, pmax, jd, pf1,
        pf2, pf3, maxit, lam2, nobs, nvars, vnames
      ),
      sum_rate = suniSRpath(
        x, y, nsum, nrate, flmins, flminr, usum,
        urate, isd, intr, beta, beta_ju, eps, pmax, jd, pf1, pf2, pf3,
        maxit, lam2, nobs, nvars, vnames
      )
    ),
    binomial = switch(method,
      zeta_tau = suniZTlogpath(
        x, y, nzeta, ntau, flminz, flmint, uzeta,
        utau, isd, intr, beta, beta_ju, eps, pmax, jd, pf1, pf2, pf3,
        maxit, lam2, nobs, nvars, vnames
      ),
      lambda_alpha = suniLAlogpath(
        x, y, nlambda, nalpha, flminl, flmina,
        ulambda, ualpha, isd, intr, beta, beta_ju, eps, pmax, jd, pf1,
        pf2, pf3, maxit, lam2, nobs, nvars, vnames
      ),
      sum_rate = suniSRlogpath(
        x, y, nsum, nrate, flmins, flminr, usum,
        urate, isd, intr, beta, beta_ju, eps, pmax, jd, pf1,
        pf2, pf3,
        maxit, lam2, nobs, nvars, vnames
      )
    )
  )

  fit$call <- this.call
  fit$family <- family
  ################################################################################
  class(fit) <- c("suni2D", class(fit))
  fit
}
