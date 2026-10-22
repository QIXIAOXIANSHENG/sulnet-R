#' Compute univariate fits using native routines
#'
#' Fit each predictor separately and return full-data coefficients together with
#'   ordinary or leave-one-out fitted values.
#'
#' @details The native \code{loofit} routine standardizes each nonconstant
#'   Gaussian predictor, fits its univariate regression and rescales the
#'   coefficients. Thus the slope is covariance divided by predictor variance.
#'   The \code{loofit_binom} routine uses iteratively reweighted univariate
#'   fitting and a leave-one-out adjustment when requested.
#'
#' @param x Double-precision numeric n by p matrix; this wrapper does not coerce
#'   its storage mode.
#' @param y Length-n double-precision response for Gaussian fitting. For
#'   binomial fitting, a two-class response; factor-level order determines the
#'   internal 0/1 coding.
#' @param loo Whether to return leave-one-out fitted values; default TRUE. The
#'   returned coefficients are fitted using all observations.
#' @param family Univariate model: \code{"gaussian"} or \code{"binomial"}.
#'
#' @return A list with \code{fitted_values}, an n by p matrix; \code{beta}, a
#'   length-p slope vector; \code{b0}, a length-p intercept vector; and
#'   \code{loo}. Binomial fitted values are on the linear-predictor scale.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(200), 20, 10)
#' beta <- rep(c(1, 0), 5)
#' y <- x %*% beta + rnorm(20, 0, 0.1)
#' unifit <- uniFit(x, y)
#' @export
uniFit <- function(x, y, loo = TRUE, family = c("gaussian", "binomial")) {
  # x <- as.matrix(x)
  # storage.mode(x) <- "double"
  # loo <- as.logical(loo)
  # family <- match.arg(family)
  # if (!is.null(beta_ju)) {
  #   beta_ju <- as.logical(beta_ju)
  #   if (length(beta_ju) != ncol(x)) stop("beta and beta_ju must have same length as number of columns in x")
  #   x <- x[, !beta_ju]
  # }
  family <- match.arg(family)

  np <- dim(x)
  nobs <- as.integer(np[1])
  nvars <- as.integer(np[2])
  if (nobs != length(y)) {
    stop("X and y have different number of rows in call to uniFit", call. = FALSE)
  }
  if (family == "binomial") {
    y <- factor(y)
    if (length(levels(y)) != 2) {
      stop("There should only be two classes in y")
    }
    levels(y) <- c(0, 1)
    y <- as.double(drop(as.vector(y)))
  }
  unifit <- switch(family,
    gaussian = .Fortran("loofit", nobs, nvars, x, y, loo,
      beta0 = double(nvars),
      beta = double(nvars),
      fit = double(nobs * nvars),
      PACKAGE = "sulnet"
    ),
    binomial = .Fortran("loofit_binom", nobs, nvars, x, y, loo,
      beta0 = double(nvars),
      beta = double(nvars),
      fit = double(nobs * nvars),
      PACKAGE = "sulnet"
    )
  )

  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
  return(list(
    fitted_values = f,
    beta = unifit$beta,
    b0 = unifit$beta0,
    loo = loo
  ))
}
