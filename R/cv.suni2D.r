#' Cross-validate a suni2D fit
#'
#' Refit the model on training folds and summarize held-out prediction loss over
#'   its penalty grid.
#'
#' @details Fit once on all observations, then once per training fold. Ties at
#'   the minimum favor the larger second-coordinate value, then the larger
#'   first-coordinate value. The one-standard-error rule searches the eligible
#'   region at least as strongly penalized in both coordinates, favoring the
#'   larger second coordinate first.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}.
#' @param pred.loss Loss selector. The current grid scorer uses squared error;
#'   \code{"misclass"} produces a warning and falls back to \code{"loss"}.
#' @param nfolds Number of folds when \code{foldid} is omitted; default 5.
#' @param foldid Optional length-n vector of fold labels in 1 through the number
#'   of folds. Its maximum replaces \code{nfolds}. Otherwise labels are randomly
#'   sampled with approximately equal fold sizes.
#' @param zeta Optional zeta grid for the corresponding parameterization, passed
#'   to \code{suni2D}. The full-data fitted grid is reused in training folds.
#' @param tau Optional tau grid for the corresponding parameterization, passed
#'   to \code{suni2D}. The full-data fitted grid is reused in training folds.
#' @param lambda Optional lambda grid passed to the full-data fit; the resulting
#'   fitted grid is reused in each fold.
#' @param alpha Optional alpha grid for \code{method = "lambda_alpha"}; passed
#'   to \code{suni2D}.
#' @param sum Optional sum grid for the corresponding parameterization, passed
#'   to \code{suni2D}. The full-data fitted grid is reused in training folds.
#' @param rate Optional rate grid for the corresponding parameterization, passed
#'   to \code{suni2D}. The full-data fitted grid is reused in training folds.
#' @param method Grid coordinates: \code{"zeta_tau"}, \code{"lambda_alpha"}, or
#'   \code{"sum_rate"}.
#' @param ... Additional arguments passed to \code{suni2D} for the full-data and
#'   training-fold fits.
#'
#' @return A list of class \code{cv.suni2D} with components:
#'   \describe{
#'     \item{Penalty coordinates}{Two vectors named according to \code{method}:
#'       \code{zeta} and \code{tau}, \code{lambda} and \code{alpha}, or
#'       \code{sum} and \code{rate}.}
#'     \item{cvm}{Matrix of mean held-out losses. columns correspond to the
#'       first penalty coordinate and rows to the second.}
#'     \item{cvsd}{Matrix of estimated standard errors for \code{cvm}.}
#'     \item{cvupper}{Matrix equal to \code{cvm + cvsd}.}
#'     \item{cvlower}{Matrix equal to \code{cvm - cvsd}.}
#'     \item{name}{Name of the cross-validation loss.}
#'     \item{suni2D.fit}{The \code{suni2D} fit to all observations.}
#'     \item{cv.min}{A list of the penalty-coordinate pair attaining the
#'       minimum mean cross-validation loss.}
#'     \item{cv.1se}{A list of the penalty-coordinate pair selected by the
#'       one-standard-error rule.}
#'   }
#'   In each loss matrix, second-coordinate values index rows and
#'   first-coordinate values index columns.
#'
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(10)
#' x <- matrix(rnorm(nobs * nvars), nrow = nobs, ncol = nvars)
#' gamma <- c(rep(2, nvars / 2), rep(0, nvars / 2))
#' y <- as.double(rnorm(nobs)) + x %*% gamma
#' cvfit <- cv.suni2D(x, y, nfolds = 5, method = "lambda_alpha")
#'
#' @importFrom stats setNames
#' @export
cv.suni2D <- function(
  x, y, pred.loss = c("misclass", "loss"),
  nfolds = 5, foldid, zeta = NULL, tau = NULL,
  lambda = NULL, alpha = NULL, sum = NULL, rate = NULL,
  method = c("zeta_tau", "lambda_alpha", "sum_rate"), ...
) {
  if (missing(pred.loss)) {
    pred.loss <- "default"
  } else {
    pred.loss <- match.arg(pred.loss)
  }
  method <- match.arg(method)

  N <- nrow(x)
  y <- drop(y)
  suni2D.object <- suni2D(x, y,
    zeta = zeta, tau = tau,
    lambda = lambda, alpha = alpha, sum = sum, rate = rate,
    method = method, ...
  )
  method <- suni2D.object$method
  pen <- getpenalty(suni2D.object)
  tol <- sqrt(.Machine$double.eps)
  if (method == "lambda_alpha") {
    if (abs(pen[[2]][1] - 0.5) <= tol) {
      pen[[2]][1] <- 0.5
    }
  } else if (method == "sum_rate") {
    if (abs(pen[[2]][1] - 1) <= tol) {
      pen[[2]][1] <- 1
    }
  }

  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else {
    nfolds <- max(foldid)
  }
  outlist <- as.list(seq(nfolds))
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    args <- c(
      list(x = x[!which, , drop = FALSE], y = y_sub, method = method),
      setNames(pen, strsplit(method, "_")[[1]]), # Dynamic arguments
      list(...) # Additional arguments
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
  out <- c(
    setNames(pen, strsplit(method, "_")[[1]]),
    list(
      cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd,
      cvlower = cvm - cvsd, name = cvname,
      suni2D.fit = suni2D.object
    )
  )
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.suni2D"
  obj
}

getpenalty <- function(object) {
  method <- object$method
  penalty1 <- switch(method,
    zeta_tau = object$zeta,
    lambda_alpha = object$lambda,
    sum_rate = object$sum
  )
  penalty2 <- switch(method,
    zeta_tau = object$tau,
    lambda_alpha = object$alpha,
    sum_rate = object$rate
  )
  list(pen1 = penalty1, pen2 = penalty2)
}
