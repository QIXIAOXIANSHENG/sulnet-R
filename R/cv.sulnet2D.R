#' Cross-validate a sulnet2D fit
#'
#' Refit the model on training folds and summarize held-out prediction loss over
#'   its penalty grid.
#'
#' @details Fit once on all observations, then once per training fold. Ties at
#'   the minimum favor the larger second-coordinate value, then the larger
#'   first-coordinate value. The one-standard-error rule searches the eligible
#'   region at least as strongly penalized in both coordinates, favoring the
#'   larger second coordinate first. At least three folds are required.
#'   The first grid coordinate is the lambda sequence; the second coordinate
#'   is either alpha or lamPos, depending on the model fit.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}.
#' @param lambda Optional lambda grid passed to the full-data fit; the resulting
#'   fitted grid is reused in each fold.
#' @param alpha Alpha grid passed to \code{sulnet2D}, including that
#'   function's rescaling rules. Use \code{NULL} for lamPos-based
#'   fitting. When set to \code{NULL}, include a \code{lamPos} grid in the
#'   \code{...} arguments to \code{sulnet2D}.
#' @param pred.loss Prediction loss selector. When omitted, the scoring routine
#'   chooses its default. Gaussian paths use squared error, logistic paths use
#'   \code{"loss"} or \code{"misclass"}, and Cox scoring uses deviance.
#' @param nfolds Number of folds when \code{foldid} is omitted; default 5.
#' @param foldid Optional length-n vector of fold labels in 1 through the number
#'   of folds. Its maximum replaces \code{nfolds}. Otherwise labels are randomly
#'   sampled with approximately equal fold sizes.
#' @param ... Additional arguments passed to \code{sulnet2D} for the full-data
#'   and training-fold fits.
#'
#' @return A list of class \code{cv.sulnet2D} with the penalty grid, \code{cvm}
#'   (mean held-out loss), \code{cvsd} (standard error), \code{cvupper},
#'   \code{cvlower}, \code{name}, \code{sulnet2D.fit}, and the selected pairs
#'   \code{cv.min} and \code{cv.1se}. Loss matrices have second-coordinate
#'   values in rows and first-coordinate values in columns. Also includes
#'   \code{nzero}.
#' @examples
#' set.seed(1)
#' nobs <- as.integer(100)
#' nvars <- as.integer(20)
#' x <- matrix(rnorm(nobs * nvars), nobs, nvars)
#' y <- rnorm(nobs)
#' cvfit <- cv.sulnet2D(x, y, nfolds = 5)
#'
#' @export
cv.sulnet2D <- function(x, y, lambda = NULL, alpha = seq(0.01, 0.99, length.out = 21),
                        pred.loss = c("misclass", "loss"),
                        nfolds = 5, foldid, ...) {
  if (missing(pred.loss)) {
    pred.loss <- "default"
  } else {
    pred.loss <- match.arg(pred.loss)
  }

  N <- nrow(x)
  ## Fit the model once to get dimensions etc of output
  y <- drop(y)
  sulnet2D.object <- sulnet2D(x, y,
    lambda = lambda, alpha = alpha,
    ...
  )
  lambda <- sulnet2D.object$lambda

  if (sulnet2D.object$use_alpha) {
    alpha <- sulnet2D.object$alpha
  } else {
    alpha <- sulnet2D.object$lamPos
  }

  ## predict -> coef
  nz <- sapply(
    coef(sulnet2D.object, type = "nonzero"),
    function(i) {
      lapply(i, function(x) {
        sapply(x, length)
      })
    }
  )
  if (missing(foldid)) {
    foldid <- sample(rep(seq(nfolds), length = N))
  } else {
    nfolds <- max(foldid)
  }
  if (nfolds < 3) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  outlist <- as.list(seq(nfolds))
  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    if (sulnet2D.object$use_alpha) {
      outlist[[i]] <- sulnet2D(
        x = x[!which, , drop = FALSE],
        y = y_sub, lambda = lambda, alpha = alpha,
        ...
      )
    } else {
      outlist[[i]] <- sulnet2D(
        x = x[!which, , drop = FALSE],
        y = y_sub, lambda = lambda, alpha = NULL,
        ...
      )
    }
  }
  ## What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(sulnet2D.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(
    outlist, lambda, alpha, x, y, foldid,
    pred.loss, sulnet2D.object
  ))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  if (sulnet2D.object$use_alpha) {
    out <- list(
      lambda = lambda, alpha = alpha, cvm = cvm, cvsd = cvsd,
      cvupper = cvm + cvsd, cvlower = cvm - cvsd,
      nzero = nz,
      name = cvname, sulnet2D.fit = sulnet2D.object
    )
  } else {
    out <- list(
      lambda = lambda, lamPos = alpha, cvm = cvm, cvsd = cvsd,
      cvupper = cvm + cvsd, cvlower = cvm - cvsd,
      nzero = nz,
      name = cvname, sulnet2D.fit = sulnet2D.object
    )
  }
  method <- ifelse(sulnet2D.object$use_alpha, "lambda_alpha", "lambda_lamPos")

  lamin <- getmin2D(lambda, alpha, cvm, cvsd, method)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sulnet2D"
  obj
  # cvstuff
  # out
}
