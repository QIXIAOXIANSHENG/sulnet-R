#' Cross-validate a sulnet fit
#'
#' Refit the model on training folds and summarize held-out prediction loss over
#'   its penalty grid.
#'
#' @details Fit once on all observations, then once per training fold. The
#'   minimum-error selection breaks ties toward the largest lambda. The
#'   one-standard-error choice uses the largest eligible lambda. At least three
#'   folds are required.
#'
#' @param x Predictor matrix with observations in rows and predictors in
#'   columns; coerced to a dense matrix by the fitting function.
#' @param y Response with one observation per row of \code{x}.
#' @param lambda Optional lambda grid passed to the full-data fit; the resulting
#'   fitted grid is reused in each fold.
#' @param pred.loss Prediction loss selector. The least-squares scorer uses
#'   squared error when omitted or set to \code{"loss"}. It warns and uses
#'   squared error if \code{"misclass"} is requested.
#' @param nfolds Number of folds when \code{foldid} is omitted; default 5.
#' @param foldid Optional length-n vector of fold labels in 1 through the number
#'   of folds. Its maximum replaces \code{nfolds}. Otherwise labels are randomly
#'   sampled with approximately equal fold sizes.
#' @param ... Additional arguments passed to \code{sulnet} for the full-data and
#'   training-fold fits.
#'
#' @return A list of class \code{cv.sulnet} with the penalty grid, \code{cvm}
#'   (mean held-out loss), \code{cvsd} (standard error), \code{cvupper},
#'   \code{cvlo}, \code{nzero}, \code{name}, \code{sulnet.fit},
#'   \code{lambda.min}, and \code{lambda.1se}. Loss summaries are vectors
#'   indexed by lambda.
#' @export
cv.sulnet <- function(x, y, lambda = NULL, pred.loss = c("misclass", "loss"),
                      nfolds = 5, foldid, ...) {
  if (missing(pred.loss)) {
    pred.loss <- "default"
  } else {
    pred.loss <- match.arg(pred.loss)
  }
  N <- nrow(x)
  ## Fit the model once to get dimensions etc of output
  y <- drop(y)
  sulnet.object <- sulnet(x, y,
    lambda = lambda,
    ...
  )
  lambda <- sulnet.object$lambda
  ## predict -> coef
  nz <- sapply(coef(sulnet.object, type = "nonzero"), length)
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
    outlist[[i]] <- sulnet(
      x = x[!which, , drop = FALSE],
      y = y_sub, lambda = lambda,
      ...
    )
  }
  ## What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(sulnet.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(
    outlist, lambda, x, y, foldid,
    pred.loss
  ))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(
    lambda = lambda, cvm = cvm, cvsd = cvsd,
    cvupper = cvm + cvsd, cvlo = cvm - cvsd,
    nzero = nz, name = cvname, sulnet.fit = sulnet.object
  )
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.sulnet"
  obj
}
