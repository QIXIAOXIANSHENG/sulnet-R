#' Score held-out predictions for sunipath fits
#'
#' Aggregate squared-error validation losses from already fitted training-fold
#' models over a two-dimensional penalty grid.
#'
#' @details Missing predictions from truncated paths are omitted from the loss
#'   averages and standard-error calculations.
#'
#' @param outlist List of fitted training-fold models, indexed by fold number.
#' @param method Character string naming the two penalty coordinates, separated
#'   by an underscore (for example, \code{"lambda_alpha"}).
#' @param pen List containing the two full-data penalty vectors. The first
#'   vector defines output columns and the second defines output rows.
#' @param x Full predictor matrix in the same observation order as the fold
#'   labels.
#' @param y Numeric response for all observations.
#' @param foldid Length-n vector of fold labels, numbered consecutively from 1.
#' @param pred.loss Loss selector; \code{"default"} selects squared error and
#'   unsupported values fall back to \code{"loss"} with a warning.
#'
#' @return A list containing \code{cvm}, \code{cvsd}, and the loss label
#'   \code{name}. The first two components are matrices with second-coordinate
#'   values in rows and first-coordinate values in columns.
#'
#' @export
cv.sunipath <- function(outlist, method, pen, x, y, foldid, pred.loss) {
  typenames <- c(
    misclass = "Misclassification Error",
    loss = "Least Square Loss"
  )
  if (pred.loss == "default") {
    pred.loss <- "loss"
  }
  if (!match(pred.loss, c("loss"), FALSE)) {
    warning("Only 'loss' available for least squares regression; 'loss' used")
    pred.loss <- "loss"
  }
  y <- as.double(y)
  nfolds <- max(foldid)
  pen1name <- strsplit(method, "_")[[1]][1]

  predmatlist <- vector("list", length = length(pen$pen2))
  predmat <- matrix(NA, length(y), length(pen$pen1))
  for (i in seq_along(predmatlist)) {
    predmatlist[[i]] <- predmat
  }
  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]][[pen1name]])
    for (a in seq_along(predmatlist)) {
      predmatlist[[a]][which, seq(nlami)] <- preds[[a]]
    }
    # predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- lapply(predmatlist, function(predmat) {
    (y - predmat)^2
  })
  N <- lapply(predmatlist, function(predmat) {
    length(y) - apply(is.na(predmat), 2, sum)
  })
  cvm <- lapply(cvraw, function(x) {
    apply(x, 2, mean, na.rm = TRUE)
  })
  cvsd <- lapply(seq_along(cvraw), function(x) {
    sqrt(apply(scale(cvraw[[x]], cvm[[x]], FALSE)^2,
      2, mean,
      na.rm = TRUE
    ) / (N[[x]] - 1))
  })
  cvm <- do.call(rbind, cvm)
  cvsd <- do.call(rbind, cvsd)
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
