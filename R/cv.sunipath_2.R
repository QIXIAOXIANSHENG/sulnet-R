#' Score held-out predictions for sunipath_2 paths
#'
#' Aggregate validation losses from already fitted training-fold models.
#'
#' @details Missing predictions from truncated paths are omitted from the loss
#'   averages and standard-error calculations.
#'
#' @param outlist List of fitted training-fold models, indexed by fold number.
#' @param lambda Full-data lambda vector defining the output columns.
#' @param alpha Second-coordinate values defining output rows; these may be
#'   alpha or lamPos values.
#' @param x Full predictor matrix in the same observation order as the fold
#'   labels.
#' @param y Numeric response for all observations.
#' @param foldid Length-n vector of fold labels, numbered consecutively from 1.
#' @param pred.loss Loss selector; \code{"default"} selects squared error and
#'   unsupported values fall back to \code{"loss"} with a warning.
#' @param delta Unused compatibility argument; it may be omitted.
#'
#' @return A list containing \code{cvm}, \code{cvsd}, and the loss label
#'   \code{name}. The first two components are matrices with second-coordinate
#'   values in rows and lambda values in columns.
#' @export
cv.sunipath_2 <- function(outlist, lambda, alpha, x, y, foldid, pred.loss, delta) {
  typenames <- c(misclass = "Misclassification Error",
                 loss = "Least Square Loss")
  if (pred.loss == "default")
    pred.loss <- "loss"
  if (!match(pred.loss, c("loss"), FALSE)) {
    warning("Only 'loss' available for least squares regression; 'loss' used")
    pred.loss <- "loss"
  }
  ## Turn y into c(0,1)
  y <- as.double(y)
  nfolds <- max(foldid)

  predmatlist = vector("list",length = length(alpha))
  predmat <- matrix(NA, length(y), length(lambda))
  for (i in seq_along(predmatlist)) {
    predmatlist[[i]] <- predmat
  }

  nlams <- double(nfolds)
  for (i in seq(nfolds)) {
    which <- foldid == i
    fitobj <- outlist[[i]]
    preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
    nlami <- length(outlist[[i]]$lambda)
    for(a in seq_along(predmatlist)){
      predmatlist[[a]][which, seq(nlami)] <- preds[[a]]
    }
    #predmat[which, seq(nlami)] <- preds
    nlams[i] <- nlami
  }
  cvraw <- lapply(predmatlist, function(predmat){(y-predmat)^2})
  N <- lapply(predmatlist, function(predmat){length(y) - apply(is.na(predmat), 2, sum)})
  cvm <- lapply(cvraw,function(x){apply(x, 2, mean, na.rm = TRUE)})
  cvsd <- lapply(seq_along(cvraw),function(x){
    sqrt(apply(scale(cvraw[[x]], cvm[[x]], FALSE)^2,
                     2, mean, na.rm = TRUE)/(N[[x]] - 1))})
  cvm <- do.call(rbind, cvm)
  cvsd <- do.call(rbind, cvsd)
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
