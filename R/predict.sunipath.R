##' @export
predict.sunipath <- function(
  object, newx,
  penalty1 = NULL, penalty2 = NULL, type = c("link"), ...
) {
  type <- match.arg(type)
  gammas <- coef(object, penalty1, penalty2)
  nnewx <- as.matrix(cbind2(1, newx))
  nfit <- lapply(seq_along(gammas), function(x) {
    as.matrix(nnewx %*% gammas[[x]])
  })
  names(nfit) <- names(gammas)
  nfit
}
