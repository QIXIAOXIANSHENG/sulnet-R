#' @method predict sunipath
#' @export
predict.sunilogpath <- function(
  object, newx, penalty1 = NULL,
  penalty2 = NULL, type = c("class", "link"), ...
) {
  type <- match.arg(type)
  gammas <- coef(object, penalty1, penalty2)
  nnewx <- as.matrix(cbind2(1, newx))
  nfit <- lapply(seq_along(gammas), function(x) {
    as.matrix(nnewx %*% gammas[[x]])
  })
  nfit <- lapply(nfit, function(x) {
    switch(type,
      link = x,
      class = ifelse(x > 0, 1, -1)
    )
  })
  names(nfit) <- names(gammas)
  nfit
}
