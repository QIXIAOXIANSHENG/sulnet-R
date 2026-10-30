#' @method predict sunilogpath
#' @export
predict.sunilogpath <- function(
  object, newx, penalty1 = NULL,
  penalty2 = NULL, type = c("link", "class"), ...
) {
  type <- match.arg(type)
  gammas <- coef(object, penalty1, penalty2)
  nnewx <- as.matrix(cbind2(1, newx))
  nfit <- lapply(seq_along(gammas), function(x) {
    as.matrix(nnewx %*% gammas[[x]])
  })
  ## The link is positive for the second level of the training response.
  classnames <- object$classnames
  if (is.null(classnames)) classnames <- c(-1, 1)
  nfit <- lapply(nfit, function(x) {
    switch(type,
      link = x,
      class = matrix(classnames[(x > 0) + 1L],
        nrow(x), ncol(x),
        dimnames = dimnames(x)
      )
    )
  })
  names(nfit) <- names(gammas)
  nfit
}
