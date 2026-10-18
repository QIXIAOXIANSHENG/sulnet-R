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
  method <- suni2D.object$use
  pen <- getpenalty(suni2D.object)
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
  method <- object$use
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
