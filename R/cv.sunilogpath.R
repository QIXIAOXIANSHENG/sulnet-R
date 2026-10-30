cv.sunilogpath <- function(outlist, method, pen, x, y, foldid, pred.loss) {
  typenames <- c(
    misclass = "Misclassification Error",
    loss = "Margin Based Loss"
  )
  if (pred.loss == "default") {
    pred.loss <- "loss"
  }
  if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
    warning("Only 'misclass' and 'loss' available for logistic regression; 'loss' used")
    pred.loss <- "loss"
  }
  prob_min <- 1e-05
  fmax <- log(1 / prob_min - 1)
  fmin <- -fmax
  y <- as.factor(y)
  y <- c(-1, 1)[as.numeric(y)]
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
  for (i in seq_along(predmatlist)) {
    predmatlist[[i]] <- pmin(pmax(predmatlist[[i]], fmin), fmax)
  }
  cvraw <- lapply(predmatlist, function(predmat) {
    switch(pred.loss,
      loss = 2 * log(1 + exp(-y * predmat)),
      misclass = (y != ifelse(predmat > 0, 1, -1))
    )
  })
  N <- lapply(predmatlist, function(predmat) {
    length(y) - apply(is.na(predmat), 2, sum)
  })
  ## A cell missing from any fold gets an NA loss instead of an average over
  ## the folds that reached it, so it cannot be compared with complete cells.
  cvm <- lapply(cvraw, function(x) {
    apply(x, 2, mean)
  })
  cvsd <- lapply(seq_along(cvraw), function(x) {
    sqrt(apply(scale(cvraw[[x]], cvm[[x]], FALSE)^2, 2, mean) / (N[[x]] - 1))
  })
  cvm <- do.call(rbind, cvm)
  cvsd <- do.call(rbind, cvsd)
  list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
