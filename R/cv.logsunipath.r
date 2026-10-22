#' Score held-out predictions for logsunipath paths
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
#' @param y Response for all observations; logistic responses are converted to
#'   two factor levels.
#' @param foldid Length-n vector of fold labels, numbered consecutively from 1.
#' @param pred.loss \code{"loss"} for twice the logistic negative log likelihood
#'   or \code{"misclass"} for classification error; \code{"default"} selects
#'   loss.
#' @param delta Unused compatibility argument; it may be omitted.
#'
#' @return A list containing \code{cvm}, \code{cvsd}, and the loss label
#'   \code{name}. The first two components are matrices with second-coordinate
#'   values in rows and lambda values in columns.
#' @export
cv.logsunipath <- function(outlist, lambda, alpha, x, y, foldid, pred.loss, delta){
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default")
        pred.loss <- "loss"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for logistic regression; 'loss' used")
        pred.loss <- "loss"
    }
    prob_min <- 1e-05
    fmax <- log(1/prob_min - 1)
    fmin <- -fmax
    ## Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
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
            predtemp <- preds[[a]]
            lpredtemp <- dim(predtemp)[2]
            if(lpredtemp == nlami){
                predmatlist[[a]][which, seq(nlami)] <- predtemp
            }else{
                predtemp <- cbind(predtemp, matrix(rep(predtemp[,lpredtemp],nlami - lpredtemp),ncol = nlami - lpredtemp))
                predmatlist[[a]][which, seq(nlami)] <- predtemp
            }
            
        }
        nlams[i] <- nlami
    }
    predmatlist <- lapply(predmatlist, function(predmat){pmin(pmax(predmat, fmin), fmax)})
    # predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- lapply(predmatlist, function(predmat){
        switch(pred.loss, loss = 2 * log(1 + exp(-y * predmat)),
                    misclass = (y != ifelse(predmat > 0, 1, -1)))
    })
    # cvraw <- switch(pred.loss, loss = 2 * log(1 + exp(-y * predmat)),
    #                 misclass = (y != ifelse(predmat > 0, 1, -1)))
    N <- lapply(predmatlist, function(predmat){
        length(y) - apply(is.na(predmat), 2, sum)
    })
    # N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- lapply(cvraw,function(x){apply(x, 2, mean, na.rm = TRUE)})
    cvsd <- lapply(seq_along(cvraw), function(x){
            sqrt(apply(scale(cvraw[[x]], cvm[[x]], FALSE)^2,
                2, mean, na.rm = TRUE)/(N[[x]] - 1))})
    cvm <- do.call(rbind, cvm)
    cvsd <- do.call(rbind, cvsd)
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}