##' @export
cv.coxsunipath <- function(outlist, lambda, alpha, x, y, foldid, pred.loss, sulnet2D.obj) {
    weights <- rep(1, nrow(x))
    typenames <- c(deviance = "Deviance")
    if (pred.loss == "default") {
        pred.loss <- "deviance"
    }
    if (!match(pred.loss, c("deviance"), FALSE)) {
        warning("Only 'deviance' available for logistic regression; 'deviance' used")
        pred.loss <- "deviance"
    }
    grouped <- TRUE
    if (is.null(alpha)) alpha <- outlist[[1]]$lamPos
    cvm <- matrix(NA, length(alpha), length(lambda))
    cvsd <- matrix(NA, length(alpha), length(lambda))
    for (a in seq_along(alpha)) {
        lamlen <- sulnet2D.obj$dim[[a]][2]
        alpha_lambda <- sulnet2D.obj$lambdamat[[a]]
        new_outlist <- lapply(outlist, function(x) {
            output <- list(
                beta = list(x$beta[[a]][, 1:lamlen]),
                df = list(sulnet2D.obj$df[[a]]),
                dim = list(sulnet2D.obj$dim[[a]]),
                lambdamat = list(alpha_lambda),
                use_alpha = x$use_alpha
            )
            if (x$use_alpha) {
                output[["alpha"]] <- x[["alpha"]][a]
            } else {
                output[["lamPos"]] <- x[["lamPos"]][a]
            }
            class(output) <- c("sulnet2D", "coxsunipath")
            output
        })
        # print("newout")
        # print(class(new_outlist[[1]]))
        # print(names(new_outlist[[1]]))
        # return(new_outlist[[1]])
        predmat <- buildPredmat(new_outlist, alpha_lambda, x, NULL, foldid, y = y, weights = weights, grouped = grouped, type.measure = pred.loss)
        # print("predmat")
        cvstuff <- cvcox(predmat, y, weights, foldid, grouped)
        # print("cvstuff")
        cvmv <- with(cvstuff, apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE))
        cvsdv <- with(cvstuff, sqrt(apply(scale(cvraw, cvmv, FALSE)^2, 2, weighted.mean,
            w = weights, na.rm = TRUE
        ) / (N - 1)))
        nas <- is.na(cvsdv)
        if (any(nas)) {
            alpha_lambda <- alpha_lambda[!nas]
            cvmv <- cvmv[!nas]
            cvsdv <- cvsdv[!nas]
        }
        cvmv <- c(cvmv, rep(NA, length(lambda) - lamlen))
        cvsdv <- c(cvsdv, rep(NA, length(lambda) - lamlen))
        cvm[a, ] <- cvmv
        cvsd[a, ] <- cvsdv
    }


    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}
