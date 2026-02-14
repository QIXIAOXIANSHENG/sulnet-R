##' @export
predict.coxsunipath <- function(
    object, newx, s = NULL,
    type = c("link", "response"),
    alpha,
    ...) {
    type <- match.arg(type)
    nbeta <- object$beta
    name <- names(nbeta)
    if (!is.null(s)) {
        nbeta <- lapply(seq_along(nbeta), function(x) {
            lambda <- object$lambdamat[[x]]
            lamlist <- lambda.interp(lambda, s)
            xuse <- nbeta[[x]]
            vnames <- dimnames(xuse)[[1]]
            dimnames(xuse) <- list(NULL, NULL)

            xuse <- xuse[, lamlist$left, drop = FALSE] %*%
                Diagonal(x = lamlist$frac) +
                xuse[, lamlist$right, drop = FALSE] %*%
                Diagonal(x = 1 - lamlist$frac)
            dimnames(xuse) <- list(vnames, paste(seq(along = s)))
            xuse
        })
    }

    if (!missing(alpha)) {
        nbeta_2 <- list()
        if (object$use_alpha) {
            for (a in alpha) {
                if (a %in% object$alpha) {
                    nbeta_2[[paste0("alpha_", a)]] <- nbeta[[paste0("alpha_", a)]]
                }
            }
        } else {
            for (a in alpha) {
                if (a %in% object$lamPos) {
                    nbeta_2[[paste0("lamPos_", a)]] <- nbeta[[paste0("lamPos_", a)]]
                }
            }
        }
        if (length(nbeta_2) != 0) {
            nbeta <- nbeta_2
        }
    }


    nfit <- lapply(seq_along(nbeta), function(x) {
        as.matrix(newx %*% nbeta[[x]])
    })

    nfit <- lapply(nfit, function(x) {
        switch(type,
            link = x,
            response = exp(x)
        )
    })
    names(nfit) <- name
    nfit
}