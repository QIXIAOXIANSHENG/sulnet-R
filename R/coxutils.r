
buildPredmat <-
    function(outlist, lambda, x, offset, foldid,  y, weights, grouped, type.measure = "deviance", ...) {
        nfolds <- max(foldid)
        if ((length(weights) / nfolds < 10) && !grouped) grouped <- TRUE
        devtrue <- type.measure == "deviance"
        cvraw <- if (devtrue) matrix(NA, nfolds, length(lambda)) else NULL
        nlambda <- length(lambda)
        predmat <- matrix(NA, nrow(x), length(lambda))
        rn <- rownames(x)
        sn <- paste("s", seq(0, length = nlambda), sep = "")
        dimnames(predmat) <- list(rn, sn)

        for (i in seq(nfolds)) {
            which <- foldid == i
            fitobj <- outlist[[i]]
            # print("before coef")
            coefmat <- coef(fitobj, s = lambda, ...)[[1]]
            # print("after coef")
            nlami <- min(ncol(coefmat), nlambda)
            if (devtrue) {
                if (grouped) {
                    plfull <- coxnet.deviance(
                        x = x, y = y, offset = offset,
                        weights = weights, beta = coefmat,
                        std.weights = FALSE
                    )
                    plminusk <- coxnet.deviance(
                        x = x[!which, ], y = y[!which, ],
                        offset = offset[!which],
                        weights = weights[!which],
                        beta = coefmat,
                        std.weights = FALSE
                    )
                    cvraw[i, seq(nlami)] <- (plfull - plminusk)[seq(nlami)]
                } else {
                    plk <- coxnet.deviance(
                        x = x[which, ], y = y[which, ],
                        offset = offset[which],
                        weights = weights[which], beta = coefmat,
                        std.weights = FALSE
                    )
                    cvraw[i, seq(nlami)] <- plk[seq(nlami)]
                }
            }
            predmat[which, seq(nlami)] <- as.matrix(x[which, ] %*% coefmat[, seq(nlami)])
            if (nlami < nlambda) {
                if (devtrue) cvraw[i, seq(from = nlami, to = nlambda)] <- cvraw[i, nlami]
                predmat[which, seq(from = nlami, to = nlambda)] <- predmat[which, nlami]
            }
        }
        if (devtrue) attr(predmat, "cvraw") <- cvraw
        predmat
    }

coxnet.deviance <- function(pred = NULL, y, x = NULL, offset = NULL, 
                            weights = NULL, std.weights = TRUE, beta = NULL) {
  y <- response.coxnet(y)
  
  # if y has 2 columns, it is right-censored data
  # if y has 3 columns, it is (start, stop] data
  # otherwise, throw error
  if (ncol(y) == 2) {
    return(coxnet.deviance0(pred = pred, y = y, x = x, offset = offset,
                            weights = weights, std.weights = std.weights,
                            beta = beta))
#   } else if (ncol(y) == 3) {
#     return(coxnet.deviance3(pred = pred, y = y, x = x, offset = offset, 
#                             weights = weights, std.weights = std.weights,
#                             beta = beta))
  } else {
    stop("Response y should have 2 columns")
  }
}

# coxnet.deviance routine for right-censored data
coxnet.deviance0 <- function(pred = NULL, y, x = NULL, offset = NULL, 
                             weights = NULL, std.weights = TRUE, beta = NULL) {
  ty <- y[, "time"]
  tevent <- y[, "status"]
  ty <- ty + (1 - tevent) * 100 * .Machine$double.eps
  nobs <- as.integer(length(ty))
  
  # hack for the case where user passes in x as sparse matrix
  if (!is.null(x) && inherits(x, "sparseMatrix")) {
    if (is.null(beta))
      stop("if x is passed, beta must also be passed")
    pred <- as.matrix(x %*% beta)
    return(coxnet.deviance0(pred = pred, y = y, offset = offset,
                            weights = weights, std.weights = std.weights))
  }
  
  # Sort out the pred, x and beta options.
  # If user provided `pred`, we let x = pred and beta = identity matrix.
  # This allows us to use the loglike Fortran routine to compute the
  # partial log likelihood.
  # In the end, only x and beta are passed to the Fortran routine.
  if (!is.null(pred)) {
    x <- as.matrix(pred)
    nvec <- ncol(x)
    beta <- diag(nvec)
    nvars <- as.integer(nvec)
  } else if (is.null(x) && is.null(beta)) {
    x <- matrix(0, nrow = nobs, ncol = 1)
    beta <- double(0)
    nvec <- 1
    nvars <- as.integer(0)
  } else if (!is.null(x) && !is.null(beta)) {
    x <- as.matrix(x)
    beta <- as.matrix(beta)
    nvec <- ncol(beta)
    nvars <- nrow(beta)
  } else {
    stop("user must pass either `pred`, or both `x` and `beta`")
  }
  storage.mode(x) <- "double"
  storage.mode(beta) <- "double"
  nvec <- as.integer(nvec)
  nvars <- as.integer(nvars)
  
  # normalize weights to sum to nobs
  if (is.null(weights))
    weights <- rep(1, nobs)
  else {
    if (std.weights) weights <- nobs * weights / sum(weights)
    weights <- as.double(weights)
  }
  
  if (is.null(offset))
    offset <- rep(0, nobs)
  else offset <- as.double(offset)
  
  # extract strata (if any)
  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  } else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) stop("length of strata != nobs")
  
  # if all in same strata, do the deviance computation
  # if not, take the sum of the strata-level deviances
  if (length(unique(strata)) == 1) {
    ### Compute saturated loglikelihood
    wd <- weights[tevent == 1]
    tyd <- ty[tevent == 1]
    if (any(duplicated(tyd))) {
      wd <- tapply(wd, tyd, sum)
    }
    wd <- wd[wd > 0]
    lsat <- -sum(wd * log(wd))
    ####
    
    fit <- .Fortran("loglike", nobs, nvars, x, ty, tevent, offset,
                    weights, nvec, beta, flog = double(nvec), jerr = integer(1),
                    PACKAGE = "sulnet")
    if (fit$jerr != 0) {
      errmsg <- jerr(fit$jerr, maxit = 0, pmax = 0, family = "cox")
      if (errmsg$fatal)
        stop(errmsg$msg, call. = FALSE)
      else warning(errmsg$msg, call. = FALSE)
    }
    return(2 * (lsat - fit$flog))
  } else {
    # more than one strata provided: return the sum of strata-level deviances
    tot_dev <- 0
    for (i in unique(strata)) {
      ii <- which(strata == i)
      tot_dev <- tot_dev + 
        coxnet.deviance0(y = y[ii, , drop = FALSE], x = x[ii, , drop = FALSE], 
                         beta = beta, offset = offset[ii], 
                         weights = weights[ii], std.weights = FALSE)
    }
    return(tot_dev)
  }
}

cvcox <- function(predmat,y,weights,foldid,grouped){
    ## Note, all the work got done in buildPredmat.coxnetlist for deviance; a special case
    nfolds = max(foldid)

    if ((length(weights)/nfolds < 10) && !grouped) {
        warning("Option grouped=TRUE enforced for cv.coxnet, since < 10 observations per fold",
                call. = FALSE)
    }

    
    cvraw=attr(predmat,"cvraw")
    N = nfolds - apply(is.na(cvraw), 2, sum)
    weights = as.vector(tapply(weights, foldid, sum))
    cvraw = cvraw / weights
    
    list(cvraw=cvraw,weights=weights,N=N,grouped=FALSE)

}

