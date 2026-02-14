coxsuninet <- function(x, is.sparse, y, weights, offset, parm, nobs, nvars, jd, vp, ne, nx, nlam, flmin, ulam, thresh, isd, vnames, maxit, lamPos, alpha, iglamPos, negOnly, loo) {
  y <- response.coxnet(y)
  iglamPos <- as.integer(iglamPos)
  alpha <- as.double(alpha)

  if (!is.matrix(y) || !all(match(c("time", "status"), dimnames(y)[[2]], 0))) stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) as a response; a 'Surv' object suffices", call. = FALSE)
  ty <- as.double(y[, "time"])
  tevent <- as.double(y[, "status"])
  ty <- ty + (1 - tevent) * 100 * .Machine$double.eps ## ties issue
  if (any(ty <= 0)) stop("negative event times encountered;  not permitted for Cox family")
  maxit <- as.integer(maxit)
  weights <- as.double(weights)
  if (is.null(offset)) {
    offset <- ty * 0
  } else {
    storage.mode(offset) <- "double"
  }
  loo <- as.logical(loo)

  if (negOnly) {
    if (iglamPos) {
      n_alpha <- length(alpha)
      b0list <- list()
      betamat <- vector("list", length = n_alpha)
      dfmat <- vector("list", length = n_alpha)
      npassesmat <- vector("list", length = n_alpha)
      jerrmat <- vector("list", length = n_alpha)
      lambdamat <- vector("list", length = n_alpha)
      dimmat <- vector("list", length = n_alpha)

      for (a in seq_along(alpha)) {
        fit <- if (is.sparse) {
          stop("Cox model mot implemented for sparse x in glmnet", call. = FALSE)
        } else {
          .Fortran("coxsuninet",
            parm = parm, nobs, nvars, x, ty, tevent, offset, weights, jd, vp, ne, nx, nlam, flmin, ulam, thresh,
            maxit, isd, # need to get JHF to reverse these
            lmu = integer(1),
            ca = double(nx * nlam),
            ia = integer(nx),
            nin = integer(nlam),
            nulldev = double(1),
            dev = double(nlam),
            alm = double(nlam),
            nlp = integer(1),
            jerr = integer(1),
            iglamPos,
            as.double(0),
            as.double(alpha[a]), PACKAGE = "sulnet"
          )
        }
        if (fit$jerr != 0) {
          errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
          if (errmsg$fatal) {
            stop(errmsg$msg, call. = FALSE)
          } else {
            warning(errmsg$msg, call. = FALSE)
          }
        }
        outlist <- coxgetcoef(fit, nvars, nx, vnames)
        outlist <- c(outlist, list(npass = fit$nlp, jerr = fit$jerr))
        b0list[[paste0("alpha_", alpha[a])]] <- outlist$a0
        betamat[[a]] <- outlist$beta
        dfmat[[a]] <- outlist$df
        npassesmat[[a]] <- outlist$npass
        jerrmat[[a]] <- fit$jerr
        if (flmin >= 1) {
          lambdamat[[a]] <- outlist$lambda
        } else {
          lambdamat[[a]] <- lamfix(outlist$lambda)
        }
      }

      dimmat[[a]] <- outlist$dim
      alphaname <- paste0("alpha_", alpha)
      names(betamat) <- alphaname
      names(dfmat) <- alphaname
      names(npassesmat) <- alphaname
      names(jerrmat) <- alphaname
      names(lambdamat) <- alphaname
      names(dimmat) <- alphaname
      outlist <- list(
        b0 = b0list,
        beta = betamat,
        df = dfmat,
        dim = dimmat,
        lambda = lamextend(lambdamat[[1]], flmin, nlam, ulam),
        lambdamat = lambdamat,
        npasses = npassesmat,
        jerr = jerrmat,
        alpha = alpha,
        use_alpha = as.logical(iglamPos),
        negOnly = negOnly
      )
    } else {
      n_alpha <- length(lamPos)
      b0list <- list()
      betamat <- vector("list", length = n_alpha)
      dfmat <- vector("list", length = n_alpha)
      npassesmat <- vector("list", length = n_alpha)
      jerrmat <- vector("list", length = n_alpha)
      lambdamat <- vector("list", length = n_alpha)
      dimmat <- vector("list", length = n_alpha)

      for (a in seq_along(lamPos)) {
        fit <- if (is.sparse) {
          stop("Cox model mot implemented for sparse x in glmnet", call. = FALSE)
        } else {
          .Fortran("coxsuninet",
            parm = parm, nobs, nvars, x, ty, tevent, offset, weights, jd, vp, ne, nx, nlam, flmin, ulam, thresh,
            maxit, isd, # need to get JHF to reverse these
            lmu = integer(1),
            ca = double(nx * nlam),
            ia = integer(nx),
            nin = integer(nlam),
            nulldev = double(1),
            dev = double(nlam),
            alm = double(nlam),
            nlp = integer(1),
            jerr = integer(1),
            iglamPos,
            as.double(lamPos[a]),
            as.double(0), PACKAGE = "sulnet"
          )
        }
        if (fit$jerr != 0) {
          errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
          if (errmsg$fatal) {
            stop(errmsg$msg, call. = FALSE)
          } else {
            warning(errmsg$msg, call. = FALSE)
          }
        }
        outlist <- coxgetcoef(fit, nvars, nx, vnames)
        outlist <- c(outlist, list(npass = fit$nlp, jerr = fit$jerr))
        b0list[[paste0("lamPos_", lamPos[a])]] <- outlist$a0
        betamat[[a]] <- outlist$beta
        dfmat[[a]] <- outlist$df
        npassesmat[[a]] <- outlist$npass
        jerrmat[[a]] <- fit$jerr
        if (flmin >= 1) {
          lambdamat[[a]] <- outlist$lambda
        } else {
          lambdamat[[a]] <- lamfix(outlist$lambda)
        }
        dimmat[[a]] <- outlist$dim
      }
      alphaname <- paste0("lamPos_", lamPos)
      names(betamat) <- alphaname
      names(dfmat) <- alphaname
      names(npassesmat) <- alphaname
      names(jerrmat) <- alphaname
      names(lambdamat) <- alphaname
      names(dimmat) <- alphaname
      outlist <- list(
        b0 = b0list,
        beta = betamat,
        df = dfmat,
        dim = dimmat,
        lambda = lamextend(lambdamat[[1]], flmin, nlam, ulam),
        lambdamat = lambdamat,
        npasses = npassesmat,
        jerr = jerrmat,
        lamPos = lamPos,
        use_alpha = as.logical(iglamPos),
        negOnly = negOnly
      )
    }

    class(outlist) <- c("coxsunipath")
    return(outlist)
  }





  unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
    beta0 = double(nvars),
    beta = double(nvars),
    fit = double(nobs * nvars),
    PACKAGE = "sulnet"
  )
  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
  storage.mode(f) <- "double"

  if (iglamPos) {
    n_alpha <- length(alpha)
    fb0list <- list()
    fbetamat <- vector("list", length = n_alpha)
    dfmat <- vector("list", length = n_alpha)
    npassesmat <- vector("list", length = n_alpha)
    jerrmat <- vector("list", length = n_alpha)
    lambdamat <- vector("list", length = n_alpha)
    dimmat <- vector("list", length = n_alpha)

    b0list <- list()
    betamat <- vector("list", length = n_alpha)

    for (a in seq_along(alpha)) {
      fit <- if (is.sparse) {
        stop("Cox model mot implemented for sparse x in glmnet", call. = FALSE)
      } else {
        .Fortran("coxsuninet",
          parm = parm, nobs, nvars, f, ty, tevent, offset, weights, jd, vp, ne, nx, nlam, flmin, ulam, thresh,
          maxit, isd, # need to get JHF to reverse these
          lmu = integer(1),
          ca = double(nx * nlam),
          ia = integer(nx),
          nin = integer(nlam),
          nulldev = double(1),
          dev = double(nlam),
          alm = double(nlam),
          nlp = integer(1),
          jerr = integer(1),
          iglamPos,
          as.double(0),
          as.double(alpha[a]), PACKAGE = "sulnet"
        )
      }
      if (fit$jerr != 0) {
        errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
        if (errmsg$fatal) {
          stop(errmsg$msg, call. = FALSE)
        } else {
          warning(errmsg$msg, call. = FALSE)
        }
      }
      outlist <- coxgetcoef(fit, nvars, nx, vnames)
      outlist <- c(outlist, list(npass = fit$nlp, jerr = fit$jerr))


      ones <- rep(1, length(outlist$lambda))
      unibeta <- outer(unifit$beta, ones)
      unibeta0 <- outer(unifit$beta0, ones)

      beta_temp <- outlist$beta
      beta0_temp <- outlist$b0

      b0list[[paste0("alpha_", alpha[a])]] <- outlist$a0
      betamat[[a]] <- outlist$beta
      dfmat[[a]] <- outlist$df
      npassesmat[[a]] <- outlist$npass
      jerrmat[[a]] <- fit$jerr
      if (flmin >= 1) {
        lambdamat[[a]] <- outlist$lambda
      } else {
        lambdamat[[a]] <- lamfix(outlist$lambda)
      }
      dimmat[[a]] <- outlist$dim


      row_idx <- beta_temp@i + 1
      col_ptrs <- beta_temp@p
      col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
      beta_result <- beta_temp
      beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]

      betamat[[a]] <- beta_result
      b0list[[paste0("alpha_", alpha[a])]] <- beta0_temp + colSums(unibeta0 * beta_temp)
    }
    alphaname <- paste0("alpha_", alpha)
    names(fbetamat) <- alphaname
    names(betamat) <- alphaname
    names(dfmat) <- alphaname
    names(npassesmat) <- alphaname
    names(jerrmat) <- alphaname
    names(lambdamat) <- alphaname
    names(dimmat) <- alphaname

    outlist <- list(
      b0 = b0list,
      beta = betamat,
      df = dfmat,
      dim = dimmat,
      lambda = lamextend(lambdamat[[1]], flmin, nlam, ulam),
      lambdamat = lambdamat,
      npasses = npassesmat,
      jerr = jerrmat,
      alpha = alpha,
      use_alpha = as.logical(iglamPos),
      negOnly = negOnly,
      LOO = loo,
      univariate.fit = list(
        beta = unifit$beta,
        beta0 = unifit$beta0,
        fitted.values = f
      ),
      fbeta = fbetamat,
      fb0 = fb0list
    )
  } else {
    n_alpha <- length(lamPos)
    fb0list <- list()
    fbetamat <- vector("list", length = n_alpha)
    dfmat <- vector("list", length = n_alpha)
    npassesmat <- vector("list", length = n_alpha)
    jerrmat <- vector("list", length = n_alpha)
    lambdamat <- vector("list", length = n_alpha)
    dimmat <- vector("list", length = n_alpha)

    b0list <- list()
    betamat <- vector("list", length = n_alpha)

    for (a in seq_along(lamPos)) {
      fit <- if (is.sparse) {
        stop("Cox model mot implemented for sparse x in glmnet", call. = FALSE)
      } else {
        .Fortran("coxsuninet",
          parm = parm, nobs, nvars, x, ty, tevent, offset, weights, jd, vp, ne, nx, nlam, flmin, ulam, thresh,
          maxit, isd, # need to get JHF to reverse these
          lmu = integer(1),
          ca = double(nx * nlam),
          ia = integer(nx),
          nin = integer(nlam),
          nulldev = double(1),
          dev = double(nlam),
          alm = double(nlam),
          nlp = integer(1),
          jerr = integer(1),
          iglamPos,
          as.double(lamPos[a]),
          as.double(0), PACKAGE = "sulnet"
        )
      }
      if (fit$jerr != 0) {
        errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
        if (errmsg$fatal) {
          stop(errmsg$msg, call. = FALSE)
        } else {
          warning(errmsg$msg, call. = FALSE)
        }
      }
      outlist <- coxgetcoef(fit, nvars, nx, vnames)
      outlist <- c(outlist, list(npass = fit$nlp, jerr = fit$jerr))
      ones <- rep(1, length(outlist$lambda))
      unibeta <- outer(unifit$beta, ones)
      unibeta0 <- outer(unifit$beta0, ones)

      beta_temp <- outlist$beta
      beta0_temp <- outlist$a0

      fb0list[[paste0("lamPos_", lamPos[a])]] <- beta0_temp
      fbetamat[[a]] <- beta_temp
      dfmat[[a]] <- outlist$df
      npassesmat[[a]] <- outlist$npass
      jerrmat[[a]] <- fit$jerr
      if (flmin >= 1) {
        lambdamat[[a]] <- outlist$lambda
      } else {
        lambdamat[[a]] <- lamfix(outlist$lambda)
      }
      dimmat[[a]] <- outlist$dim


      row_idx <- beta_temp@i + 1
      col_ptrs <- beta_temp@p
      col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
      beta_result <- beta_temp
      beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]

      betamat[[a]] <- beta_result
      b0list[[paste0("lamPos_", lamPos[a])]] <- beta0_temp + colSums(unibeta0 * beta_temp)
    }
    alphaname <- paste0("lamPos_", lamPos)
    names(fbetamat) <- alphaname
    names(betamat) <- alphaname
    names(dfmat) <- alphaname
    names(npassesmat) <- alphaname
    names(jerrmat) <- alphaname
    names(lambdamat) <- alphaname
    names(dimmat) <- alphaname

    outlist <- list(
      b0 = b0list,
      beta = betamat,
      df = dfmat,
      dim = dimmat,
      lambda = lamextend(lambdamat[[1]], flmin, nlam, ulam),
      lambdamat = lambdamat,
      npasses = npassesmat,
      jerr = jerrmat,
      lamPos = lamPos,
      use_alpha = as.logical(iglamPos),
      negOnly = negOnly,
      LOO = loo,
      univariate.fit = list(
        beta = unifit$beta,
        beta0 = unifit$beta0,
        fitted.values = f
      ),
      fbeta = fbetamat,
      fb0 = fb0list
    )
  }


  class(outlist) <- c("coxsunipath")
  return(outlist)

  # fit <- if (is.sparse) {
  #   stop("Cox model mot implemented for sparse x in glmnet", call. = FALSE)
  # } else {
  #   .Fortran("coxsuninet",
  #     parm = parm, nobs, nvars, x, ty, tevent, offset, weights, jd, vp, ne, nx, nlam, flmin, ulam, thresh,
  #     maxit, isd, # need to get JHF to reverse these
  #     lmu = integer(1),
  #     ca = double(nx * nlam),
  #     ia = integer(nx),
  #     nin = integer(nlam),
  #     nulldev = double(1),
  #     dev = double(nlam),
  #     alm = double(nlam),
  #     nlp = integer(1),
  #     jerr = integer(1),
  #     iglampos = iglamPos,
  #     olampos = as.double(lamPos),
  #     alpha = alpha, PACKAGE = "sulnet"
  #   )
  # }
  # if (fit$jerr != 0) {
  #   errmsg <- jerr(fit$jerr, maxit, pmax = nx, family = "cox")
  #   if (errmsg$fatal) {
  #     stop(errmsg$msg, call. = FALSE)
  #   } else {
  #     warning(errmsg$msg, call. = FALSE)
  #   }
  # }
  # outlist <- coxgetcoef(fit, nvars, nx, vnames)
  # dev <- fit$dev[seq(fit$lmu)]
  # outlist <- c(outlist, list(dev.ratio = dev, nulldev = fit$nulldev, npasses = fit$nlp, jerr = fit$jerr, offset = is.offset, lamPos = fit$olampos, iglamPos = fit$iglampos, alpha = fit$alpha))
  # class(outlist) <- "coxnet"
  # outlist
}
