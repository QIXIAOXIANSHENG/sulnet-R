##' @import Matrix

sunipath_2 <- function(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, pf,
                     pf2, maxit, lam2, lamPos, loo, negOnly, nobs, nvars, vnames,
                     alpha, ignore_lamPos) {
  ################################################################################
  ## data setup
  y <- as.double(y)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)


  ################################################################################
  ## if only computing the negative steps

  if(negOnly){
    if(!is.null(alpha)){
      n_alpha = length(alpha)
      b0list = list()
      betamat = vector("list", length = n_alpha)
      dfmat = vector("list", length = n_alpha)
      npassesmat = vector("list", length = n_alpha)
      jerrmat = vector("list", length = n_alpha)

      for(a in seq_along(alpha)){
        fit <- .Fortran("suniwalpha", lam2, lamPos, nobs, nvars, x, as.double(y), jd, pf, pf2,
                        dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                        nalam = integer(1), b0 = double(nlam),
                        beta = double(pmax * nlam), ibeta = integer(pmax),
                        nbeta = integer(nlam), alam = double(nlam),
                        npass = integer(1), jerr = integer(1),
                        alpha = as.double(alpha[a]), iglamPos = as.logical(ignore_lamPos),
                        PACKAGE = "sulnet")
        outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
        b0list[[paste0("alpha_", alpha[a])]] = outlist$b0
        betamat[[a]] <- outlist$beta
        dfmat[[a]] <- outlist$df
        npassesmat[[a]] <- fit$npass
        jerrmat[[a]] <- fit$jerr
      }
      alphaname = paste0("alpha_", alpha)
      names(betamat) <- alphaname
      names(dfmat) <- alphaname
      names(npassesmat) <- alphaname
      names(jerrmat) <- alphaname
      dim = c(n_alpha, outlist$dim)
      lambda_total = outlist$lambda

      outlist <- list(b0 = b0list,
                      beta = betamat,
                      df = dfmat,
                      dim = dim,
                      lambda = lambda_total,
                      npasses = npassesmat,
                      jerr = jerrmat,
                      alpha = alpha,
                      use_alpha = ignore_lamPos,
                      negOnly = negOnly)

    }else{
      n_alpha = length(lamPos)
      b0list = list()
      betamat = vector("list", length = n_alpha)
      dfmat = vector("list", length = n_alpha)
      npassesmat = vector("list", length = n_alpha)
      jerrmat = vector("list", length = n_alpha)

      for(a in seq_along(lamPos)){
        fit <- .Fortran("suniwalpha", lam2,lamPos[a], nobs, nvars, x, as.double(y), jd, pf, pf2,
                        dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                        nalam = integer(1), b0 = double(nlam),
                        beta = double(pmax * nlam), ibeta = integer(pmax),
                        nbeta = integer(nlam), alam = double(nlam),
                        npass = integer(1), jerr = integer(1),
                        alpha = as.double(0), iglamPos = as.logical(ignore_lamPos),
                        PACKAGE = "sulnet")
        outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
        b0list[[paste0("lamPos_", lamPos[a])]] = outlist$b0
        betamat[[a]] <- outlist$beta
        dfmat[[a]] <- outlist$df
        npassesmat[[a]] <- fit$npass
        jerrmat[[a]] <- fit$jerr
      }
      alphaname = paste0("lamPos_", lamPos)
      names(betamat) <- alphaname
      names(dfmat) <- alphaname
      names(npassesmat) <- alphaname
      names(jerrmat) <- alphaname
      dim = c(n_alpha, outlist$dim)
      lambda_total = outlist$lambda

      outlist <- list(b0 = b0list,
                      beta = betamat,
                      df = dfmat,
                      dim = dim,
                      lambda = lambda_total,
                      npasses = npassesmat,
                      jerr = jerrmat,
                      lamPos = lamPos,
                      use_alpha = ignore_lamPos,
                      negOnly = negOnly)
    }


    class(outlist) <- c("sunipath_2")
    return(outlist)

  }

  ################################################################################
  ## univariate fit

  unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
                     beta0 = double(nvars),
                     beta  = double(nvars),
                     fit   = double(nobs * nvars),
                     PACKAGE = "sulnet")
  f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)




  if(!is.null(alpha)){
    n_alpha = length(alpha)
    fb0list = list()
    fbetamat = vector("list", length = n_alpha)
    dfmat = vector("list", length = n_alpha)
    npassesmat = vector("list", length = n_alpha)
    jerrmat = vector("list", length = n_alpha)

    b0list = list()
    betamat = vector("list", length = n_alpha)

    for(a in seq_along(alpha)){
      fit <- .Fortran("suniwalpha", lam2, lamPos, nobs, nvars, f, as.double(y), jd, pf, pf2,
                      dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                      nalam = integer(1), b0 = double(nlam),
                      beta = double(pmax * nlam), ibeta = integer(pmax),
                      nbeta = integer(nlam), alam = double(nlam),
                      npass = integer(1), jerr = integer(1),
                      alpha = as.double(alpha[a]), iglamPos = as.logical(ignore_lamPos),
                      PACKAGE = "sulnet")
      outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
      ones = rep(1,fit$nalam)
      unibeta <- outer(unifit$beta, ones)
      unibeta0 <- outer(unifit$beta0, ones)

      beta_temp = outlist$beta
      beta0_temp = outlist$b0

      fb0list[[paste0("alpha_", alpha[a])]] = beta0_temp
      fbetamat[[a]] <- beta_temp
      dfmat[[a]] <- outlist$df
      npassesmat[[a]] <- fit$npass
      jerrmat[[a]] <- fit$jerr


      row_idx <- beta_temp@i + 1
      col_ptrs <- beta_temp@p
      col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
      beta_result <- beta_temp
      beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]

      betamat[[a]] = beta_result
      b0list[[paste0("alpha_", alpha[a])]] = beta0_temp + colSums(unibeta0 * beta_temp)
    }
    alphaname = paste0("alpha_", alpha)
    names(fbetamat) <- alphaname
    names(betamat) <- alphaname
    names(dfmat) <- alphaname
    names(npassesmat) <- alphaname
    names(jerrmat) <- alphaname
    dim = c(n_alpha, outlist$dim)
    lambda_total = outlist$lambda

    outlist <- list(b0 = b0list,
                    beta = betamat,
                    df = dfmat,
                    dim = dim,
                    lambda = lambda_total,
                    npasses = npassesmat,
                    jerr = jerrmat,
                    alpha = alpha,
                    use_alpha = ignore_lamPos,
                    negOnly = negOnly,
                    LOO = loo,
                    univariate.fit = list(beta = unifit$beta,
                                          beta0 = unifit$beta0,
                                          fitted.values = f),
                    fbeta = fbetamat,
                    fb0 = fb0list
                    )

  }else{
    n_alpha = length(lamPos)
    fb0list = list()
    fbetamat = vector("list", length = n_alpha)
    dfmat = vector("list", length = n_alpha)
    npassesmat = vector("list", length = n_alpha)
    jerrmat = vector("list", length = n_alpha)

    b0list = list()
    betamat = vector("list", length = n_alpha)

    for(a in seq_along(lamPos)){
      fit <- .Fortran("suniwalpha", lam2, lamPos[a], nobs, nvars, f, as.double(y), jd, pf, pf2,
                      dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                      nalam = integer(1), b0 = double(nlam),
                      beta = double(pmax * nlam), ibeta = integer(pmax),
                      nbeta = integer(nlam), alam = double(nlam),
                      npass = integer(1), jerr = integer(1),
                      alpha = as.double(0), iglamPos = as.logical(ignore_lamPos),
                      PACKAGE = "sulnet")
      outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
      ones = rep(1,fit$nalam)
      unibeta <- outer(unifit$beta, ones)
      unibeta0 <- outer(unifit$beta0, ones)

      beta_temp = outlist$beta
      beta0_temp = outlist$b0

      fb0list[[paste0("lamPos_", lamPos[a])]] = beta0_temp
      fbetamat[[a]] <- beta_temp
      dfmat[[a]] <- outlist$df
      npassesmat[[a]] <- fit$npass
      jerrmat[[a]] <- fit$jerr


      row_idx <- beta_temp@i + 1
      col_ptrs <- beta_temp@p
      col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
      beta_result <- beta_temp
      beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]

      betamat[[a]] = beta_result
      b0list[[paste0("lamPos_", lamPos[a])]] = beta0_temp + colSums(unibeta0 * beta_temp)
    }
    alphaname = paste0("lamPos_", lamPos)
    names(fbetamat) <- alphaname
    names(betamat) <- alphaname
    names(dfmat) <- alphaname
    names(npassesmat) <- alphaname
    names(jerrmat) <- alphaname
    dim = c(n_alpha, outlist$dim)
    lambda_total = outlist$lambda

    outlist <- list(b0 = b0list,
                    beta = betamat,
                    df = dfmat,
                    dim = dim,
                    lambda = lambda_total,
                    npasses = npassesmat,
                    jerr = jerrmat,
                    lamPos = lamPos,
                    use_alpha = ignore_lamPos,
                    negOnly = negOnly,
                    LOO = loo,
                    univariate.fit = list(beta = unifit$beta,
                                          beta0 = unifit$beta0,
                                          fitted.values = f),
                    fbeta = fbetamat,
                    fb0 = fb0list
    )
  }


  class(outlist) <- c("sunipath_2")
  return(outlist)


  # ################################################################################
  # ## call Fortran core
  # fit <- .Fortran("soft_unilassoNET", lam2,lamPos, nobs, nvars, f, as.double(y), jd, pf, pf2,
  #                 dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
  #                 nalam = integer(1), b0 = double(nlam),
  #                 beta = double(pmax * nlam), ibeta = integer(pmax),
  #                 nbeta = integer(nlam), alam = double(nlam),
  #                 npass = integer(1), jerr = integer(1), PACKAGE = "sulnet")
  # ################################################################################
  # ## output
  # outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  # outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr,
  #                            lamPos = lamPos, LOO = loo,
  #                            univariate.fit = list(beta = unifit$beta,
  #                                                  beta0 = unifit$beta0,
  #                                                  fitted.values = f)
  # ))
  # ones = rep(1,fit$nalam)
  # unibeta <- outer(unifit$beta, ones)
  # unibeta0 <- outer(unifit$beta0, ones)
  # beta_temp = outlist$beta
  # beta0_temp = outlist$b0
  #
  # row_idx <- beta_temp@i + 1
  # col_ptrs <- beta_temp@p
  # col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
  # beta_result <- beta_temp
  # beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]
  #
  # outlist$beta = beta_result
  # outlist$b0 = beta0_temp + colSums(unibeta0 * beta_temp)
  #
  # outlist$fbeta = beta_temp; outlist$fb0 = beta0_temp
  # class(outlist) <- c("sunipath_2")
  # outlist
}


