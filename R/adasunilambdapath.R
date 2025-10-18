##' @import Matrix

adasunilambdapath  <- function(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax,  missexc, jd, pf,
                        pf2, maxit, lam2, lamPos, loo, negOnly, nobs, nvars, vnames,
                        alpha, ignore_lamPos, asuweight) {
  ################################################################################
  ## data setup
  y <- as.double(y)
  storage.mode(x) <- "double"
  loo <- as.logical(loo)

  ################################################################################
  ## adaptive lasso for choosing coefficients
  if(nvars > nobs & asuweight == "ols") asuweight <- "lasso_ols"
  if(missexc){
    adapf <- switch(asuweight,
                    ols = lm.fit(cbind(rep(1,nobs),x),y)$coefficients[-1],
                    lasso = coef(cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE),
                                 s = fit_ls$lambda.min)[-1],
                    lasso_ols = lasso_ols_supp(x,y),
                    univar = uniFit(x,y)$beta)
    adapf <- 1/pmax.int(abs(adapf), 1e-6)
    adafit <- cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE,
                        pf = adapf, lambda2 = lam2)
  }else{
    adapf <- switch(asuweight,
                    ols = lm.fit(cbind(rep(1,nobs),x[,-jd[-1]]),y)$coefficients[-1],
                    lasso = coef(cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE,
                                           exclude = jd[-1]),
                                 s = fit_ls$lambda.min)[-1],
                    lasso_ols = lasso_ols_supp(x,y,jd = jd),
                    univar = uniFit(x,y)$beta)
    adapf <- 1/pmax.int(abs(adapf), 1e-6)
    adafit <- cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE,
                        pf = adapf, lambda2 = lam2,
                        exclude = jd[-1])
  }

  selected_vars <- which(coef(adafit)[-1] != 0)
  new_jd <- seq(nvars)[-selected_vars]
  new_jd <- append(new_jd, length(new_jd), after = 0)
  new_jd <- as.integer(new_jd)

  ################################################################################
  ## lambda setup
  if(negOnly){
    getlambda <- .Fortran("getlambda", nobs, nvars, nlam, ulam = ulam, x,
                          y, pf, flmin, PACKAGE = "sulnet")
  }else{
    unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
                       beta0 = double(nvars),
                       beta  = double(nvars),
                       fit   = double(nobs * nvars),
                       PACKAGE = "sulnet")
    f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
    storage.mode(f) <- "double"
    getlambda <- .Fortran("getlambda", nobs, nvars, nlam, ulam = ulam, f,
                          y, pf, flmin, PACKAGE = "sulnet")
  }
  ulam <- getlambda$ulam
  flmin = as.double(1)

  ################################################################################
  ## if only computing the negative steps

  lam2 <- as.double(0)

  if(negOnly){
    if(!is.null(alpha)){
      n_alpha = length(alpha)
      b0list = list()
      betamat = vector("list", length = n_alpha)
      dfmat = vector("list", length = n_alpha)
      npassesmat = vector("list", length = n_alpha)
      jerrmat = vector("list", length = n_alpha)

      for(a in seq_along(alpha)){
        fit <- .Fortran("suniwalpha", lam2, lamPos, nobs,nvars, x, as.double(y), new_jd, pf, pf2,
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
        fit <- .Fortran("suniwalpha", lam2,lamPos[a], nobs, nvars, x, as.double(y), new_jd, pf, pf2,
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
      fit <- .Fortran("suniwalpha", lam2, lamPos, nobs, nvars, f, as.double(y), new_jd, pf, pf2,
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
      fit <- .Fortran("suniwalpha", lam2, lamPos[a], nobs, nvars, f, as.double(y), new_jd, pf, pf2,
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
}


