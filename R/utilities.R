##' @importFrom stats approx
##' @importFrom methods new
##' @import Matrix
##' @importFrom graphics segments
#####################################################################
## These functions are minor modifications or direct copies
## from the `glmnet` package:
##
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate Descent.
##   Journal of Statistical Software, 33(1), 1-22.
##   URL: https://www.jstatsoft.org/v33/i01/.
##
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.

canvas_structure <- function(number_of_plots) {
  # original function, intended for giving canvas distribution
  # for S3 plot method for sulnet2D and related.

  sqrt_result <- sqrt(number_of_plots)
  col <- floor(sqrt_result)
  diff <- sqrt_result - col

  if (diff > 0 && diff < 0.5) {
    row <- ceiling(sqrt_result)
  } else if (diff >= 0.5) {
    row <- ceiling(sqrt_result) + 1
  } else {
    row <- col
  }

  return(c(row, col))
}

coxgetcoef=function(fit,nvars,nx,vnames){
  lmu=fit$lmu
  if(lmu<1){
    ## changed this to a warning message, and return an empty model
    warning("an empty model has been returned; probably a convergence issue")
    coefob=list(a0=fit$a0,beta=zeromat(nvars,as.integer(1),vnames,"s0"),df=0,dim=c(nvars,1),lambda=Inf)
    return(coefob)
  }
  nin=fit$nin[seq(lmu)]
  ninmax=max(nin)
  lam=fit$alm[seq(lmu)]

  stepnames=paste("s",seq(lmu)-1,sep="")
  dd=c(nvars,lmu)
  if(ninmax>0){
           ca=matrix(fit$ca[seq(nx*lmu)],nx,lmu)[seq(ninmax),,drop=FALSE]
           df=apply(abs(ca)>0,2,sum)
           ja=fit$ia[seq(ninmax)]
####confusing but too hard to change
###glmnet builds a list of ever active variables which is nondecreasing
###Since ca was initialized to zero, no harm is done in passing a square matrix
###to new(); then when we do a drop0 it makes it really sparse
           oja=order(ja)
           ja=rep(ja[oja],lmu)
           ia=cumsum(c(1,rep(ninmax,lmu)))
           beta=drop0(new("dgCMatrix",Dim=dd,Dimnames=list(vnames,stepnames),x=as.vector(ca[oja,]),p=as.integer(ia-1),i=as.integer(ja-1)))
         }else {
           beta = zeromat(nvars,lmu,vnames,stepnames)
           df=rep(0,lmu)
         }
  a0=fit$a0
  if(!is.null(a0)){#for Cox model
    a0=a0[seq(lmu)]
  names(a0)=stepnames
  }
  list(a0=a0,beta=beta,df=df,dim=dd,lambda=lam)
}


err <- function(n, maxit, pmax) {
  if (n == 0)
    msg <- ""
  if (n > 0) {
    if (n < 7777)
      msg <- "Memory allocation error"
    if (n == 7777)
      msg <- "All used predictors have zero variance"
    if (n == 10000)
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in sulnet fortran code -", msg)
  }
  if (n < 0) {
    if (n > -10000)
      msg <- paste("Convergence for ", -n,
                   "th lambda value not reached after maxit=", maxit,
                   " iterations; solutions for larger lambdas returned",
                   sep = "")
    if (n < -10000)
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
                   pmax, " at ", -n - 10000,
                   "th lambda value; solutions for larger lambdas returned",
                   sep = "")
    n <- -1
    msg <- paste("from sulnet fortran code -", msg)
  }
  list(n = n, msg = msg)
}



error.bars <- function(x, upper, lower, width = 0.02,
                       ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

getmin <- function(lambda, cvm, cvsd) {
  cvmin <- min(cvm)
  idmin <- cvm <= cvmin
  lambda.min <- max(lambda[idmin])
  idmin <- match(lambda.min, lambda)
  semin <- (cvm + cvsd)[idmin]
  idmin <- cvm <= semin
  ## cat('\n\nidmin\n\n',idmin)
  ## cat('\n\nlambda[idmin]\n\n',lambda[idmin])
  ## cat('\n\nmax\n\n',max(lambda[idmin]))
  lambda.1se <- max(lambda[idmin])
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}

getmin2D <- function(lambda, alpha, cvm, cvsd, use_alpha) {
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- which(cvm <= cvmin, arr.ind = TRUE)

  # first pick alpha/lamPos such that it reaches maximum (preserve most sign consistency)
  alpha.id <- which(max(alpha[idmin[,1]], na.rm = TRUE) == alpha[idmin[,1]])
  lambda.id <- which.max(lambda[idmin[alpha.id,2]])
  lambda.min <- max(lambda[idmin[alpha.id,2]])
  alpha.min <- alpha[idmin[alpha.id[lambda.id],1]]

  lambdaidmin <- match(lambda.min, lambda)
  alphaidmin <- match(alpha.min, alpha)

  # 1se, same rule for choosing alpha/lamPos and lambda
  semin <- (cvm + cvsd)[alphaidmin,lambdaidmin]
  idse <- which(cvm <= semin, arr.ind = TRUE)
  idse <- matrix(idse[idse[,1]>=alphaidmin,],ncol = 2)
  idse <- matrix(idse[idse[,2]<=lambdaidmin,],ncol = 2)

  ## method 1, in the selected region choose the row/col index s.t. cvm max
  #rowse <- which((cvm + cvsd)[idse] == max((cvm + cvsd)[idse]),arr.ind = TRUE)
  ## method 2, in the selected region choose the row/col index s.t. alpha/lamPos max, then lambda max
  rowse <- which(idse[,1] == max(idse[,1]))
  rowse <- rowse[which(idse[rowse,2] == min(idse[rowse,2]))]

  idse <- matrix(idse[rowse,],ncol = 2)

  alpha.id <- which(max(alpha[idse[,1]]) == alpha[idse[,1]])
  lambda.id <- which.max(lambda[idse[alpha.id,2]])
  lambda.se <- max(lambda[idse[alpha.id,2]])
  alpha.se <- alpha[idse[alpha.id[lambda.id],1]]

  if(use_alpha){
    return(list(cv.min = list(lambda.min = lambda.min,
                              alpha.min = alpha.min),
                cv.1se = list(lambda.1se = lambda.se,
                              alpha.1se = alpha.se)))
  }else{
    return(list(cv.min = list(lambda.min = lambda.min,
                              lamPos.min = alpha.min),
                cv.1se = list(lambda.1se = lambda.se,
                              lamPos.1se = alpha.se)))
  }

}



getoutput <- function(fit, maxit, pmax, nvars, vnames) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- err(fit$jerr, maxit, pmax)
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE),
         `-1` = print(errmsg$msg, call. = FALSE))
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(fit$beta[seq(pmax * nalam)],
                   pmax, nalam)[seq(nbetamax), , drop = FALSE]
    df <- apply(abs(beta) > 0, 2, sum)
    ja <- fit$ibeta[seq(nbetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nalam)
    ibeta <- cumsum(c(1, rep(nbetamax, nalam)))
    beta <- new("dgCMatrix", Dim = dd,
                Dimnames = list(vnames, stepnames),
                x = as.vector(beta[oja, ]),
                p = as.integer(ibeta - 1),
                i = as.integer(ja - 1))
  } else {
    beta <- zeromat(nvars, nalam, vnames, stepnames)
    df <- rep(0, nalam)
  }
  b0 <- fit$b0
  if (!is.null(b0)) {
    b0 <- b0[seq(nalam)]
    names(b0) <- stepnames
  }
  list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
}



lambda.interp <- function(lambda, s) {
  ## lambda is the index sequence that is produced by the model
  ## s is the new vector at which evaluations are required.
  ## the value is a vector of left and right indices, and a vector of fractions.
  ## the new values are interpolated bewteen the two using the fraction
  ## Note: lambda decreases. you take: sfrac*left+(1-sfrac*right)
  if (length(lambda) == 1) {
    nums <- length(s)
    left <- rep(1, nums)
    right <- left
    sfrac <- rep(1, nums)
  } else {
    s[s > max(lambda)] <- max(lambda)
    s[s < min(lambda)] <- min(lambda)
    k <- length(lambda)
    sfrac <- (lambda[1] - s)/(lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda)/(lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right])/(lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

lamextend <- function(lambda, flmin, nlam, ulam){
  if(flmin >= 1){
    ulam
  }else{
    lamstart <- log(lambda[1])
    lamend <- lamstart + log(flmin)
    lambda <- exp(seq(lamstart, lamend, length.out = nlam))
    lambda
  }
}

lamfix <- function(lam) {
  llam <- log(lam)
  lam[1] <- exp(2 * llam[2] - llam[3])
  lam
}

lasso_w_supp <- function(x,y,jd = NULL){
  if(is.null(jd)){
    fit_ls <- cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE)
    lasso <- coef(fit_ls,
                 s = fit_ls$lambda.min)[-1]
  }else{
    fit_ls <- cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE,
                        exclude = jd[-1])
    lasso <- coef(fit_ls,
                 s = fit_ls$lambda.min)[-1]
  }
  pmax(abs(lasso),1e-6)
}

lasso_ols_supp <- function(x,y,jd = NULL){
  if(is.null(jd)){
    fit_ls = cv.sulnet(x,y)
    beta_ls = coef(fit_ls, s = fit_ls$lambda.min)[-1]
    x_new = x[, which(as.vector(beta_ls) != 0)]
    ols = lm.fit(cbind(1, x_new), y)$coefficients[-1]
    ols_pf = rep(1e-6,ncol(x))
    ols_pf[which(as.vector(beta_ls) != 0)] = ols
    ols_pf
  }else{
    beta_ls = coef(cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE,
                           exclude = jd[-1]),
                 s = fit_ls$lambda.min)[-1]
    x_new = x[, which(as.vector(beta_ls) != 0)]
    ols = lm.fit(cbind(1, x_new), y)$coefficients[-1]
    ols_pf = rep(1e-6,ncol(x))
    ols_pf[which(as.vector(beta_ls) != 0)] = ols
    ols_pf
  }
}

multuni <- function(beta_temp, beta0_temp, unibeta, unibeta0){
  row_idx <- beta_temp@i + 1
  col_ptrs <- beta_temp@p
  col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
  beta_result <- beta_temp
  beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]
  list(beta = beta_result,
       b0 = beta0_temp + colSums(unibeta0 * beta_temp)
  )
}

nonzero <- function(beta, bystep = FALSE) {
  ns <- ncol(beta)
  ##beta should be in 'dgCMatrix' format
  if (nrow(beta) == 1) {
    if (bystep)
      apply(beta, 2,
            function(x) if (abs(x) > 0) 1 else NULL) else {
                                                       if (any(abs(beta) > 0))
                                                         1 else NULL
                                                     }
  } else {
    beta <- t(beta)
    which <- diff(beta@p)
    which <- seq(which)[which > 0]
    if (bystep) {
      nzel <- function(x, which) if (any(x))
                                   which[x] else NULL
      beta <- abs(as.matrix(beta[, which])) > 0
      if (ns == 1)
        apply(beta, 2, nzel, which) else apply(beta, 1, nzel, which)
    } else which
  }
}



zeromat <- function(nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames),
      x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 1))
}
