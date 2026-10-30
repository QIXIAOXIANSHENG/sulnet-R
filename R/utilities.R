##' @importFrom stats approx
##' @importFrom methods new
##' @import Matrix
##' @importFrom graphics segments
##' @importFrom stats cor lm.fit sd
#####################################################################
## Some of the functions are minor modifications or direct copies
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
  # for S3 plot method for sulnet2D, suni2D and related.

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

coxgetcoef <- function(fit, nvars, nx, vnames) {
  lmu <- fit$lmu
  if (lmu < 1) {
    ## changed this to a warning message, and return an empty model
    warning("an empty model has been returned; probably a convergence issue")
    coefob <- list(a0 = fit$a0, beta = zeromat(nvars, as.integer(1), vnames, "s0"), df = 0, dim = c(nvars, 1), lambda = Inf)
    return(coefob)
  }
  nin <- fit$nin[seq(lmu)]
  ninmax <- max(nin)
  lam <- fit$alm[seq(lmu)]

  stepnames <- paste("s", seq(lmu) - 1, sep = "")
  dd <- c(nvars, lmu)
  if (ninmax > 0) {
    ca <- matrix(fit$ca[seq(nx * lmu)], nx, lmu)[seq(ninmax), , drop = FALSE]
    df <- apply(abs(ca) > 0, 2, sum)
    ja <- fit$ia[seq(ninmax)]
    #### confusing but too hard to change
    ### glmnet builds a list of ever active variables which is nondecreasing
    ### Since ca was initialized to zero, no harm is done in passing a square matrix
    ### to new(); then when we do a drop0 it makes it really sparse
    oja <- order(ja)
    ja <- rep(ja[oja], lmu)
    ia <- cumsum(c(1, rep(ninmax, lmu)))
    beta <- drop0(new("dgCMatrix", Dim = dd, Dimnames = list(vnames, stepnames), x = as.vector(ca[oja, ]), p = as.integer(ia - 1), i = as.integer(ja - 1)))
  } else {
    beta <- zeromat(nvars, lmu, vnames, stepnames)
    df <- rep(0, lmu)
  }
  a0 <- fit$a0
  if (!is.null(a0)) { # for Cox model
    a0 <- a0[seq(lmu)]
    names(a0) <- stepnames
  }
  list(a0 = a0, beta = beta, df = df, dim = dd, lambda = lam)
}


err <- function(n, maxit, pmax) {
  if (n == 0) {
    msg <- ""
  }
  if (n > 0) {
    if (n < 7777) {
      msg <- "Memory allocation error"
    }
    if (n == 7777) {
      msg <- "All used predictors have zero variance"
    }
    if (n == 10000) {
      msg <- "All penalty factors are <= 0"
    }
    n <- 1
    msg <- paste("in sulnet fortran code -", msg)
  }
  if (n < 0) {
    if (n > -10000) {
      msg <- paste("Convergence for ", -n,
        "th lambda value not reached after maxit=", maxit,
        " iterations; solutions for larger lambdas returned",
        sep = ""
      )
    }
    if (n < -10000) {
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=",
        pmax, " at ", -n - 10000,
        "th lambda value; solutions for larger lambdas returned",
        sep = ""
      )
    }
    n <- -1
    msg <- paste("from sulnet fortran code -", msg)
  }
  list(n = n, msg = msg)
}

## Stop on fatal native errors and warn on incomplete two-dimensional grids.
## Cells exceeding pmax are skipped and the traversal continues; running out
## of passes ends the traversal. jerr reports the first such grid cell in
## traversal order.
grid_err <- function(jerr, maxit, pmax) {
  errmsg <- err(jerr, maxit, pmax)
  if (errmsg$n == 1L) {
    stop(errmsg$msg, call. = FALSE)
  }
  if (jerr < -10000L) {
    warning("from sulnet fortran code - Number of nonzero coefficients ",
      "exceeds pmax=", pmax, " from penalty combination ", -jerr - 10000L,
      " of the traversal on; combinations exceeding pmax are not computed",
      call. = FALSE
    )
  } else if (jerr < 0L) {
    warning("from sulnet fortran code - Convergence for penalty ",
      "combination ", -jerr, " of the traversal not reached within the ",
      "total pass limit ", maxit, "; it and later combinations are not ",
      "computed",
      call. = FALSE
    )
  }
  invisible(NULL)
}

## Act on the per-predictor codes from the univariate logistic fits in
## betalogFIT (see unilogFIT for their meaning).
betalogcheck <- function(jerr, vnames) {
  if (any(jerr == 8002L)) {
    stop("y has a single class; the univariate logistic fits do not exist",
      call. = FALSE
    )
  }
  if (any(jerr == 8001L)) {
    stop("y must be coded -1 or 1 for the univariate logistic fits",
      call. = FALSE
    )
  }
  if (any(jerr == -1L)) {
    warning("univariate logistic fit did not converge for ",
      paste(vnames[jerr == -1L], collapse = ", "),
      " (possibly separation); the last iterate is used",
      call. = FALSE
    )
  }
  invisible(NULL)
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

getmin2D <- function(pen1, pen2, cvm, cvsd, method) {
  # pen1 <- rev(sort(pen1))
  npen2 <- sort(pen2)
  if (!all(pen2 == npen2)) {
    pen2 <- npen2
    cvm <- cvm[rev(seq_len(nrow(cvm))), , drop = FALSE]
    cvsd <- cvsd[rev(seq_len(nrow(cvsd))), , drop = FALSE]
  }
  cvmin <- min(cvm, na.rm = TRUE)
  idmin <- which(cvm <= cvmin, arr.ind = TRUE)

  # first pick pen2/lamPos such that it reaches maximum (preserve most sign consistency)
  pen2.id <- which(max(pen2[idmin[, 1]], na.rm = TRUE) == pen2[idmin[, 1]])
  pen1.id <- which.max(pen1[idmin[pen2.id, 2]])
  pen1.min <- max(pen1[idmin[pen2.id, 2]])
  pen2.min <- pen2[idmin[pen2.id[pen1.id], 1]]

  pen1idmin <- match(pen1.min, pen1)
  pen2idmin <- match(pen2.min, pen2)

  # 1se, same rule for choosing pen2/lamPos and pen1
  semin <- (cvm + cvsd)[pen2idmin, pen1idmin]
  idse <- which(cvm <= semin, arr.ind = TRUE)
  idse <- matrix(idse[idse[, 1] >= pen2idmin, ], ncol = 2)
  idse <- matrix(idse[idse[, 2] <= pen1idmin, ], ncol = 2)

  ## method 1, in the selected region choose the row/col index s.t. cvm max
  # rowse <- which((cvm + cvsd)[idse] == max((cvm + cvsd)[idse]),arr.ind = TRUE)
  ## method 2, in the selected region choose the row/col index s.t. pen2/lamPos max, then pen1 max
  rowse <- which(idse[, 1] == max(idse[, 1]))
  rowse <- rowse[which(idse[rowse, 2] == min(idse[rowse, 2]))]

  idse <- matrix(idse[rowse, ], ncol = 2)

  pen2.id <- which(max(pen2[idse[, 1]]) == pen2[idse[, 1]])
  pen1.id <- which.max(pen1[idse[pen2.id, 2]])
  pen1.se <- max(pen1[idse[pen2.id, 2]])
  pen2.se <- pen2[idse[pen2.id[pen1.id], 1]]

  cv.min <- list(pen1.min, pen2.min)
  names(cv.min) <- strsplit(method, "_")[[1]]
  cv.1se <- list(pen1.se, pen2.se)
  names(cv.1se) <- strsplit(method, "_")[[1]]
  return(list(cv.min = cv.min, cv.1se = cv.1se))
}



#' Unpack a native zeta-tau coefficient grid
#'
#' Convert compressed native theta arrays to sparse coefficient matrices and
#'   preserve grid coordinates and the computed traversal path.
#'
#' @details Native arrays are column-major: active coefficients vary fastest,
#'   then zeta, then tau. Predictor indices in itheta are one-based. Only the
#'   first napenalty rows of apenalty identify computed grid cells, except
#'   cells with a negative ntheta, which the solver skipped for exceeding pmax
#'   (they are also dropped from penalty.path); traversal
#'   order need not match the rectangular storage order. Invalid counts, indices
#'   or grid matches stop. Positive native errors stop and negative errors warn
#'   before partial output is returned. Uncomputed cells have NA intercepts,
#'   degrees of freedom and active-set counts, but zero sparse slope
#'   placeholders.
#'
#' @param fit List returned by the corresponding native grid solver, including
#'   grid dimensions and vectors, napenalty, theta0, theta, itheta, ntheta,
#'   apenalty, and jerr.
#' @param maxit Pass limit used to format native convergence messages.
#' @param pmax Leading dimension of the compressed native coefficient arrays.
#' @param nvars Number of predictors in the full design.
#' @param vnames Predictor names for coefficient matrix rows.
#' @return A list containing theta (one nvars by first-grid-length sparse matrix
#'   per second-coordinate value), theta0 and df matrices, dim, the zeta and tau
#'   vectors, ntheta, napenalty, a logical computed matrix, and penalty.path.
#'   Theta coefficients refer to the transformed predictors, before
#'   reconstruction.
#' @noRd
getoutputZT <- function(fit, maxit, pmax, nvars, vnames) {
  nzeta <- as.integer(fit$nzeta[1L])
  ntau <- as.integer(fit$ntau[1L])

  zeta <- as.numeric(fit$uzeta[seq_len(nzeta)])
  tau <- as.numeric(fit$utau[seq_len(ntau)])

  npenalty <- nzeta * ntau
  napenalty <- as.integer(fit$napenalty[1L])

  if (napenalty < 0L || napenalty > npenalty) {
    stop("invalid napenalty returned by suniZTNET", call. = FALSE)
  }

  ## ------------------------------------------------------------
  ## Check Fortran error
  ## ------------------------------------------------------------
  grid_err(fit$jerr, maxit, pmax)

  ## ------------------------------------------------------------
  ## Names
  ## ------------------------------------------------------------
  zetanames <- paste0("z", seq_len(nzeta) - 1L)
  taunames <- paste0("t", seq_len(ntau) - 1L)

  ## ------------------------------------------------------------
  ## ntheta
  ##
  ## Fortran:
  ##   ntheta(nzeta, ntau)
  ##
  ## Note that ntheta is the number of variables in the active set,
  ## not necessarily the number of nonzero theta values.
  ## ------------------------------------------------------------
  ntheta <- matrix(
    as.integer(fit$ntheta[seq_len(npenalty)]),
    nrow = nzeta,
    ncol = ntau,
    dimnames = list(zetanames, taunames)
  )

  ## ------------------------------------------------------------
  ## theta0
  ##
  ## Fortran:
  ##   theta0(1, nzeta, ntau)
  ##
  ## Since the first dimension is 1, its R vectorization is
  ## equivalent to an nzeta x ntau matrix.
  ## ------------------------------------------------------------
  theta0 <- matrix(
    as.numeric(fit$theta0[seq_len(npenalty)]),
    nrow = nzeta,
    ncol = ntau,
    dimnames = list(zetanames, taunames)
  )

  ## ------------------------------------------------------------
  ## Penalty traversal path
  ##
  ## apenalty(l, 1) = zeta
  ## apenalty(l, 2) = tau
  ##
  ## Only rows 1:napenalty are valid.
  ## ------------------------------------------------------------
  apenalty <- matrix(
    as.numeric(fit$apenalty[seq_len(2L * npenalty)]),
    nrow = npenalty,
    ncol = 2L
  )
  colnames(apenalty) <- c("zeta", "tau")

  if (napenalty > 0L) {
    penalty.path <- apenalty[seq_len(napenalty), , drop = FALSE]
    rownames(penalty.path) <- paste0(
      "p",
      seq_len(napenalty) - 1L
    )
  } else {
    penalty.path <- apenalty[FALSE, , drop = FALSE]
  }

  ## ------------------------------------------------------------
  ## Determine which grid points were actually computed.
  ##
  ## Normally all entries are TRUE. This is mainly needed when
  ## suniZTNET returns partial results after a non-fatal error.
  ## ------------------------------------------------------------
  computed <- matrix(
    FALSE,
    nrow = nzeta,
    ncol = ntau,
    dimnames = list(zetanames, taunames)
  )

  if (napenalty > 0L) {
    iz.path <- match(penalty.path[, 1L], zeta)
    it.path <- match(penalty.path[, 2L], tau)

    if (anyNA(iz.path) || anyNA(it.path)) {
      stop(
        "apenalty returned by suniZTNET does not match uzeta/utau",
        call. = FALSE
      )
    }

    ## Cells skipped for exceeding pmax carry a negative ntheta; they stay
    ## uncomputed and are dropped from the reported path.
    cells <- cbind(iz.path, it.path)
    done <- ntheta[cells] >= 0L
    computed[cells[done, , drop = FALSE]] <- TRUE
    penalty.path <- penalty.path[done, , drop = FALSE]
  }

  ## Uncomputed intercepts should not look like valid zero models.
  theta0[!computed] <- NA_real_

  ## ------------------------------------------------------------
  ## df:
  ##   rows    = zeta
  ##   columns = tau
  ##
  ## Important:
  ## ntheta is the active-set size, whereas df is the actual
  ## number of nonzero theta coefficients.
  ## ------------------------------------------------------------
  df <- matrix(
    NA_integer_,
    nrow = nzeta,
    ncol = ntau,
    dimnames = list(zetanames, taunames)
  )

  ## ------------------------------------------------------------
  ## theta:
  ##
  ## theta[[it]] is an nvars x nzeta dgCMatrix for tau[it].
  ##
  ## Rows    = variables
  ## Columns = zeta values
  ##
  ## We intentionally do NOT first convert fit$theta and fit$itheta
  ## into dense 3D R arrays. Direct linear indexing avoids an
  ## unnecessary large memory copy.
  ## ------------------------------------------------------------
  theta <- vector("list", ntau)
  names(theta) <- taunames

  for (it in seq_len(ntau)) {
    i.parts <- vector("list", nzeta)
    j.parts <- vector("list", nzeta)
    x.parts <- vector("list", nzeta)

    for (iz in seq_len(nzeta)) {
      if (!computed[iz, it]) {
        next
      }

      nk <- ntheta[iz, it]

      if (is.na(nk) || nk < 0L || nk > pmax) {
        stop(
          "invalid ntheta returned by suniZTNET",
          call. = FALSE
        )
      }

      if (nk == 0L) {
        df[iz, it] <- 0L
        next
      }

      ## Fortran array:
      ##
      ##   theta(pmax, nzeta, ntau)
      ##   itheta(pmax, nzeta, ntau)
      ##
      ## R/Fortran both use column-major ordering, so the beginning
      ## of cell (iz, it) is:
      ##
      ##   pmax * ((iz - 1) + nzeta * (it - 1))
      ##
      offset <- pmax * (
        (iz - 1L) +
          nzeta * (it - 1L)
      )

      pos <- offset + seq_len(nk)

      rows <- as.integer(fit$itheta[pos])
      vals <- as.numeric(fit$theta[pos])

      if (any(rows < 1L | rows > nvars)) {
        stop(
          "invalid itheta returned by suniZTNET",
          call. = FALSE
        )
      }

      ## ntheta is active-set size. Some variables in that active
      ## set may have returned to exactly zero.
      keep <- abs(vals) > 0

      df[iz, it] <- sum(keep)

      if (any(keep)) {
        i.parts[[iz]] <- rows[keep]
        j.parts[[iz]] <- rep.int(iz, sum(keep))
        x.parts[[iz]] <- vals[keep]
      }
    }

    ii <- as.integer(
      unlist(i.parts, use.names = FALSE)
    )
    jj <- as.integer(
      unlist(j.parts, use.names = FALSE)
    )
    xx <- as.numeric(
      unlist(x.parts, use.names = FALSE)
    )

    theta[[it]] <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = xx,
      dims = c(nvars, nzeta),
      dimnames = list(vnames, zetanames)
    )
  }

  ## ntheta from an uncomputed grid point should also not be
  ## interpreted as a valid active-set size.
  ntheta[!computed] <- NA_integer_

  list(
    theta0 = theta0,
    theta = theta,
    df = df,
    dim = c(nvars, nzeta, ntau),
    zeta = zeta,
    tau = tau,
    ntheta = ntheta,
    napenalty = napenalty,
    computed = computed,
    penalty.path = penalty.path
  )
}


getoutput <- function(fit, maxit, pmax, nvars, vnames) {
  nalam <- fit$nalam
  nbeta <- fit$nbeta[seq(nalam)]
  nbetamax <- max(nbeta)
  lam <- fit$alam[seq(nalam)]
  stepnames <- paste("s", seq(nalam) - 1, sep = "")
  errmsg <- err(fit$jerr, maxit, pmax)
  switch(paste(errmsg$n),
    `1` = stop(errmsg$msg, call. = FALSE),
    `-1` = print(errmsg$msg, call. = FALSE)
  )
  dd <- c(nvars, nalam)
  if (nbetamax > 0) {
    beta <- matrix(
      fit$beta[seq(pmax * nalam)],
      pmax, nalam
    )[seq(nbetamax), , drop = FALSE]
    df <- apply(abs(beta) > 0, 2, sum)
    ja <- fit$ibeta[seq(nbetamax)]
    oja <- order(ja)
    ja <- rep(ja[oja], nalam)
    ibeta <- cumsum(c(1, rep(nbetamax, nalam)))
    beta <- new("dgCMatrix",
      Dim = dd,
      Dimnames = list(vnames, stepnames),
      x = as.vector(beta[oja, ]),
      p = as.integer(ibeta - 1),
      i = as.integer(ja - 1)
    )
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


#' Unpack a native lambda-alpha coefficient grid
#'
#' Convert compressed native theta arrays to sparse coefficient matrices and
#'   preserve grid coordinates and the computed traversal path.
#'
#' @details Native arrays are column-major: active coefficients vary fastest,
#'   then lambda, then alpha. Predictor indices in itheta are one-based. Only
#'   the first napenalty rows of apenalty identify computed grid cells, except
#'   cells with a negative ntheta, which the solver skipped for exceeding pmax
#'   (they are also dropped from penalty.path);
#'   traversal order need not match the rectangular storage order. Invalid
#'   counts, indices or grid matches stop. Positive native errors stop and
#'   negative errors warn before partial output is returned. Uncomputed cells
#'   have NA intercepts, degrees of freedom and active-set counts, but zero
#'   sparse slope placeholders. The traversal metadata appends zeta = lambda *
#'   (1 - alpha) and tau = lambda * (2 * alpha - 1); these are the common-factor
#'   L1 and extra negative-penalty coordinates, not the two native side
#'   thresholds.
#'
#' @param fit List returned by the corresponding native grid solver, including
#'   grid dimensions and vectors, napenalty, theta0, theta, itheta, ntheta,
#'   apenalty, and jerr.
#' @param maxit Pass limit used to format native convergence messages.
#' @param pmax Leading dimension of the compressed native coefficient arrays.
#' @param nvars Number of predictors in the full design.
#' @param vnames Predictor names for coefficient matrix rows.
#' @return A list containing theta (one nvars by first-grid-length sparse matrix
#'   per second-coordinate value), theta0 and df matrices, dim, the lambda and
#'   alpha vectors, ntheta, napenalty, a logical computed matrix, and
#'   penalty.path. Theta coefficients refer to the transformed predictors,
#'   before reconstruction.
#' @noRd
getoutputLA <- function(fit, maxit, pmax, nvars, vnames) {
  nlambda <- as.integer(fit$nlambda[1L])
  nalpha <- as.integer(fit$nalpha[1L])

  lambda <- as.numeric(fit$ulambda[seq_len(nlambda)])
  alpha <- as.numeric(fit$ualpha[seq_len(nalpha)])

  npenalty <- nlambda * nalpha
  napenalty <- as.integer(fit$napenalty[1L])

  if (napenalty < 0L || napenalty > npenalty) {
    stop("invalid napenalty returned by suniZTNET", call. = FALSE)
  }

  ## ------------------------------------------------------------
  ## Check Fortran error
  ## ------------------------------------------------------------
  grid_err(fit$jerr, maxit, pmax)

  ## ------------------------------------------------------------
  ## Names
  ## ------------------------------------------------------------
  lambdanames <- paste0("l", seq_len(nlambda) - 1L)
  alphanames <- paste0("a", seq_len(nalpha) - 1L)

  ## ------------------------------------------------------------
  ## ntheta
  ##
  ## Fortran:
  ##   ntheta(nlambda, nalpha)
  ##
  ## Note that ntheta is the number of variables in the active set,
  ## not necessarily the number of nonzero theta values.
  ## ------------------------------------------------------------
  ntheta <- matrix(
    as.integer(fit$ntheta[seq_len(npenalty)]),
    nrow = nlambda,
    ncol = nalpha,
    dimnames = list(lambdanames, alphanames)
  )

  ## ------------------------------------------------------------
  ## theta0
  ##
  ## Fortran:
  ##   theta0(1, nlambda, nalpha)
  ##
  ## Since the first dimension is 1, its R vectorization is
  ## equivalent to an nlambda x nalpha matrix.
  ## ------------------------------------------------------------
  theta0 <- matrix(
    as.numeric(fit$theta0[seq_len(npenalty)]),
    nrow = nlambda,
    ncol = nalpha,
    dimnames = list(lambdanames, alphanames)
  )

  ## ------------------------------------------------------------
  ## Penalty traversal path
  ##
  ## apenalty(l, 1) = lambda
  ## apenalty(l, 2) = alpha
  ##
  ## Only rows 1:napenalty are valid.
  ## ------------------------------------------------------------
  apenalty <- matrix(
    as.numeric(fit$apenalty[seq_len(2L * npenalty)]),
    nrow = npenalty,
    ncol = 2L
  )
  colnames(apenalty) <- c("lambda", "alpha")

  if (napenalty > 0L) {
    penalty.path <- apenalty[seq_len(napenalty), , drop = FALSE]
    rownames(penalty.path) <- paste0(
      "p",
      seq_len(napenalty) - 1L
    )
  } else {
    penalty.path <- apenalty[FALSE, , drop = FALSE]
  }
  penalty.path <- cbind(
    penalty.path,
    zeta = penalty.path[, "lambda"] *
      (1 - penalty.path[, "alpha"]),
    tau = penalty.path[, "lambda"] *
      (2 * penalty.path[, "alpha"] - 1)
  )

  ## ------------------------------------------------------------
  ## Determine which grid points were actually computed.
  ##
  ## Normally all entries are TRUE. This is mainly needed when
  ## suniZTNET returns partial results after a non-fatal error.
  ## ------------------------------------------------------------
  computed <- matrix(
    FALSE,
    nrow = nlambda,
    ncol = nalpha,
    dimnames = list(lambdanames, alphanames)
  )

  if (napenalty > 0L) {
    iz.path <- match(penalty.path[, 1L], lambda)
    it.path <- match(penalty.path[, 2L], alpha)

    if (anyNA(iz.path) || anyNA(it.path)) {
      stop(
        "apenalty returned by suniZTNET does not match ulambda/ualpha",
        call. = FALSE
      )
    }

    ## Cells skipped for exceeding pmax carry a negative ntheta; they stay
    ## uncomputed and are dropped from the reported path.
    cells <- cbind(iz.path, it.path)
    done <- ntheta[cells] >= 0L
    computed[cells[done, , drop = FALSE]] <- TRUE
    penalty.path <- penalty.path[done, , drop = FALSE]
  }

  ## Uncomputed intercepts should not look like valid zero models.
  theta0[!computed] <- NA_real_

  ## ------------------------------------------------------------
  ## df:
  ##   rows    = lambda
  ##   columns = alpha
  ##
  ## Important:
  ## ntheta is the active-set size, whereas df is the actual
  ## number of nonzero theta coefficients.
  ## ------------------------------------------------------------
  df <- matrix(
    NA_integer_,
    nrow = nlambda,
    ncol = nalpha,
    dimnames = list(lambdanames, alphanames)
  )

  ## ------------------------------------------------------------
  ## theta:
  ##
  ## theta[[it]] is an nvars x nlambda dgCMatrix for alpha[it].
  ##
  ## Rows    = variables
  ## Columns = lambda values
  ##
  ## We intentionally do NOT first convert fit$theta and fit$itheta
  ## into dense 3D R arrays. Direct linear indexing avoids an
  ## unnecessary large memory copy.
  ## ------------------------------------------------------------
  theta <- vector("list", nalpha)
  names(theta) <- alphanames

  for (it in seq_len(nalpha)) {
    i.parts <- vector("list", nlambda)
    j.parts <- vector("list", nlambda)
    x.parts <- vector("list", nlambda)

    for (iz in seq_len(nlambda)) {
      if (!computed[iz, it]) {
        next
      }

      nk <- ntheta[iz, it]

      if (is.na(nk) || nk < 0L || nk > pmax) {
        stop(
          "invalid ntheta returned by suniZTNET",
          call. = FALSE
        )
      }

      if (nk == 0L) {
        df[iz, it] <- 0L
        next
      }

      ## Fortran array:
      ##
      ##   theta(pmax, nlambda, nalpha)
      ##   itheta(pmax, nlambda, nalpha)
      ##
      ## R/Fortran both use column-major ordering, so the beginning
      ## of cell (iz, it) is:
      ##
      ##   pmax * ((iz - 1) + nlambda * (it - 1))
      ##
      offset <- pmax * (
        (iz - 1L) +
          nlambda * (it - 1L)
      )

      pos <- offset + seq_len(nk)

      rows <- as.integer(fit$itheta[pos])
      vals <- as.numeric(fit$theta[pos])

      if (any(rows < 1L | rows > nvars)) {
        stop(
          "invalid itheta returned by suniZTNET",
          call. = FALSE
        )
      }

      ## ntheta is active-set size. Some variables in that active
      ## set may have returned to exactly zero.
      keep <- abs(vals) > 0

      df[iz, it] <- sum(keep)

      if (any(keep)) {
        i.parts[[iz]] <- rows[keep]
        j.parts[[iz]] <- rep.int(iz, sum(keep))
        x.parts[[iz]] <- vals[keep]
      }
    }

    ii <- as.integer(
      unlist(i.parts, use.names = FALSE)
    )
    jj <- as.integer(
      unlist(j.parts, use.names = FALSE)
    )
    xx <- as.numeric(
      unlist(x.parts, use.names = FALSE)
    )

    theta[[it]] <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = xx,
      dims = c(nvars, nlambda),
      dimnames = list(vnames, lambdanames)
    )
  }

  ## ntheta from an uncomputed grid point should also not be
  ## interpreted as a valid active-set size.
  ntheta[!computed] <- NA_integer_

  list(
    theta0 = theta0,
    theta = theta,
    df = df,
    dim = c(nvars, nlambda, nalpha),
    lambda = lambda,
    alpha = alpha,
    ntheta = ntheta,
    napenalty = napenalty,
    computed = computed,
    penalty.path = penalty.path
  )
}


#' Unpack a native sum-rate coefficient grid
#'
#' Convert compressed native theta arrays to sparse coefficient matrices and
#'   preserve grid coordinates and the computed traversal path.
#'
#' @details Native arrays are column-major: active coefficients vary fastest,
#'   then sum, then rate. Predictor indices in itheta are one-based. Only the
#'   first napenalty rows of apenalty identify computed grid cells, except
#'   cells with a negative ntheta, which the solver skipped for exceeding pmax
#'   (they are also dropped from penalty.path); traversal
#'   order need not match the rectangular storage order. Invalid counts, indices
#'   or grid matches stop. Positive native errors stop and negative errors warn
#'   before partial output is returned. Uncomputed cells have NA intercepts,
#'   degrees of freedom and active-set counts, but zero sparse slope
#'   placeholders. The traversal metadata appends zeta = sum/(1 + rate) and tau
#'   = sum * (rate - 1)/(1 + rate); these are the common-factor L1 and extra
#'   negative-penalty coordinates.
#'
#' @param fit List returned by the corresponding native grid solver, including
#'   grid dimensions and vectors, napenalty, theta0, theta, itheta, ntheta,
#'   apenalty, and jerr.
#' @param maxit Pass limit used to format native convergence messages.
#' @param pmax Leading dimension of the compressed native coefficient arrays.
#' @param nvars Number of predictors in the full design.
#' @param vnames Predictor names for coefficient matrix rows.
#' @return A list containing theta (one nvars by first-grid-length sparse matrix
#'   per second-coordinate value), theta0 and df matrices, dim, the sum and rate
#'   vectors, ntheta, napenalty, a logical computed matrix, and penalty.path.
#'   Theta coefficients refer to the transformed predictors, before
#'   reconstruction.
#' @noRd
getoutputSR <- function(fit, maxit, pmax, nvars, vnames) {
  nsum <- as.integer(fit$nsum[1L])
  nrate <- as.integer(fit$nrate[1L])

  sum <- as.numeric(fit$usum[seq_len(nsum)])
  rate <- as.numeric(fit$urate[seq_len(nrate)])

  npenalty <- nsum * nrate
  napenalty <- as.integer(fit$napenalty[1L])

  if (napenalty < 0L || napenalty > npenalty) {
    stop("invalid napenalty returned by suniZTNET", call. = FALSE)
  }

  ## ------------------------------------------------------------
  ## Check Fortran error
  ## ------------------------------------------------------------
  grid_err(fit$jerr, maxit, pmax)

  ## ------------------------------------------------------------
  ## Names
  ## ------------------------------------------------------------
  sumnames <- paste0("l", seq_len(nsum) - 1L)
  ratenames <- paste0("a", seq_len(nrate) - 1L)

  ## ------------------------------------------------------------
  ## ntheta
  ##
  ## Fortran:
  ##   ntheta(nsum, nrate)
  ##
  ## Note that ntheta is the number of variables in the active set,
  ## not necessarily the number of nonzero theta values.
  ## ------------------------------------------------------------
  ntheta <- matrix(
    as.integer(fit$ntheta[seq_len(npenalty)]),
    nrow = nsum,
    ncol = nrate,
    dimnames = list(sumnames, ratenames)
  )

  ## ------------------------------------------------------------
  ## theta0
  ##
  ## Fortran:
  ##   theta0(1, nsum, nrate)
  ##
  ## Since the first dimension is 1, its R vectorization is
  ## equivalent to an nsum x nrate matrix.
  ## ------------------------------------------------------------
  theta0 <- matrix(
    as.numeric(fit$theta0[seq_len(npenalty)]),
    nrow = nsum,
    ncol = nrate,
    dimnames = list(sumnames, ratenames)
  )

  ## ------------------------------------------------------------
  ## Penalty traversal path
  ##
  ## apenalty(l, 1) = sum
  ## apenalty(l, 2) = rate
  ##
  ## Only rows 1:napenalty are valid.
  ## ------------------------------------------------------------
  apenalty <- matrix(
    as.numeric(fit$apenalty[seq_len(2L * npenalty)]),
    nrow = npenalty,
    ncol = 2L
  )
  colnames(apenalty) <- c("sum", "rate")

  if (napenalty > 0L) {
    penalty.path <- apenalty[seq_len(napenalty), , drop = FALSE]
    rownames(penalty.path) <- paste0(
      "p",
      seq_len(napenalty) - 1L
    )
  } else {
    penalty.path <- apenalty[FALSE, , drop = FALSE]
  }
  penalty.path <- cbind(
    penalty.path,
    zeta = penalty.path[, "sum"] /
      (1 + penalty.path[, "rate"]),
    tau = penalty.path[, "sum"] *
      (penalty.path[, "rate"] - 1) / (1 + penalty.path[, "rate"])
  )

  ## ------------------------------------------------------------
  ## Determine which grid points were actually computed.
  ##
  ## Normally all entries are TRUE. This is mainly needed when
  ## suniZTNET returns partial results after a non-fatal error.
  ## ------------------------------------------------------------
  computed <- matrix(
    FALSE,
    nrow = nsum,
    ncol = nrate,
    dimnames = list(sumnames, ratenames)
  )

  if (napenalty > 0L) {
    iz.path <- match(penalty.path[, 1L], sum)
    it.path <- match(penalty.path[, 2L], rate)

    if (anyNA(iz.path) || anyNA(it.path)) {
      stop(
        "apenalty returned by suniZTNET does not match usum/urate",
        call. = FALSE
      )
    }

    ## Cells skipped for exceeding pmax carry a negative ntheta; they stay
    ## uncomputed and are dropped from the reported path.
    cells <- cbind(iz.path, it.path)
    done <- ntheta[cells] >= 0L
    computed[cells[done, , drop = FALSE]] <- TRUE
    penalty.path <- penalty.path[done, , drop = FALSE]
  }

  ## Uncomputed intercepts should not look like valid zero models.
  theta0[!computed] <- NA_real_

  ## ------------------------------------------------------------
  ## df:
  ##   rows    = sum
  ##   columns = rate
  ##
  ## Important:
  ## ntheta is the active-set size, whereas df is the actual
  ## number of nonzero theta coefficients.
  ## ------------------------------------------------------------
  df <- matrix(
    NA_integer_,
    nrow = nsum,
    ncol = nrate,
    dimnames = list(sumnames, ratenames)
  )

  ## ------------------------------------------------------------
  ## theta:
  ##
  ## theta[[it]] is an nvars x nsum dgCMatrix for rate[it].
  ##
  ## Rows    = variables
  ## Columns = sum values
  ##
  ## We intentionally do NOT first convert fit$theta and fit$itheta
  ## into dense 3D R arrays. Direct linear indexing avoids an
  ## unnecessary large memory copy.
  ## ------------------------------------------------------------
  theta <- vector("list", nrate)
  names(theta) <- ratenames

  for (it in seq_len(nrate)) {
    i.parts <- vector("list", nsum)
    j.parts <- vector("list", nsum)
    x.parts <- vector("list", nsum)

    for (iz in seq_len(nsum)) {
      if (!computed[iz, it]) {
        next
      }

      nk <- ntheta[iz, it]

      if (is.na(nk) || nk < 0L || nk > pmax) {
        stop(
          "invalid ntheta returned by suniZTNET",
          call. = FALSE
        )
      }

      if (nk == 0L) {
        df[iz, it] <- 0L
        next
      }

      ## Fortran array:
      ##
      ##   theta(pmax, nsum, nrate)
      ##   itheta(pmax, nsum, nrate)
      ##
      ## R/Fortran both use column-major ordering, so the beginning
      ## of cell (iz, it) is:
      ##
      ##   pmax * ((iz - 1) + nsum * (it - 1))
      ##
      offset <- pmax * (
        (iz - 1L) +
          nsum * (it - 1L)
      )

      pos <- offset + seq_len(nk)

      rows <- as.integer(fit$itheta[pos])
      vals <- as.numeric(fit$theta[pos])

      if (any(rows < 1L | rows > nvars)) {
        stop(
          "invalid itheta returned by suniZTNET",
          call. = FALSE
        )
      }

      ## ntheta is active-set size. Some variables in that active
      ## set may have returned to exactly zero.
      keep <- abs(vals) > 0

      df[iz, it] <- sum(keep)

      if (any(keep)) {
        i.parts[[iz]] <- rows[keep]
        j.parts[[iz]] <- rep.int(iz, sum(keep))
        x.parts[[iz]] <- vals[keep]
      }
    }

    ii <- as.integer(
      unlist(i.parts, use.names = FALSE)
    )
    jj <- as.integer(
      unlist(j.parts, use.names = FALSE)
    )
    xx <- as.numeric(
      unlist(x.parts, use.names = FALSE)
    )

    theta[[it]] <- Matrix::sparseMatrix(
      i = ii,
      j = jj,
      x = xx,
      dims = c(nvars, nsum),
      dimnames = list(vnames, sumnames)
    )
  }

  ## ntheta from an uncomputed grid point should also not be
  ## interpreted as a valid active-set size.
  ntheta[!computed] <- NA_integer_

  list(
    theta0 = theta0,
    theta = theta,
    df = df,
    dim = c(nvars, nsum, nrate),
    sum = sum,
    rate = rate,
    ntheta = ntheta,
    napenalty = napenalty,
    computed = computed,
    penalty.path = penalty.path
  )
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
    sfrac <- (lambda[1] - s) / (lambda[1] - lambda[k])
    lambda <- (lambda[1] - lambda) / (lambda[1] - lambda[k])
    coord <- approx(lambda, seq(lambda), sfrac)$y
    left <- floor(coord)
    right <- ceiling(coord)
    sfrac <- (sfrac - lambda[right]) / (lambda[left] - lambda[right])
    sfrac[left == right] <- 1
  }
  list(left = left, right = right, frac = sfrac)
}

lamextend <- function(lambda, flmin, nlam, ulam) {
  if (flmin >= 1) {
    ulam
  } else {
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

lasso_w_supp <- function(x, y, jd = NULL) {
  if (is.null(jd)) {
    fit_ls <- cv.sulnet(x, y, method = "ls", standardize = TRUE, intercept = TRUE)
    lasso <- coef(fit_ls,
      s = fit_ls$lambda.min
    )[-1]
  } else {
    fit_ls <- cv.sulnet(x, y,
      method = "ls", standardize = TRUE, intercept = TRUE,
      exclude = jd[-1]
    )
    lasso <- coef(fit_ls,
      s = fit_ls$lambda.min
    )[-1]
  }
  pmax(abs(lasso), 1e-6)
}

lasso_ols_supp <- function(x, y, jd = NULL) {
  if (is.null(jd)) {
    fit_ls <- cv.sulnet(x, y)
    beta_ls <- coef(fit_ls, s = fit_ls$lambda.min)[-1]
    x_new <- x[, which(as.vector(beta_ls) != 0)]
    ols <- lm.fit(cbind(1, x_new), y)$coefficients[-1]
    ols_pf <- rep(1e-6, ncol(x))
    ols_pf[which(as.vector(beta_ls) != 0)] <- ols
    ols_pf
  } else {
    beta_ls <- coef(
      cv.sulnet(x, y,
        method = "ls", standardize = TRUE, intercept = TRUE,
        exclude = jd[-1]
      ),
      s = fit_ls$lambda.min
    )[-1]
    x_new <- x[, which(as.vector(beta_ls) != 0)]
    ols <- lm.fit(cbind(1, x_new), y)$coefficients[-1]
    ols_pf <- rep(1e-6, ncol(x))
    ols_pf[which(as.vector(beta_ls) != 0)] <- ols
    ols_pf
  }
}

multuni <- function(beta_temp, beta0_temp, unibeta, unibeta0) {
  row_idx <- beta_temp@i + 1
  col_ptrs <- beta_temp@p
  col_idx <- rep(seq_along(col_ptrs[-1]), diff(col_ptrs))
  beta_result <- beta_temp
  beta_result@x <- beta_temp@x * unibeta[cbind(row_idx, col_idx)]
  list(
    beta = beta_result,
    b0 = beta0_temp + colSums(unibeta0 * beta_temp)
  )
}

nonzero <- function(beta, bystep = FALSE) {
  ns <- ncol(beta)
  ## beta should be in 'dgCMatrix' format
  if (nrow(beta) == 1) {
    if (bystep) {
      apply(
        beta, 2,
        function(x) if (abs(x) > 0) 1 else NULL
      )
    } else {
      if (any(abs(beta) > 0)) {
        1
      } else {
        NULL
      }
    }
  } else {
    beta <- t(beta)
    which <- diff(beta@p)
    which <- seq(which)[which > 0]
    if (bystep) {
      nzel <- function(x, which) {
        if (any(x)) {
          which[x]
        } else {
          NULL
        }
      }
      beta <- abs(as.matrix(beta[, which])) > 0
      if (ns == 1) {
        apply(beta, 2, nzel, which)
      } else {
        apply(beta, 1, nzel, which)
      }
    } else {
      which
    }
  }
}

#' Reconstruct coefficients on the original predictor scale
#'
#' Combine transformed-design coefficients with univariate slopes and offsets.
#'
#' @details Uses \eqn{\gamma_j = \beta_j\theta_j} and \eqn{\gamma_0 = \theta_0 +
#'   \sum_j \beta_{0j}\theta_j}. The intercept adjustment uses theta, not gamma.
#'   Nonzero counts are recomputed from gamma; consequently zero placeholders at
#'   uncomputed grid cells are also counted as zero.
#' @param out Converted grid output with theta as a nonempty list of equally
#'   sized sparse matrices and theta0 as a first-grid by second-grid matrix.
#' @param betaFIT First-stage native result containing length-p beta and beta0.
#' @return The input metadata with theta and theta0 removed, gamma and gamma0
#'   added in the same grid layout, and df replaced by counts of nonzero gamma.
#'   Inconsistent coefficient dimensions produce an error.
#' @noRd
#' @importFrom methods as
theta2gamma <- function(out, betaFIT) {
  if (is.null(out$theta) || is.null(out$theta0)) {
    stop("out must contain both 'theta' and 'theta0'.", call. = FALSE)
  }

  if (!is.list(out$theta) || length(out$theta) == 0L) {
    stop("out$theta must be a non-empty list of sparse matrices.",
      call. = FALSE
    )
  }

  nvars <- nrow(out$theta[[1L]])
  ngrid1 <- ncol(out$theta[[1L]])
  ngrid2 <- length(out$theta)

  beta <- as.numeric(betaFIT$beta[seq_len(nvars)])
  beta0 <- as.numeric(betaFIT$beta0[seq_len(nvars)])

  if (length(beta) != nvars) {
    stop("length of betaFIT$beta does not match nvars.", call. = FALSE)
  }

  if (length(beta0) != nvars) {
    stop("length of betaFIT$beta0 does not match nvars.", call. = FALSE)
  }

  if (!all(dim(out$theta0) == c(ngrid1, ngrid2))) {
    stop(
      "dimensions of out$theta0 are inconsistent with out$theta.",
      call. = FALSE
    )
  }

  gamma <- vector("list", ngrid2)
  gamma0 <- out$theta0

  if (!is.null(names(out$theta))) {
    names(gamma) <- names(out$theta)
  }

  for (k in seq_len(ngrid2)) {
    theta_k <- out$theta[[k]]

    if (!inherits(theta_k, "dgCMatrix")) {
      theta_k <- as(theta_k, "dgCMatrix")
    }

    if (!all(dim(theta_k) == c(nvars, ngrid1))) {
      stop(
        "all elements of out$theta must have the same dimensions.",
        call. = FALSE
      )
    }

    ## ----------------------------------------------------------
    ## gamma_j = theta_j * beta_j
    ##
    ## dgCMatrix stores row index in @i (0-based).
    ## Therefore each stored nonzero theta can be multiplied
    ## directly by the corresponding beta without creating
    ## any dense matrices.
    ## ----------------------------------------------------------
    gamma_k <- theta_k

    if (length(theta_k@x) > 0L) {
      row_idx <- theta_k@i + 1L
      gamma_k@x <- theta_k@x * beta[row_idx]
    }

    gamma[[k]] <- Matrix::drop0(gamma_k)

    ## ----------------------------------------------------------
    ## gamma0 = theta0 + sum_j theta_j * beta0_j
    ##
    ## This uses theta, not gamma.
    ##
    ## crossprod(beta0, theta_k) returns one value for each
    ## column / first penalty dimension.
    ## ----------------------------------------------------------
    intercept_adjustment <- as.numeric(
      Matrix::crossprod(beta0, theta_k)
    )

    gamma0[, k] <- out$theta0[, k] + intercept_adjustment
  }

  ## ------------------------------------------------------------
  ## Number of nonzero gamma coefficients for each model
  ## ------------------------------------------------------------
  df <- matrix(
    0L,
    nrow = ngrid1,
    ncol = ngrid2,
    dimnames = dimnames(out$theta0)
  )

  for (k in seq_len(ngrid2)) {
    df[, k] <- Matrix::colSums(gamma[[k]] != 0)
  }

  ## ------------------------------------------------------------
  ## Preserve the metadata in out, but replace theta objects
  ## by final gamma objects.
  ## ------------------------------------------------------------
  result <- out

  result$gamma0 <- gamma0
  result$gamma <- gamma
  result$df <- df

  ## theta/theta0 are intermediate coefficients and are removed
  ## from the final returned object.
  result$theta0 <- NULL
  result$theta <- NULL

  result
}


zeromat <- function(nvars, nalam, vnames, stepnames) {
  ca <- rep(0, nalam)
  ia <- seq(nalam + 1)
  ja <- rep(1, nalam)
  dd <- c(nvars, nalam)
  new("dgCMatrix",
    Dim = dd, Dimnames = list(vnames, stepnames),
    x = as.vector(ca), p = as.integer(ia - 1), i = as.integer(ja - 1)
  )
}
