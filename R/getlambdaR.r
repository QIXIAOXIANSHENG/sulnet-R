getlambdaR <- function(
  nobs, nvars, nlam, ulam = ulam,
  x, y, pf, flmin, loo, family, negOnly, iglamPos
) {
  if (family == "gaussian") {
    if (negOnly) {
      while (TRUE) {
        getlambdagauss <- .Fortran("getlambdagauss", nobs, nvars, nlam,
          ulam = ulam, x,
          y, pf, flmin, PACKAGE = "sulnet"
        )
        ulamtemp <- as.double(getlambdagauss$ulam)
        if (!anyNA(ulamtemp)) {
          flmin <- as.double(1)
          if (iglamPos) {
            ulam <- as.double(2 * ulamtemp)
          } else {
            ulam <- as.double(ulamtemp)
          }
          break
        }
      }
    } else {
      unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
        beta0 = double(nvars),
        beta = double(nvars),
        fit = double(nobs * nvars),
        PACKAGE = "sulnet"
      )
      f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
      storage.mode(f) <- "double"
      while (TRUE) {
        getlambdagauss <- .Fortran("getlambdagauss", nobs, nvars, nlam,
          ulam = ulam, f,
          y, pf, flmin, PACKAGE = "sulnet"
        )
        ulamtemp <- as.double(getlambdagauss$ulam)
        if (!anyNA(ulamtemp)) {
          flmin <- as.double(1)
          if (iglamPos) {
            ulam <- as.double(2 * ulamtemp)
          } else {
            ulam <- as.double(ulamtemp)
          }
          break
        }
      }
    }
  } else if (family == "binomial") {
    if (negOnly) {
      while (TRUE) {
        getlambdabinom <- .Fortran("getlambdabinom", nobs, nvars, nlam,
          ulam = ulam, x,
          y, pf, flmin, PACKAGE = "sulnet"
        )
        ulamtemp <- as.double(getlambdabinom$ulam)
        if (!anyNA(ulamtemp)) {
          flmin <- as.double(1)
          if (iglamPos) {
            ulam <- as.double(2 * ulamtemp)
          } else {
            ulam <- as.double(ulamtemp)
          }
          break
        }
      }
    } else {
      unifit <- .Fortran("loofit", nobs, nvars, x, y, loo,
        beta0 = double(nvars),
        beta = double(nvars),
        fit = double(nobs * nvars),
        PACKAGE = "sulnet"
      )
      f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
      storage.mode(f) <- "double"
      while (TRUE) {
        getlambdabinom <- .Fortran("getlambdabinom", nobs, nvars, nlam,
          ulam = ulam, f,
          y, pf, flmin, PACKAGE = "sulnet"
        )
        ulamtemp <- as.double(getlambdabinom$ulam)
        if (!anyNA(ulamtemp)) {
          flmin <- as.double(1)
          if (iglamPos) {
            ulam <- as.double(2 * ulamtemp)
          } else {
            ulam <- as.double(ulamtemp)
          }
          break
        }
      }
    }
  }
  ulam
}

getztR <- function(
  nobs, nvars, nzeta, ntau, x, y,
  uzeta, utau, pf, flminz, flmint, isd, intr, family
) {
  if (family == "gaussian") {
    getzt <- .Fortran(
      "getzt", nobs, nvars, nzeta, ntau,
      uzeta = uzeta, utau = utau,
      x, y, pf, flminz, flmint, isd, intr,
      diffz = double(1), difft = double(1),
      PACKAGE = "sulnet"
    )
    uzeta <- as.double(getzt$uzeta)
    utau <- as.double(getzt$utau)
    diffz <- as.double(getzt$diffz)
    difft <- as.double(getzt$difft)
    return(list(
      uzeta = uzeta, utau = utau,
      diffz = diffz, difft = difft
    ))
  }
}

getlaR <- function(
  nobs, nvars, nlambda, nalpha, x, y,
  ulambda, ualpha, pf, flmin, flmina, isd, intr
) {
  getla <- .Fortran("getla", nobs, nvars, nlambda, nalpha,
    ulambda = ulambda, ualpha = ualpha, x,
    y, pf, flmin, flmina, isd, intr, PACKAGE = "sulnet"
  )
  ulambda <- as.double(getla$ulambda)
  ualpha <- as.double(getla$ualpha)
  return(list(
    ulambda = ulambda, ualpha = ualpha
  ))
}

getsrR <- function(
  nobs, nvars, nsum, nrate, x, y,
  usum, urate, pf, flmins, flminr, isd, intr
) {
  getsr <- .Fortran("getsr", nobs, nvars, nsum, nrate,
    usum = usum, urate = urate, x,
    y, pf, flmins, flminr, isd, intr, PACKAGE = "sulnet"
  )
  usum <- as.double(getsr$usum)
  urate <- as.double(getsr$urate)
  return(list(
    usum = usum, urate = urate
  ))
}
