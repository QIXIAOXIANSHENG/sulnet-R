getlambdaR <- function(nobs, nvars, nlam, ulam = ulam, x, y, pf, flmin, loo, family, negOnly, iglamPos) {
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


getztR <- function(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family) {
    if (family == "gaussian") {
        getzt <- .Fortran(
            "getzt", nobs, nvars, nzeta, ntau,
            uzeta = uzeta, utau = utau,
            x, y, pf, flminz, flmint, as.double(1), as.double(1),
            as.integer(0),
            PACKAGE = "sulnet"
        )
        uzeta <- as.double(getzt$uzeta)
        utau <- as.double(getzt$utau)
        return(list(uzeta = uzeta, utau = utau))
    }
}
