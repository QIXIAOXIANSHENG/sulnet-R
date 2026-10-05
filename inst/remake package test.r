library(devtools)
clean_dll()
document()
load_all()
set.seed(1)
nobs = as.integer(100)
nvars = as.integer(10)
nzeta = as.integer(4)
ntau = as.integer(10)
x = matrix(rnorm(nobs * nvars), nrow = nobs, ncol = nvars)
storage.mode(x) = "double"
y = as.double(rnorm(nobs))
pf = as.double(rep(1, nobs))
uzeta = double(nzeta)
utau = double(ntau)
flminz = 0.01
flmint = 0.01
family = "gaussian"
alpha = NULL
getztR(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family, alpha)
diff(log(getztR(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family, alpha)$uzeta))

max(t(x) %*% y)
min(t(x) %*% y)
getlambdaR(nobs, nvars, nzeta, double(nzeta), x, y, pf, flminz, FALSE, family, TRUE, TRUE)
