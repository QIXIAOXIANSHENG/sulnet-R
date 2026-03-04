suni = sulnet2D(x,y,family="cox",method = "septhresh",alpha = seq(0.01,0.99,length.out = 21)[6])
y <- drop(y3)
x <- as.matrix(x3)
vnames <- colnames(x)
np <- dim(x)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
lambda = suni$lambdamat[[1]]
nlambda = length(lambda)
method = "septhresh"
family = "cox"
lambda.factor = 0.01
lambda2 = 0
pf = rep(1,nvars)
pf2 = rep(1, nvars)
dfmax = nvars + 1
pmax = min(dfmax * 1.2, nvars)
standardize = FALSE
intercept = TRUE
eps = 1e-08
maxit = 1e+09
lamPos = 0.1
loo = TRUE
alpha = 0.255
negOnly = FALSE
lam2 <- as.double(lambda2)
pf <- as.double(pf)
pf2 <- as.double(pf2)
isd <- as.integer(standardize)
intr <- as.integer(intercept)
eps <- as.double(eps)
dfmax <- as.integer(dfmax)
pmax <- as.integer(pmax)
lamPos <- as.double(sort(lamPos))
ignore_lamPos <- FALSE
if (!is.null(alpha)) {
alpha <- as.double(sort(alpha))
ignore_lamPos <- TRUE
}
jd <- as.integer(0)
nlam <- as.integer(nlambda)
flmin <- as.double(lambda.factor)
ulam <- as.double(lambda)
offset = NULL
y<- response.coxnet(y)
iglamPos <- as.integer(ignore_lamPos)
alpha <- as.double(alpha)
if (!is.matrix(y) || !all(match(c("time", "status"), dimnames(y)[[2]], 0))) stop("Cox model requires a matrix with columns 'time' (>0) and 'status'  (binary) as a response; a 'Surv' object suffices", call. = FALSE)
ty <- as.double(y[, "time"])
tevent <- as.double(y[, "status"])
ty <- ty + (1 - tevent) * 100 * .Machine$double.eps ## ties issue
if (any(ty <= 0)) stop("negative event times encountered;  not permitted for Cox family")
maxit <- as.integer(maxit)
weights <- rep(1,nobs)
weights <- as.double(weights)
if (is.null(offset)) {
offset <- ty * 0
} else {
storage.mode(offset) <- "double"
}
loo <- as.logical(loo)
save.image("./inst/cox tryout.RData")
unifit <- .Fortran("loofit", nobs, nvars, x, ty, loo,
    beta0 = double(nvars),
    beta = double(nvars),
    fit = double(nobs * nvars),
    PACKAGE = "sulnet"
)
f <- matrix(unifit$fit, nrow = nobs, ncol = nvars)
f[which(f == Inf)] <- 9e10
f[which(f == -Inf)] <- -9e10
storage.mode(f) <- "double"
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
fit <- .Fortran("coxsuninet",
          parm = as.double(1), nobs, nvars, f, ty, tevent, offset, as.double(rep(1, nobs)), jd, pf2, dfmax, as.integer(pmax), nlam, flmin, ulam, eps,
          maxit, isd, # need to get JHF to reverse these
          lmu = integer(1),
          ca = double(pmax * nlam),
          ia = integer(pmax),
          nin = integer(nlam),
          nulldev = double(1),
          dev = double(nlam),
          alm = double(nlam),
          nlp = integer(1),
          jerr = integer(1),
          iglamPos,
          as.double(0),
          as.double(alpha[1]), PACKAGE = "sulnet"
        )
fit$jerr
outlist <- coxgetcoef(fit, nvars, pmax, vnames)
save.image("./inst/cox tryout.RData")
load("./inst/cox tryout.RData")
