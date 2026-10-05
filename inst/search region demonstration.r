library(devtools)
clean_dll()
document()
load_all()

# data generation
set.seed(1)
nobs <- as.integer(1000)
p = as.integer(1000)
nvars = as.integer(10)
nzeta <- as.integer(100)
ntau <- as.integer(20)
x <- matrix(rnorm(nobs * p), nrow = nobs, ncol = p)
storage.mode(x) = "double"
y <- as.double(rnorm(nobs)) + x %*% c(rep(1, nvars), rep(0, p - nvars))
pf = as.double(rep(1, nobs))
uzeta = double(nzeta)
utau = double(ntau)
flminz <- 0.01
flmint <- 0.01
family = "gaussian"

# 1
alpha = NULL
lams = getztR(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family, alpha)
grid1 = as.matrix(expand.grid(lams$uzeta, lams$utau))
zeta1 = log(grid1[, 1])
tau1 = log(grid1[, 2] + grid1[, 1])

# 2
alpha = NULL
lams = getztR(nobs, nvars, nzeta, ntau, x, y, uzeta, utau, pf, flminz, flmint, family, alpha)
grid2 <- as.matrix(expand.grid(lams$uzeta, lams$utau))
zeta2 <- log(grid2[, 1])
tau2 = log(grid2[, 2] * grid2[, 1] / max(grid2[, 1]) + grid2[, 1])

# 3
alpha = seq(0.5, 0.99, length.out = ntau)
grid3 = as.matrix(expand.grid(lams$uzeta, alpha))
zeta3 <- log(grid3[, 1] * (1 - grid3[, 2]) * 2)
tau3 <- log(grid3[, 1] * grid3[, 2] * 2)

# on log scale
plot(zeta1, tau1, pch = 20, col = "red", xlab = "log(zeta)", ylab = "log(tau)", main = "Search region demonstration", xlim = range(c(zeta1, zeta2, zeta3)), ylim = range(c(tau1, tau2, tau3)))
points(zeta2, tau2, pch = 20, col = "blue")
points(zeta3, tau3, pch = 20, col = "green")
legend("bottomright", legend = c("1", "2", "3"), col = c("red", "blue", "green"), pch = 20)

# on original scale
plot(exp(zeta1), exp(tau1), pch = 20, col = "red", xlab = "zeta", ylab = "tau", main = "Search region demonstration", xlim = range(c(exp(zeta1), exp(zeta2), exp(zeta3))), ylim = range(c(exp(tau1), exp(tau2), exp(tau3))))
points(exp(zeta2), exp(tau2), pch = 20, col = "blue")
points(exp(zeta3), exp(tau3), pch = 20, col = "green")
legend("bottomright", legend = c("1", "2", "3"), col = c("red", "blue", "green"), pch = 20)



max(exp(zeta1) * exp(tau1))
max(exp(zeta2) * exp(tau2))
max(exp(zeta3) * exp(tau3))
