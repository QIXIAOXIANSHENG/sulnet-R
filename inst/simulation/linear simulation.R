library(sulnet)
library(MASS)
library(glmnet)
library(uniLasso)
source("inst/simulation/compare_linear.R")
set.seed(1)
nobs <- as.integer(100)
nvars <- as.integer(10)
nzeta <- as.integer(10)
ntau <- as.integer(5)
nlambda <- as.integer(10)
nalpha <- as.integer(5)
nsum <- as.integer(10)
nrate <- as.integer(5)
x <- matrix(rnorm(nobs * nvars), nrow = nobs, ncol = nvars)
storage.mode(x) <- "double"
y <- as.double(rnorm(nobs)) + x %*% rep(c(2, 0), each = nvars / 2)
ybinom <- rbinom(nobs, 1, 1 / (1 + exp(-y)))

unifit <- cv.uniLasso(x, y)
sunifit <- cv.suni2D(x, y, method = "lambda_alpha")


# model 1 Huang, AR(1)
set.seed(1)
p <- 200
p_true <- 15
sigma <- 1.5
rho <- 0.5
n <- 100
ntest <- 30
C <- matrix(NA, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    C[i, j] <- rho^(abs(i - j))
  }
}
C[(p_true + 1):p, 1:p_true] <- 0
C[1:p_true, (p_true + 1):p] <- 0
gamma <- c(rep(2.5, 5), rep(1.5, 5), rep(0.5, 5), rep(0, p - p_true))
x <- mvrnorm(n, rep(0, p), Sigma = C)
xtest <- mvrnorm(ntest, rep(0, p), Sigma = C)
y <- drop(x %*% gamma) + sigma * rnorm(n)
res_model1 <- run_simulation(x, y, gamma, sigma, xtest, B = 200)
res_model1$table
saveRDS(res_model1, "./inst/simulation/linear/res_model1.rds")


# model 2 Huang, AR(1) with larger rho
set.seed(2)
p <- 200
p_true <- 15
sigma <- 1.5
rho <- 0.95
n <- 100
ntest <- 30
C <- matrix(NA, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    C[i, j] <- rho^(abs(i - j))
  }
}
C[(p_true + 1):p, 1:p_true] <- 0
C[1:p_true, (p_true + 1):p] <- 0
gamma <- c(rep(2.5, 5), rep(1.5, 5), rep(0.5, 5), rep(0, p - p_true))
x <- mvrnorm(n, rep(0, p), Sigma = C)
xtest <- mvrnorm(ntest, rep(0, p), Sigma = C)
y <- drop(x %*% gamma) + sigma * rnorm(n)
res_model2 <- run_simulation(x, y, gamma, sigma, xtest, B = 200)
res_model2$table
saveRDS(res_model2, "./inst/simulation/linear/res_model2.rds")

# model 3 Huang, different setup
set.seed(3)
p <- 200
p_true <- 15
sigma <- 1.5
rho <- 0.5
n <- 100
ntest <- 30
z1 <- rnorm(n)
z2 <- rnorm(n)
z3 <- rnorm(n)
x1 <- z1 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
x2 <- z2 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
x3 <- z3 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
x4 <- matrix(rnorm(n * (p - 15)), n, p - 15)
x <- cbind(x1, x2, x3, x4)
z1test <- rnorm(ntest)
z2test <- rnorm(ntest)
z3test <- rnorm(ntest)
x1test <- z1test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
x2test <- z2test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
x3test <- z3test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
x4test <- matrix(rnorm(ntest * (p - 15)), ntest, p - 15)
xtest <- cbind(x1test, x2test, x3test, x4test)
gamma <- c(rep(1.5, p_true), rep(0, p - p_true))
y <- drop(x %*% gamma) + sigma * rnorm(n)
res_model3 <- run_simulation(x, y, gamma, sigma, xtest, B = 200)
res_model3$table
saveRDS(res_model3, "./inst/simulation/linear/res_model3.rds")

# model 4, uniLasso
set.seed(4)
n <- 300
ntest <- 100
p <- 1000
rho <- 0.5
p_true <- 30 # number of non-zero coefficients
sigma <- 0.1
Sigma <- matrix(rho, nrow = p, ncol = p)
diag(Sigma) <- 1
x <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
xtest <- mvrnorm(n = ntest, mu = rep(0, p), Sigma = Sigma)
B <- 200
res <- list()
for (i in seq(B)) {
  gamma <- c(runif(p_true), rep(0, p - p_true))
  y <- x %*% gamma + rnorm(n, 0, sigma)
  r <- run_simulation(x, y, gamma, sigma, xtest, B = 1)
  res[[i]] <- r$raw[1, , ]
}
res_model4 <- summarize_sim(res)
res_model4$table
saveRDS(res_model4, "./inst/simulation/linear/res_model4.rds")

res4 <- readRDS("./inst/simulation/linear/res_model4.rds")
