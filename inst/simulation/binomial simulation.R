library(sulnet)
library(glmnet)
library(uniLasso)
library(MASS)
source("inst/simulation/compare_binom.R")
source("inst/simulation/binomial_simulation_setup.R")


# model 1. uniLasso case
n <- 100
ntest <- 30
p <- 300
rho <- 0.8
Sigma <- matrix(0, p, p)
for (i in 1:p) {
  for (j in 1:p) {
    Sigma[i, j] <- rho^abs(i - j)
  }
}
X <- matrix(0, n, p)
xtest <- matrix(0, ntest, p)
gamma <- rep(0, p)
s <- floor(0.1 * p)
gamma[1:s] <- 1
res <- list()
for (b in 1:200) {
  set.seed(b)
  ybase <- rbinom(n, 1, 0.5)
  for (i in 1:n) {
    mu <- rep(0, p)
    if (ybase[i] == 1) {
      mu[1:20] <- 0.5
    }
    X[i, ] <- mvrnorm(1, mu = mu, Sigma = Sigma)
  }
  xtest <- mvrnorm(ntest, rep(0, p), Sigma = Sigma)
  eta <- X %*% gamma
  y <- rbinom(n, 1, plogis(eta))
  res[[b]] <- run_simulation_binom(X, y, gamma, xtest, B = 1)$raw[1, , ]
}
res_model1 <- summarize_sim(res)
res_model1$table
saveRDS(res_model1, "./inst/simulation/binomial/res_model1.rds")

# model 2, independent data
n <- 200
p <- 1000
ntest <- 100
res <- list()
for (b in 1:200) {
  dat <- simulation_setting(
    "independent",
    n = n,
    p = p, n_test = ntest, seed = b
  )
  x <- dat$train$x
  y <- dat$train$y
  xtest <- dat$test$x
  gamma <- dat$beta
  res[[b]] <- run_simulation_binom(x, y, gamma, xtest, B = 1)$raw[1, , ]
}
res_model2 <- summarize_sim(res)
res_model2$table
saveRDS(res_model2, "./inst/simulation/binomial/res_model2.rds")


# model 3, AR(1)
n <- 200
p <- 1000
ntest <- 100
res <- list()
for (b in 1:200) {
  dat <- simulation_setting(
    "ar1",
    n = n,
    p = p, n_test = ntest, seed = b
  )
  x <- dat$train$x
  y <- dat$train$y
  xtest <- dat$test$x
  gamma <- dat$beta
  res[[b]] <- run_simulation_binom(x, y, gamma, xtest, B = 1)$raw[1, , ]
}
res_model3 <- summarize_sim(res)
res_model3$table
saveRDS(res_model3, "./inst/simulation/binomial/res_model3.rds")

# model 4, blocked
n <- 200
p <- 1000
ntest <- 100
res <- list()
for (b in 1:200) {
  dat <- simulation_setting(
    "block",
    n = n,
    p = p, n_test = ntest, seed = b
  )
  x <- dat$train$x
  y <- dat$train$y
  xtest <- dat$test$x
  gamma <- dat$beta
  res[[b]] <- run_simulation_binom(x, y, gamma, xtest, B = 1)$raw[1, , ]
}
res_model4 <- summarize_sim(res)
res_model4$table
saveRDS(res_model4, "./inst/simulation/binomial/res_model4.rds")


# model 5, imbalance
n <- 200
p <- 1000
ntest <- 100
res <- list()
for (b in 1:200) {
  dat <- simulation_setting(
    "imbalanced",
    n = n,
    p = p, n_test = ntest, seed = b
  )
  x <- dat$train$x
  y <- dat$train$y
  xtest <- dat$test$x
  gamma <- dat$beta
  res[[b]] <- run_simulation_binom(x, y, gamma, xtest, B = 1)$raw[1, , ]
}
res_model5 <- summarize_sim(res)
res_model5$table
saveRDS(res_model5, "./inst/simulation/binomial/res_model5.rds")
