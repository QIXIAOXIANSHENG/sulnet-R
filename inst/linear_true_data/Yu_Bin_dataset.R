library(sulnet)
library(MASS)
library(glmnet)
library(uniLasso)
## fit_methods() and summarize_sim() fit and summarize the same seven methods
## as the simulation: suni2D (three parameterizations), uniLasso, lasso,
## adalasso_uni and adalasso_lasso.
source("inst/simulation/compare_linear.R")
csv_list = list.files(path = "./inst/linear_true_data/data", pattern = "csv", full.names = TRUE)
names(csv_list) <- sub("\\.csv$", "", basename(csv_list))

## Read one dataset. The response is the "target" column; every other column
## is a predictor, expanded by model.matrix (character columns become dummies)
## with the intercept column removed.
read_real_data <- function(path) {
  dat <- na.omit(read.csv(path))
  x <- model.matrix(target ~ ., data = dat)[, -1, drop = FALSE]
  list(x = x, y = dat$target)
}

## Fit every method on the training rows and evaluate on the remaining rows.
## Predictors constant on the training rows are dropped first: their
## univariate slopes are undefined, and they count as unselected.
## With scale_y = TRUE the methods are fit to y standardized by the training
## mean and sd, and predictions are mapped back, so mse is on the original
## scale. The glmnet-based methods are equivariant to this, but the
## automatic suni2D penalty grid scales with var(y) rather than sd(y), so on
## an unscaled y it lies above (large var(y)) or below (small var(y)) the
## useful range. Returns an M x 2 matrix of
##   mse     mean squared prediction error on the validation rows
##   support number of nonzero slopes
eval_split <- function(
  x, y, train, foldid, s = c("min", "1se"), scale_y = TRUE
) {
  s <- match.arg(s)
  xtrain <- x[train, , drop = FALSE]
  keep <- apply(xtrain, 2, function(v) any(v != v[1]))
  center <- if (scale_y) mean(y[train]) else 0
  scale <- if (scale_y) sd(y[train]) else 1
  fit <- fit_methods(
    xtrain[, keep, drop = FALSE], (y[train] - center) / scale,
    x[-train, keep, drop = FALSE], foldid, s
  )
  pred <- center + scale * fit$pred
  cbind(
    mse = colMeans((pred - y[-train])^2),
    support = colSums(fit$coef != 0)
  )
}

## nsplit random 50-50 splits of (x, y): floor(n / 2) training rows, the rest
## for validation. Within a split all methods share one random fold
## assignment of the training rows. Every split is drawn up front, so the
## result depends on the seed set before calling but not on ncores. With
## ncores > 1 the splits run on a local cluster, and progress messages are
## not shown.
##
## Returns a list with
##   mean, se  M x 2 matrices; se is sd / sqrt(number of finished splits)
##   table     the same as "mean (se)" strings
##   raw       list of the nsplit per-split M x 2 matrices
compare_real_data <- function(
  x, y, nsplit = 100, s = c("min", "1se"), nfolds = 10, scale_y = TRUE,
  ncores = 1
) {
  s <- match.arg(s)
  x <- as.matrix(x)
  y <- drop(y)
  n <- nrow(x)
  stopifnot(length(y) == n, nsplit >= 2)
  ntrain <- n %/% 2

  splits <- lapply(seq_len(nsplit), function(b) {
    list(
      train = sample(n, ntrain),
      foldid = sample(rep(seq_len(nfolds), length.out = ntrain))
    )
  })
  run_split <- function(b) {
    warn <- character()
    res <- withCallingHandlers(
      eval_split(x, y, splits[[b]]$train, splits[[b]]$foldid, s, scale_y),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (b %% 10L == 0L) message("finished split ", b, " / ", nsplit)
    list(res = res, warn = warn)
  }
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, {
      library(sulnet)
      library(glmnet)
      library(uniLasso)
    })
    parallel::clusterExport(cl, c("method_names", "fit_methods", "eval_split"))
    out <- parallel::parLapplyLB(cl, seq_len(nsplit), run_split)
  } else {
    out <- lapply(seq_len(nsplit), run_split)
  }
  for (b in seq_len(nsplit)) {
    for (w in out[[b]]$warn) message("split ", b, ": ", w)
  }
  raw <- lapply(out, `[[`, "res")
  c(summarize_sim(raw), list(raw = raw))
}

## Run every dataset and save each result as it finishes.
ncores <- max(1L, parallel::detectCores() - 2L)
out_dir <- "./inst/linear_true_data/results"
dir.create(out_dir, showWarnings = FALSE)
res_real <- list()
for (i in seq_along(csv_list)) {
  name <- names(csv_list)[i]
  message("dataset ", name)
  dat <- read_real_data(csv_list[[i]])
  set.seed(i)
  res_real[[name]] <- compare_real_data(
    dat$x, dat$y, nsplit = 100, ncores = ncores
  )
  print(res_real[[name]]$table)
  saveRDS(res_real[[name]], file.path(out_dir, paste0(name, ".rds")))
}
