library(sulnet)
library(glmnet)
library(uniLasso)
## fit_methods_binom() fits the same seven methods as the binomial simulation:
## suni2D (three parameterizations), uniLasso, lasso, adalasso_uni and
## adalasso_lasso. summarize_sim() returns their means and standard errors.
source("inst/simulation/compare_binom.R")

## The CSVs written by prepare_uci_data.R and prepare_high_dim_data.R.
csv_list <- list.files(
  path = "./inst/binomial_true_data/data/csv", pattern = "\\.csv$",
  full.names = TRUE
)
names(csv_list) <- sub("\\.csv$", "", basename(csv_list))

## Read one dataset. The response is the 0/1 column "target"; every other
## column is a predictor. Numeric columns are used as they are. Text columns
## are expanded by model.matrix into treatment-contrast dummies: one per
## level, except the first level in alphabetical order.
read_binary_data <- function(path) {
  dat <- na.omit(read.csv(path))
  y <- dat$target
  dat$target <- NULL
  is_text <- vapply(dat, is.character, logical(1))
  x <- as.matrix(dat[!is_text])
  if (any(is_text)) {
    x <- cbind(x, model.matrix(~., data = dat[is_text])[, -1, drop = FALSE])
  }
  ## uniFit hands x to Fortran without coercing it.
  storage.mode(x) <- "double"
  stopifnot(all(y %in% c(0, 1)))
  list(x = x, y = as.numeric(y))
}

## Fit every method on the training rows and evaluate on the remaining rows.
## Predictors constant on the training rows are dropped first: their
## univariate fits are undefined, and they count as unselected. Returns an
## M x 2 matrix of
##   misclass  misclassification rate of yhat = 1{eta_hat > 0} on the
##             validation rows; eta_hat = 0 counts as half an error, as in
##             rpe_class()
##   support   number of nonzero slopes
## A method that fails gives NA in both columns.
eval_split_binom <- function(x, y, train, foldid, s = c("min", "1se")) {
  s <- match.arg(s)
  xtrain <- x[train, , drop = FALSE]
  keep <- apply(xtrain, 2, function(v) any(v != v[1]))
  fit <- fit_methods_binom(
    xtrain[, keep, drop = FALSE], y[train], x[-train, keep, drop = FALSE],
    foldid, s
  )
  yhat <- (1 + sign(fit$link)) / 2
  cbind(
    misclass = colMeans(abs(y[-train] - yhat)),
    support = colSums(fit$coef != 0)
  )
}

## nsplit random 50-50 splits of (x, y), stratified by class: half of each
## class (rounded down) is used for training and the rest for validation, so
## both halves keep the class balance. Within a split all methods share one
## fold assignment of the training rows, also stratified by class as in
## run_simulation_binom(), so every CV training set keeps both classes. Every
## split is drawn up front, so the result depends on the seed set before
## calling but not on ncores. With ncores > 1 the splits run on a local
## cluster, and progress messages are not shown.
##
## Returns a list with
##   mean, se  M x 2 matrices; se is sd / sqrt(number of finished splits)
##   table     the same as "mean (se)" strings
##   n_ok      number of splits in which each method finished
##   raw       list of the nsplit per-split M x 2 matrices
##   info      n, p (after dummy coding) and the number of rows with y = 1
compare_real_data_binom <- function(
  x, y, nsplit = 100, s = c("min", "1se"), nfolds = 10, ncores = 1
) {
  s <- match.arg(s)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  y <- as.numeric(drop(y))
  n <- nrow(x)
  stopifnot(length(y) == n, all(y %in% c(0, 1)), nsplit >= 2)

  class_rows <- split(seq_len(n), y)
  splits <- lapply(seq_len(nsplit), function(b) {
    train <- sort(unlist(lapply(class_rows, function(i) {
      i[sample.int(length(i), length(i) %/% 2)]
    }), use.names = FALSE))
    ntrain <- length(train)
    foldid <- integer(ntrain)
    foldid[order(y[train], runif(ntrain))] <- rep_len(seq_len(nfolds), ntrain)
    list(train = train, foldid = foldid)
  })
  run_split <- function(b) {
    warn <- character()
    res <- withCallingHandlers(
      eval_split_binom(x, y, splits[[b]]$train, splits[[b]]$foldid, s),
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
    parallel::clusterExport(
      cl, c("method_names_binom", "fit_methods_binom", "eval_split_binom")
    )
    out <- parallel::parLapplyLB(cl, seq_len(nsplit), run_split)
  } else {
    out <- lapply(seq_len(nsplit), run_split)
  }
  ## Each distinct warning is shown once per split. Rare dummy levels make
  ## sulnet's univariate logistic fits warn about separation once per CV fold.
  for (b in seq_len(nsplit)) {
    for (w in unique(out[[b]]$warn)) message("split ", b, ": ", w)
  }
  raw <- lapply(out, `[[`, "res")
  n_ok <- rowSums(!is.na(sapply(raw, function(r) r[, "misclass"])))
  c(summarize_sim(raw), list(
    n_ok = n_ok, raw = raw, info = c(n = n, p = ncol(x), n_y1 = sum(y))
  ))
}

## Run every dataset and save each result as it finishes; datasets with a
## saved result are skipped, so an interrupted run can be resumed.
## The seed depends only on the dataset's name, so adding or removing a CSV
## does not change the splits of the others.
## A dataset is skipped if the training half would have fewer than 3 rows of
## either class: glmnet's binomial fit needs at least 2 rows of each class in
## every CV training set. This skips Gode2011 (8 arrays, 4 per class).
ncores <- max(1L, parallel::detectCores() - 2L)
out_dir <- "./inst/binomial_true_data/results"
dir.create(out_dir, showWarnings = FALSE)
res_real <- list()
for (name in names(csv_list)) {
  out_file <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(out_file)) {
    next
  }
  dat <- read_binary_data(csv_list[[name]])
  n_min <- min(table(dat$y))
  if (n_min %/% 2 < 3) {
    message(
      "skipping ", name, ": the smaller class has only ", n_min,
      " rows, too few for a 50-50 split with cross-validation"
    )
    next
  }
  message("dataset ", name, " (n = ", nrow(dat$x), ", p = ", ncol(dat$x), ")")
  set.seed(sum(utf8ToInt(name)))
  res_real[[name]] <- compare_real_data_binom(
    dat$x, dat$y,
    nsplit = 100, ncores = ncores
  )
  print(res_real[[name]]$table)
  saveRDS(res_real[[name]], out_file)
}

## Collect every saved result into two tables in out_dir:
##   summary.csv  one row per dataset and method: mean and standard error of
##                the misclassification rate and the support, and the number
##                of splits in which the method finished
##   raw.csv      one row per dataset, split and method: the raw
##                misclassification rate and support
res_files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
names(res_files) <- sub("\\.rds$", "", basename(res_files))
res_all <- lapply(res_files, readRDS)
summary_tab <- do.call(rbind, lapply(names(res_all), function(name) {
  r <- res_all[[name]]
  data.frame(
    dataset = name, method = rownames(r$mean),
    misclass_mean = r$mean[, "misclass"], misclass_se = r$se[, "misclass"],
    support_mean = r$mean[, "support"], support_se = r$se[, "support"],
    n_ok = r$n_ok, row.names = NULL
  )
}))
raw_tab <- do.call(rbind, lapply(names(res_all), function(name) {
  raw <- res_all[[name]]$raw
  do.call(rbind, lapply(seq_along(raw), function(b) {
    data.frame(
      dataset = name, split = b, method = rownames(raw[[b]]),
      misclass = raw[[b]][, "misclass"], support = raw[[b]][, "support"],
      row.names = NULL
    )
  }))
}))
write.csv(summary_tab, file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(raw_tab, file.path(out_dir, "raw.csv"), row.names = FALSE)
