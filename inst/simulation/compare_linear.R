## Fixed-design simulation comparing suni2D (three penalty parameterizations),
## uniLasso, lasso, and two adaptive lassos in the linear model
## y = x %*% gamma + sigma * eps (no intercept in the truth).
##
## Usage:
##   source("inst/simulation/compare_linear.R")
##   set.seed(1)
##   res <- run_simulation(x, y, gamma, sigma, newx, B = 200)
##   res$table

library(sulnet)
library(uniLasso)
library(glmnet)

method_names <- c(
  "suni2D_zeta_tau", "suni2D_lambda_alpha", "suni2D_sum_rate",
  "uniLasso", "lasso", "adalasso_uni", "adalasso_lasso"
)
criteria_names <- c("sse", "rpe", "support", "sign", "tpr", "fpr")

## Fit every method on (x, y) with shared CV folds. Returns the p x M matrix of
## slopes (intercepts dropped) and the nrow(newx) x M matrix of predictions.
fit_methods <- function(x, y, newx, foldid, s = c("min", "1se")) {
  s <- match.arg(s)
  s_glmnet <- paste0("lambda.", s)
  s_suni <- paste0("cv.", s)
  p <- ncol(x)
  na_fit <- list(coef = rep(NA_real_, p), pred = rep(NA_real_, nrow(newx)))
  safely <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      warning(name, " failed: ", conditionMessage(e), call. = FALSE)
      na_fit
    })
  }
  suni_out <- function(method) {
    cvfit <- cv.suni2D(x, y, foldid = foldid, method = method)
    list(
      coef = as.numeric(coef(cvfit, s = s_suni)[[1]])[-1],
      pred = as.numeric(predict(cvfit, newx, s = s_suni)[[1]])
    )
  }
  glmnet_out <- function(cvfit) {
    list(
      coef = as.numeric(coef(cvfit, s = s_glmnet))[-1],
      pred = as.numeric(predict(cvfit, newx, s = s_glmnet))
    )
  }
  ## Zou (2006) weights w_j = 1 / |b_j|. A zero b_j gives an infinite weight,
  ## which glmnet treats as an excluded variable. The weighted penalty is
  ## defined on the original scale, so the reweighted fit is not standardized.
  adalasso_out <- function(b) {
    if (all(b == 0)) {
      return(list(coef = double(p), pred = rep(mean(y), nrow(newx))))
    }
    glmnet_out(cv.glmnet(x, y,
      foldid = foldid, penalty.factor = 1 / abs(b), standardize = FALSE
    ))
  }

  fits <- list()
  for (method in c("zeta_tau", "lambda_alpha", "sum_rate")) {
    name <- paste0("suni2D_", method)
    fits[[name]] <- safely(name, suni_out(method))
  }
  fits$uniLasso <- safely("uniLasso", glmnet_out(
    cv.uniLasso(x, y, foldid = foldid)
  ))
  fits$lasso <- safely("lasso", glmnet_out(cv.glmnet(x, y, foldid = foldid)))
  uni_slope <- drop(cov(x, y)) / apply(x, 2, var)
  fits$adalasso_uni <- safely("adalasso_uni", adalasso_out(uni_slope))
  fits$adalasso_lasso <- safely(
    "adalasso_lasso", adalasso_out(fits$lasso$coef)
  )

  fits <- fits[method_names]
  list(
    coef = do.call(cbind, lapply(fits, `[[`, "coef")),
    pred = do.call(cbind, lapply(fits, `[[`, "pred"))
  )
}

## Criteria for each column of gamma_hat / pred. Returns an M x 6 matrix.
##   sse     ||gamma - gamma_hat||_2^2
##   rpe     relative prediction error of Zou (2006),
##           E[(yhat - x'gamma)^2] / sigma^2, averaged over the rows of newx.
##           E[(y_test - yhat)^2] / sigma^2 equals rpe + 1.
##   support number of nonzero slopes
##   sign    proportion of the p coordinates with sign(gamma_hat) == sign(gamma),
##           where sign(0) = 0, so a missed signal or a false selection counts
##           as inconsistent
##   tpr     selected true signals / true signals
##   fpr     selected noise variables / noise variables
eval_criteria <- function(gamma_hat, pred, gamma, mu_new, sigma) {
  nz_hat <- gamma_hat != 0
  nz <- gamma != 0
  cbind(
    sse = colSums((gamma_hat - gamma)^2),
    rpe = colMeans((pred - mu_new)^2) / sigma^2,
    support = colSums(nz_hat),
    sign = colMeans(sign(gamma_hat) == sign(gamma)),
    tpr = if (any(nz)) colSums(nz_hat & nz) / sum(nz) else NA_real_,
    fpr = if (any(!nz)) colSums(nz_hat & !nz) / sum(!nz) else NA_real_
  )
}

## Replicate 1 uses the supplied y; replicates 2..B redraw
## y = x %*% gamma + sigma * rnorm(n) with x held fixed. Within a replicate all
## methods share one random fold assignment. Set the seed before calling for
## reproducible results.
##
## Returns a list with
##   mean, se  M x 6 matrices; se is sd / sqrt(number of finished replicates)
##   table     the same as "mean (se)" strings
##   n         number of finished (non-NA) replicates per method
##   raw       B x M x 6 array of per-replicate criteria
run_simulation <- function(
  x, y, gamma, sigma, newx, B = 200,
  s = c("min", "1se"), nfolds = 10
) {
  s <- match.arg(s)
  x <- as.matrix(x)
  newx <- as.matrix(newx)
  y <- drop(y)
  gamma <- as.numeric(gamma)
  n <- nrow(x)
  stopifnot(
    length(y) == n, length(gamma) == ncol(x), ncol(newx) == ncol(x),
    sigma > 0, B >= 1
  )
  mu <- drop(x %*% gamma)
  mu_new <- drop(newx %*% gamma)

  raw <- array(NA_real_,
    dim = c(B, length(method_names), length(criteria_names)),
    dimnames = list(NULL, method_names, criteria_names)
  )
  for (b in seq_len(B)) {
    yb <- if (b == 1L) y else mu + sigma * rnorm(n)
    foldid <- sample(rep(seq_len(nfolds), length.out = n))
    fit <- withCallingHandlers(
      fit_methods(x, yb, newx, foldid, s),
      warning = function(w) {
        message("replicate ", b, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    raw[b, , ] <- eval_criteria(fit$coef, fit$pred, gamma, mu_new, sigma)
    if (b %% 10L == 0L) message("finished replicate ", b, " / ", B)
  }

  n_ok <- apply(!is.na(raw[, , "sse", drop = FALSE]), 2, sum)
  avg <- apply(raw, c(2, 3), mean, na.rm = TRUE)
  se <- apply(raw, c(2, 3), sd, na.rm = TRUE) / sqrt(n_ok)
  table <- noquote(matrix(sprintf("%.4f (%.4f)", avg, se),
    nrow = nrow(avg), dimnames = dimnames(avg)
  ))
  list(mean = avg, se = se, table = table, n = n_ok, raw = raw)
}

summarize_sim <- function(results, digits = 4) {
  # Convert to 3D array: method × criterion × simulation
  arr <- simplify2array(lapply(results, as.matrix))

  # Mean and standard error
  mean_mat <- apply(arr, c(1, 2), mean, na.rm = TRUE)
  se_mat <- apply(arr, c(1, 2), function(x) {
    n <- sum(!is.na(x))
    if (n <= 1) {
      return(NA_real_)
    }
    sd(x, na.rm = TRUE) / sqrt(n)
  })

  # Formatted table: mean (se)
  table <- matrix(
    sprintf(
      paste0("%.", digits, "f (%.", digits, "f)"),
      mean_mat, se_mat
    ),
    nrow = nrow(mean_mat),
    dimnames = dimnames(mean_mat)
  )

  list(
    table = noquote(table),
    mean = mean_mat,
    se = se_mat
  )
}
