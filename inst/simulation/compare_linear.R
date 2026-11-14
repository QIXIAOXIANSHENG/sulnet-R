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

## Lambda paths of the glmnet-based methods (lasso, uniLasso, both adaptive
## lassos), changed on 2026-11-12; see inst/simulation/linear_fixed/README.md.
##   * glmnet's defaults end the path once the deviance ratio passes 0.999 or
##     improves by less than 1e-5 per step (glmnet.control(devmax, fdev)). In
##     high-SNR models (model 4: sigma = 0.1, R^2 about 0.9999) the path then
##     stops before the lambda CV wants. fit_methods() sets fdev = 0 and
##     devmax = 1 while it runs and restores the old values on exit.
##   * glmnet's default lambda.min.ratio is 0.01 when n < p; lasso and uniLasso
##     use GLMNET_LMR instead.
##   * The adaptive penalty factors span up to five decades, so with the
##     default path the variables with large weights never entered. The
##     adaptive lassos use ADA_NLAMBDA values down to ADA_LMR * lambda_max.
GLMNET_LMR <- 1e-4
ADA_LMR <- 1e-6
ADA_NLAMBDA <- 200L

## Adaptive lasso with Zou (2006) weights w_j = 1 / |b_j| for a first-stage
## estimate b; b_j = 0 means w_j = Inf, which drops variable j. It is fitted as
## a plain lasso on the rescaled columns x_j |b_j| with standardize = FALSE,
## which is the same problem as penalty.factor = w on the original scale but
## keeps lambda on the scale of lambda * sum_j |beta_j| / |b_j| (glmnet
## rescales penalty factors to sum to the number of variables), so that one
## lambda sequence means the same penalty in every CV fold.
ada_path <- function(x, y, b, lambda = NULL) {
  keep <- which(b != 0)
  scale <- abs(b[keep])
  xs <- x[, keep, drop = FALSE] * rep(scale, each = nrow(x))
  ## glmnet needs two columns; a constant column is never selected
  if (length(keep) < 2L) xs <- cbind(xs, 0)
  fit <- if (is.null(lambda)) {
    glmnet(xs, y,
      standardize = FALSE, lambda.min.ratio = ADA_LMR, nlambda = ADA_NLAMBDA
    )
  } else {
    glmnet(xs, y, standardize = FALSE, lambda = lambda)
  }
  list(fit = fit, keep = keep, scale = scale)
}
ada_predict <- function(af, newx, s) {
  xs <- newx[, af$keep, drop = FALSE] * rep(af$scale, each = nrow(newx))
  if (length(af$keep) < 2L) xs <- cbind(xs, 0)
  predict(af$fit, xs, s = s)
}
ada_coef <- function(af, s, p) {
  b <- double(p)
  b[af$keep] <- as.numeric(coef(af$fit, s = s))[1L + seq_along(af$keep)] *
    af$scale
  b
}

## Adaptive lasso tuned by honest 10-fold CV: in each fold the first stage is
## refit on the training folds, so the held-out fold never influences the
## weights. first(x, y, foldid) returns the first-stage slopes; foldid gives
## the folds of an inner CV, if the first stage needs one. b is the
## first stage on the full data. The CV error of each lambda is the mean
## squared error over the folds, with its standard error, lambda.min and
## lambda.1se defined as in cv.glmnet.
ada_cv <- function(x, y, newx, foldid, b, first, s) {
  p <- ncol(x)
  if (all(b == 0)) {
    return(list(coef = double(p), pred = rep(mean(y), nrow(newx))))
  }
  full <- ada_path(x, y, b)
  lambda <- full$fit$lambda
  nfolds <- max(foldid)
  cvraw <- matrix(NA_real_, nfolds, length(lambda))
  for (k in seq_len(nfolds)) {
    tr <- foldid != k
    ## the remaining folds, relabelled 1..nfolds-1, are the inner CV folds
    inner <- as.integer(factor(foldid[tr]))
    bk <- first(x[tr, , drop = FALSE], y[tr], inner)
    pred <- if (all(bk == 0)) {
      matrix(mean(y[tr]), sum(!tr), length(lambda))
    } else {
      ada_predict(ada_path(x[tr, , drop = FALSE], y[tr], bk, lambda),
        x[!tr, , drop = FALSE], lambda)
    }
    cvraw[k, ] <- colMeans((y[!tr] - pred)^2)
  }
  nk <- tabulate(foldid, nfolds)
  cvm <- colSums(nk * cvraw) / sum(nk)
  cvsd <- sqrt(colSums(nk * sweep(cvraw, 2, cvm)^2) / sum(nk) / (nfolds - 1))
  idmin <- which.min(cvm)
  lam <- if (s == "min") {
    lambda[idmin]
  } else {
    max(lambda[cvm <= cvm[idmin] + cvsd[idmin]])
  }
  list(
    coef = ada_coef(full, lam, p),
    pred = as.numeric(ada_predict(full, newx, lam))
  )
}

## Fit every method on (x, y) with shared CV folds. Returns the p x M matrix of
## slopes (intercepts dropped) and the nrow(newx) x M matrix of predictions.
## No random numbers are drawn: every CV, including the inner CV of the
## adaptive lassos' first stage, uses folds derived from foldid.
fit_methods <- function(x, y, newx, foldid, s = c("min", "1se")) {
  s <- match.arg(s)
  s_glmnet <- paste0("lambda.", s)
  s_suni <- paste0("cv.", s)
  p <- ncol(x)
  ctrl <- glmnet.control()
  glmnet.control(fdev = 0, devmax = 1)
  on.exit(glmnet.control(fdev = ctrl$fdev, devmax = ctrl$devmax), add = TRUE)
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
  ## first stages of the adaptive lassos
  first_uni <- function(x, y, foldid) drop(cov(x, y)) / apply(x, 2, var)
  first_lasso <- function(x, y, foldid) {
    as.numeric(coef(
      cv.glmnet(x, y, foldid = foldid, lambda.min.ratio = GLMNET_LMR),
      s = s_glmnet
    ))[-1]
  }

  fits <- list()
  for (method in c("zeta_tau", "lambda_alpha", "sum_rate")) {
    name <- paste0("suni2D_", method)
    fits[[name]] <- safely(name, suni_out(method))
  }
  fits$uniLasso <- safely("uniLasso", glmnet_out(
    cv.uniLasso(x, y, foldid = foldid, lambda.min.ratio = GLMNET_LMR)
  ))
  fits$lasso <- safely("lasso", glmnet_out(
    cv.glmnet(x, y, foldid = foldid, lambda.min.ratio = GLMNET_LMR)
  ))
  fits$adalasso_uni <- safely("adalasso_uni", ada_cv(
    x, y, newx, foldid, first_uni(x, y, foldid), first_uni, s
  ))
  ## the full-data first stage is the lasso fit above
  fits$adalasso_lasso <- safely("adalasso_lasso", if (anyNA(fits$lasso$coef)) {
    stop("the lasso fit failed")
  } else {
    ada_cv(x, y, newx, foldid, fits$lasso$coef, first_lasso, s)
  })

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
