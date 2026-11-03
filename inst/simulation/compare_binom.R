## Fixed-design simulation comparing suni2D (three penalty parameterizations),
## uniLasso, lasso, and two adaptive lassos in the logistic model
## y ~ Bernoulli(plogis(gamma0 + x %*% gamma)).
##
## Usage:
##   source("inst/simulation/compare_binom.R")
##   set.seed(1)
##   y <- rbinom(nrow(x), 1, plogis(drop(x %*% gamma)))
##   res <- run_simulation_binom(x, y, gamma, newx, B = 200)
##   res$table
##
## Names carry a _binom suffix so this file can be sourced together with
## compare_linear.R.

library(sulnet)
library(uniLasso)
library(glmnet)

method_names_binom <- c(
  "suni2D_zeta_tau", "suni2D_lambda_alpha", "suni2D_sum_rate",
  "uniLasso", "lasso", "adalasso_uni", "adalasso_lasso"
)
criteria_names_binom <- c(
  "sse", "rpe_class", "rpe_kl", "support", "sign", "tpr", "fpr"
)

## Both relative errors take the true linear predictor eta (p = plogis(eta))
## and the fitted linear predictor eta_hat, a vector or a matrix with one column
## per method, and return one value per column. Working on the link scale keeps
## them finite: p is never formed as 1 - plogis(eta), which rounds to 0 once
## eta > 37, and log probabilities come from plogis(log.p = TRUE), so
## 0 * log(0) and log(0) never occur.

## Misclassification error of yhat = 1{eta_hat > 0} relative to the Bayes rule,
## sum P(y != yhat) / sum min(p, 1 - p) - 1. Where yhat and the Bayes class
## disagree the excess error is |2p - 1| = |tanh(eta / 2)|; a tie eta_hat = 0
## counts as a coin flip and adds half of that.
rpe_class <- function(eta, eta_hat) {
  eta_hat <- as.matrix(eta_hat)
  excess <- abs(tanh(eta / 2)) * (1 - sign(eta) * sign(eta_hat)) / 2
  colSums(excess) / sum(plogis(-abs(eta)))
}

## Bernoulli KL divergence relative to the entropy,
## sum KL(p || p_hat) / sum H(p). The expected test log-loss,
## -E[y log p_hat + (1 - y) log(1 - p_hat)], divided by sum H(p) equals
## rpe_kl + 1.
rpe_kl <- function(eta, eta_hat) {
  eta_hat <- as.matrix(eta_hat)
  p <- plogis(eta)
  q <- plogis(-eta)
  log_p <- plogis(eta, log.p = TRUE)
  log_q <- plogis(-eta, log.p = TRUE)
  kl <- p * (log_p - plogis(eta_hat, log.p = TRUE)) +
    q * (log_q - plogis(-eta_hat, log.p = TRUE))
  colSums(kl) / -sum(p * log_p + q * log_q)
}

## Fit every method on (x, y) with shared CV folds; y is coded 0/1. Arguments in
## ... (for example nzeta, ntau, eps) go to cv.suni2D only. Returns the p x M
## matrix of slopes (intercepts dropped) and the nrow(newx) x M matrix of
## linear predictors.
fit_methods_binom <- function(x, y, newx, foldid, s = c("min", "1se"), ...) {
  s <- match.arg(s)
  s_glmnet <- paste0("lambda.", s)
  s_suni <- paste0("cv.", s)
  p <- ncol(x)
  na_fit <- list(coef = rep(NA_real_, p), link = rep(NA_real_, nrow(newx)))
  safely <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      warning(name, " failed: ", conditionMessage(e), call. = FALSE)
      na_fit
    })
  }
  suni_out <- function(method) {
    cvfit <- cv.suni2D(x, y,
      foldid = foldid, method = method, family = "binomial", ...
    )
    list(
      coef = as.numeric(coef(cvfit, s = s_suni)[[1]])[-1],
      link = as.numeric(predict(cvfit, newx, s = s_suni, type = "link")[[1]])
    )
  }
  glmnet_out <- function(cvfit) {
    list(
      coef = as.numeric(coef(cvfit, s = s_glmnet))[-1],
      link = as.numeric(predict(cvfit, newx, s = s_glmnet, type = "link"))
    )
  }
  ## Zou (2006) weights w_j = 1 / |b_j|. A zero b_j gives an infinite weight,
  ## which glmnet treats as an excluded variable. The weighted penalty is
  ## defined on the original scale, so the reweighted fit is not standardized.
  adalasso_out <- function(b) {
    if (all(b == 0)) {
      return(list(coef = double(p), link = rep(qlogis(mean(y)), nrow(newx))))
    }
    glmnet_out(cv.glmnet(x, y,
      family = "binomial", foldid = foldid,
      penalty.factor = 1 / abs(b), standardize = FALSE
    ))
  }

  fits <- list()
  for (method in c("zeta_tau", "lambda_alpha", "sum_rate")) {
    name <- paste0("suni2D_", method)
    fits[[name]] <- safely(name, suni_out(method))
  }
  fits$uniLasso <- safely("uniLasso", glmnet_out(
    cv.uniLasso(x, y, family = "binomial", foldid = foldid)
  ))
  fits$lasso <- safely("lasso", glmnet_out(
    cv.glmnet(x, y, family = "binomial", foldid = foldid)
  ))
  fits$adalasso_uni <- safely("adalasso_uni", adalasso_out(
    uniFit(x, y, loo = FALSE, family = "binomial")$beta
  ))
  fits$adalasso_lasso <- safely(
    "adalasso_lasso", adalasso_out(fits$lasso$coef)
  )

  fits <- fits[method_names_binom]
  list(
    coef = do.call(cbind, lapply(fits, `[[`, "coef")),
    link = do.call(cbind, lapply(fits, `[[`, "link"))
  )
}

## Criteria for each column of gamma_hat / link. Returns an M x 7 matrix.
##   sse       ||gamma - gamma_hat||_2^2 (slopes only)
##   rpe_class see rpe_class() above, over the rows of newx
##   rpe_kl    see rpe_kl() above, over the rows of newx
##   support   number of nonzero slopes
##   sign      proportion of the p coordinates with
##             sign(gamma_hat) == sign(gamma), where sign(0) = 0, so a missed
##             signal or a false selection counts as inconsistent
##   tpr       selected true signals / true signals
##   fpr       selected noise variables / noise variables
eval_criteria_binom <- function(gamma_hat, link, gamma, eta_new) {
  nz_hat <- gamma_hat != 0
  nz <- gamma != 0
  cbind(
    sse = colSums((gamma_hat - gamma)^2),
    rpe_class = rpe_class(eta_new, link),
    rpe_kl = rpe_kl(eta_new, link),
    support = colSums(nz_hat),
    sign = colMeans(sign(gamma_hat) == sign(gamma)),
    tpr = if (any(nz)) colSums(nz_hat & nz) / sum(nz) else NA_real_,
    fpr = if (any(!nz)) colSums(nz_hat & !nz) / sum(!nz) else NA_real_
  )
}

## Replicate 1 uses the supplied 0/1 response y; replicates 2..B redraw
## y = rbinom(n, 1, plogis(gamma0 + x %*% gamma)) with x held fixed. Within a
## replicate all methods share one fold assignment, stratified by class:
## observations are sorted by class, shuffled within class, and dealt to the
## folds in turn, so every training fold keeps both classes. Set the seed
## before calling for reproducible results. Arguments in ... go to cv.suni2D,
## for example grid sizes (nzeta/ntau, nlambda/nalpha, nsum/nrate; each method
## uses its own pair), eps or maxit.
##
## Returns a list with
##   mean, se  M x 7 matrices; se is sd / sqrt(number of finished replicates)
##   table     the same as "mean (se)" strings
##   n         number of finished (non-NA) replicates per method
##   raw       B x M x 7 array of per-replicate criteria
run_simulation_binom <- function(
  x, y, gamma, newx, gamma0 = 0, B = 200,
  s = c("min", "1se"), nfolds = 10, ...
) {
  s <- match.arg(s)
  x <- as.matrix(x)
  ## uniFit hands x to Fortran without coercing it.
  storage.mode(x) <- "double"
  newx <- as.matrix(newx)
  y <- drop(y)
  gamma <- as.numeric(gamma)
  n <- nrow(x)
  stopifnot(
    length(y) == n, all(y %in% c(0, 1)), length(gamma) == ncol(x),
    ncol(newx) == ncol(x), length(gamma0) == 1L, B >= 1
  )
  y <- as.numeric(y)
  eta <- gamma0 + drop(x %*% gamma)
  eta_new <- gamma0 + drop(newx %*% gamma)

  raw <- array(NA_real_,
    dim = c(B, length(method_names_binom), length(criteria_names_binom)),
    dimnames = list(NULL, method_names_binom, criteria_names_binom)
  )
  for (b in seq_len(B)) {
    yb <- if (b == 1L) y else rbinom(n, 1, plogis(eta))
    foldid <- integer(n)
    foldid[order(yb, runif(n))] <- rep_len(seq_len(nfolds), n)
    fit <- withCallingHandlers(
      fit_methods_binom(x, yb, newx, foldid, s, ...),
      warning = function(w) {
        message("replicate ", b, ": ", conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    raw[b, , ] <- eval_criteria_binom(fit$coef, fit$link, gamma, eta_new)
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
