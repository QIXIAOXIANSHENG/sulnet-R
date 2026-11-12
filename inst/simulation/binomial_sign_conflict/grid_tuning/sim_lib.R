## Helpers for the binomial grid-tuning simulations (rounds 1-3, run on Argon
## in November 2026). Sourced by round*.R and final_30reps.R.
##
## Paths (environment variables, defaults for a run from the package root):
##   SIM_DIR   folder with compare_binom.R and binomial_simulation_setup.R
##             (default inst/simulation)
##   TUNE_DIR  this folder (default inst/simulation/binomial_sign_conflict/grid_tuning)
## On Argon every file was in one folder "scripts": SIM_DIR=scripts TUNE_DIR=scripts.

suppressPackageStartupMessages({
  library(sulnet)
  library(glmnet)
  library(uniLasso)
  library(MASS)
})
sim_dir <- Sys.getenv("SIM_DIR", "inst/simulation")
tune_dir <- Sys.getenv("TUNE_DIR", file.path(sim_dir, "binomial_sign_conflict", "grid_tuning"))
source(file.path(sim_dir, "compare_binom.R"))
source(file.path(sim_dir, "binomial_simulation_setup.R"))

## Command-line argument i, or a default when it is missing or empty.
cli <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(cli) >= i && nzchar(cli[i])) cli[i] else default

## ---------------------------------------------------------------------------
## Data: identical to inst/simulation/"binomial simulation.R" for replicate b,
## followed by the stratified fold draw of run_simulation_binom().
## ---------------------------------------------------------------------------
gen_data <- function(model, b, nfolds = 10, draw_folds = TRUE) {
  if (model == 1) {
    n <- 100; ntest <- 30; p <- 300; rho <- 0.8
    Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
    gamma <- rep(0, p)
    gamma[seq_len(floor(0.1 * p))] <- 1
    X <- matrix(0, n, p)
    set.seed(b)
    ybase <- rbinom(n, 1, 0.5)
    for (i in 1:n) {
      mu <- rep(0, p)
      if (ybase[i] == 1) mu[1:20] <- 0.5
      X[i, ] <- mvrnorm(1, mu = mu, Sigma = Sigma)
    }
    xtest <- mvrnorm(ntest, rep(0, p), Sigma = Sigma)
    y <- rbinom(n, 1, plogis(X %*% gamma))
    gamma0 <- 0
  } else {
    dat <- switch(as.character(model),
      "2" = simulation_setting("independent", n = 200, p = 1000, n_test = 100, seed = b),
      "3" = simulation_setting("ar1", n = 200, p = 1000, n_test = 100, seed = b),
      "4" = simulation_setting("block", n = 200, p = 1000, n_test = 100, seed = b),
      "5" = simulation_setting("imbalanced", n = 200, p = 1000, n_test = 100, seed = b),
      "6" = bh2011(n = 400, p = 1000, amplitude = 1, n_test = 100, seed = b),
      ## exploratory sign-conflict designs (not in "binomial simulation.R"):
      ## three equicorrelated blocks of 10 each carry one negative coefficient
      ## whose univariate slope is positive (beta_neg + rho * sum(beta_pos) > 0).
      ## Ids 8 and 9 here are models 8 and 9 of sign_conflict_simulation.R;
      ## id 7 here is an exploratory design that was not kept.
      "7" = simulate_binomial(n = 200, p = 1000,
        beta_active = rep(c(1, 1, 1, -1), 3), active = c(1:4, 11:14, 21:24),
        design = "block", rho = 0.7, block_size = 10L, n_test = 100, seed = b),
      "8" = simulate_binomial(n = 400, p = 1000,
        beta_active = rep(c(1, 1, 1, -1.5), 3), active = c(1:4, 11:14, 21:24),
        design = "block", rho = 0.7, block_size = 10L, n_test = 100, seed = b),
      "9" = simulate_binomial(n = 400, p = 1000,
        beta_active = rep(c(2, 2, -2), 3), active = c(1:3, 11:13, 21:23),
        design = "block", rho = 0.8, block_size = 10L, n_test = 100, seed = b)
    )
    X <- dat$train$x
    y <- dat$train$y
    xtest <- dat$test$x
    gamma <- dat$beta
    gamma0 <- dat$intercept
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  y <- as.numeric(drop(y))
  n <- nrow(X)
  foldid <- integer(n)
  ## draw_folds = FALSE leaves the RNG where run_simulation_binom() expects it
  if (draw_folds) foldid[order(y, runif(n))] <- rep_len(seq_len(nfolds), n)
  list(
    x = X, y = y, xtest = as.matrix(xtest), gamma = as.numeric(gamma),
    gamma0 = gamma0, eta_new = gamma0 + drop(xtest %*% gamma), foldid = foldid
  )
}

## ---------------------------------------------------------------------------
## Baselines (s = lambda.min): uniLasso without LOO, lasso, two adaptive lassos.
## ---------------------------------------------------------------------------
fit_baselines <- function(d) {
  x <- d$x; y <- d$y; newx <- d$xtest; foldid <- d$foldid
  p <- ncol(x)
  na_fit <- list(coef = rep(NA_real_, p), link = rep(NA_real_, nrow(newx)))
  safely <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      warning(name, " failed: ", conditionMessage(e), call. = FALSE)
      na_fit
    })
  }
  glmnet_out <- function(cvfit) {
    list(
      coef = as.numeric(coef(cvfit, s = "lambda.min"))[-1],
      link = as.numeric(predict(cvfit, newx, s = "lambda.min", type = "link"))
    )
  }
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
  fits$uniLasso <- safely("uniLasso", glmnet_out(
    cv.uniLasso(x, y, family = "binomial", foldid = foldid, loo = FALSE)
  ))
  fits$lasso <- safely("lasso", glmnet_out(
    cv.glmnet(x, y, family = "binomial", foldid = foldid)
  ))
  fits$adalasso_uni <- safely("adalasso_uni", adalasso_out(
    uniFit(x, y, loo = FALSE, family = "binomial")$beta
  ))
  fits$adalasso_lasso <- safely("adalasso_lasso", adalasso_out(fits$lasso$coef))
  eval_criteria_binom(
    do.call(cbind, lapply(fits, `[[`, "coef")),
    do.call(cbind, lapply(fits, `[[`, "link")),
    d$gamma, d$eta_new
  )
}

## ---------------------------------------------------------------------------
## cv.suni2D with the full-data and fold fits shared across CV losses.
## Mirrors sulnet::cv.suni2D (b11feb6). Returns, with rows = second coordinate
## (in the fit's order) and columns = first coordinate:
##   pen1, pen2      full-data grid
##   cvm, cvsd       per loss
##   crit            n2 x n1 x 7 array of test criteria per cell (NA if the
##                   full-data fit did not compute the cell)
##   sel             per loss, the (row, col) index of cv.min and cv.1se
## ---------------------------------------------------------------------------
cv_suni2D_multi <- function(d, method, grid = list(), opts = list(),
                            losses = c("loss", "Brier")) {
  x <- d$x; y <- d$y
  t0 <- proc.time()[["elapsed"]]
  fit <- do.call(suni2D, c(
    list(x = x, y = y, method = method, family = "binomial"), grid, opts
  ))
  method <- fit$method
  pen <- sulnet:::getpenalty(fit)
  pen.fit <- pen
  tol <- sqrt(.Machine$double.eps)
  if (method == "lambda_alpha" && abs(pen[[2]][1] - 0.5) <= tol) pen[[2]][1] <- 0.5
  if (method == "sum_rate" && abs(pen[[2]][1] - 1) <= tol) pen[[2]][1] <- 1
  pnames <- strsplit(method, "_")[[1]]
  foldid <- d$foldid
  outlist <- lapply(seq_len(max(foldid)), function(i) {
    w <- foldid == i
    do.call(suni2D, c(
      list(x = x[!w, , drop = FALSE], y = y[!w], method = method, family = "binomial"),
      setNames(pen, pnames), opts
    ))
  })
  uncomputed <- t(is.na(fit$gamma0))
  cvm <- cvsd <- sel <- list()
  for (loss in losses) {
    cs <- sulnet:::cv.sunilogpath(outlist, method, pen, x, y, foldid, loss)
    m <- cs$cvm
    s <- cs$cvsd
    m[uncomputed] <- NA
    s[uncomputed] <- NA
    cvm[[loss]] <- m
    cvsd[[loss]] <- s
    sel[[loss]] <- if (all(is.na(m))) NULL else {
      lm <- sulnet:::getmin2D(pen.fit$pen1, pen.fit$pen2, m, s, method)
      lapply(lm, function(l) c(row = match(l[[2]], pen.fit$pen2), col = match(l[[1]], pen.fit$pen1)))
    }
  }
  ## per-cell test criteria
  preds <- predict(fit, d$xtest, type = "link")
  n1 <- length(pen.fit$pen1)
  n2 <- length(pen.fit$pen2)
  crit <- array(NA_real_, c(n2, n1, length(criteria_names_binom)),
    dimnames = list(NULL, NULL, criteria_names_binom)
  )
  for (a in seq_len(n2)) {
    cr <- eval_criteria_binom(as.matrix(fit$gamma[[a]]), preds[[a]], d$gamma, d$eta_new)
    cr[is.na(fit$gamma0[, a]), ] <- NA
    crit[a, , ] <- cr
  }
  list(
    method = method, pen1 = pen.fit$pen1, pen2 = pen.fit$pen2,
    cvm = cvm, cvsd = cvsd, crit = crit, sel = sel,
    zeta1 = if (method == "zeta_tau") fit$zeta[1] else NA_real_,
    time = proc.time()[["elapsed"]] - t0
  )
}

## criteria of the selected cell
sel_crit <- function(r, loss = "loss", s = "cv.min") {
  ij <- r$sel[[loss]][[s]]
  if (is.null(ij)) return(setNames(rep(NA_real_, 7), criteria_names_binom))
  r$crit[ij[["row"]], ij[["col"]], ]
}

## Run one task (model, b) with a list of suni configurations. Each config is
## list(name, method, grid = function(d, zmax) list(...), opts = list(...)).
## zmax is the largest zeta of the default standardize = FALSE zeta grid.
run_task <- function(model, b, configs, out_file) {
  if (file.exists(out_file)) return(invisible("skip"))
  warn <- character()
  res <- withCallingHandlers(
    tryCatch({
      t0 <- proc.time()[["elapsed"]]
      d <- gen_data(model, b)
      base <- fit_baselines(d)
      f3 <- suni2D(d$x, d$y, family = "binomial", standardize = FALSE, nzeta = 3, ntau = 3)
      zmax <- f3$zeta[1] / (1 + 1e-9)
      suni <- list()
      for (cf in configs) {
        grid <- if (is.null(cf$grid)) list() else cf$grid(d, zmax)
        suni[[cf$name]] <- tryCatch(
          cv_suni2D_multi(d, cf$method, grid, cf$opts),
          error = function(e) {
            warning(cf$name, " failed: ", conditionMessage(e), call. = FALSE)
            NULL
          }
        )
      }
      list(model = model, b = b, zmax = zmax, base = base, suni = suni,
           time = proc.time()[["elapsed"]] - t0)
    }, error = function(e) list(model = model, b = b, error = conditionMessage(e))),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  res$warnings <- unique(warn)
  saveRDS(res, out_file)
  invisible(if (is.null(res$error)) "ok" else "error")
}
