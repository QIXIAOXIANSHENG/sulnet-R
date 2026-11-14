## Linear models 1-10 with the tuned methods of fit_methods_tuned.R
## (2026-11-12): suni2D with standardize = FALSE and TRUE, uniLasso with and
## without LOO, lasso and the honest-CV adaptive lassos of compare_linear.R.
## Needs sulnet with the 2026-11-12 Gaussian solver fix. See README.md.
##
## Usage, from the package root (the folder that contains inst/):
##   Rscript inst/simulation/linear_tuned/linear_tuned_simulation.R \
##     [models] [reps] [ncores] [out_dir]
## Defaults: models = 1:10, reps = 1:200, ncores = 1,
##   out_dir = inst/simulation/linear_tuned
## Example, models 6 and 7 with replicates 1-10 on 8 cores:
##   Rscript <this file> 6:7 1:10 8
##
## Models 1-5 use the data, replicate responses and CV folds of
## "linear simulation.R" and inst/simulation/linear_fixed (same draws, same
## order, same seeds). Models 6-10 are sign-conflict designs: model 6 is a
## small block design, models 7-10 have the covariates and coefficients of the
## binomial sign-conflict models 7-10 (inst/simulation/binomial_sign_conflict,
## simulate_binomial(seed = b)) with a Gaussian response, SNR = 3.
## The fits draw no random numbers. ncores > 1 uses a PSOCK cluster, which
## works on Windows.
##
## Each replicate is saved as <out_dir>/reps/m<k>_b<bbb>.rds. A replicate whose
## file already exists is skipped, so an interrupted run resumes where it
## stopped when rerun with the same arguments. After each model,
## <out_dir>/res_model<k>.rds is rebuilt from every saved replicate of it
## (mean, se, table, n, raw, reps, seconds) and <out_dir>/linear_tuned_tables.pdf
## from every res_model<k>.rds in out_dir.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) args[i] else default
}
models <- eval(parse(text = arg(1, "1:10")))
reps <- eval(parse(text = arg(2, "1:200")))
ncores <- as.integer(arg(3, "1"))
out_dir <- arg(4, "inst/simulation/linear_tuned")
sim_dir <- "inst/simulation"
code_files <- file.path(sim_dir, c("compare_linear.R",
  "linear_tuned/fit_methods_tuned.R", "binomial_simulation_setup.R"))

suppressPackageStartupMessages({
  library(MASS)
  library(sulnet)
  library(glmnet)
  library(uniLasso)
})
for (f in code_files) source(f)
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------
draw_foldid <- function(n, nfolds = 10L) {
  sample(rep(seq_len(nfolds), length.out = n))
}

## Models 1-4 (fixed x): model_data(k, B) makes the draws of
## "linear simulation.R" in their order there and returns
##   shared  data common to all replicates of model k
##   reps    one list per replicate b = 1..B with the rest of
##           x, y, gamma, sigma, xtest, foldid
## Models 1-3: replicate 1 is the y drawn with x, replicates 2..B redraw
## y = x gamma + sigma eps, as in run_simulation().
fixed_x_reps <- function(x, y, gamma, sigma, B) {
  mu <- drop(x %*% gamma)
  lapply(seq_len(B), function(b) {
    yb <- if (b == 1L) drop(y) else mu + sigma * rnorm(nrow(x))
    list(y = yb, foldid = draw_foldid(nrow(x)))
  })
}

model_data <- function(k, B) {
  if (k %in% 1:2) {
    set.seed(k)
    p <- 200; p_true <- 15; sigma <- 1.5; n <- 100; ntest <- 30
    rho <- c(0.5, 0.95)[k]
    C <- rho^abs(outer(1:p, 1:p, "-"))
    C[(p_true + 1):p, 1:p_true] <- 0
    C[1:p_true, (p_true + 1):p] <- 0
    gamma <- c(rep(2.5, 5), rep(1.5, 5), rep(0.5, 5), rep(0, p - p_true))
    x <- mvrnorm(n, rep(0, p), Sigma = C)
    xtest <- mvrnorm(ntest, rep(0, p), Sigma = C)
    y <- drop(x %*% gamma) + sigma * rnorm(n)
    return(list(
      shared = list(x = x, xtest = xtest, gamma = gamma, sigma = sigma),
      reps = fixed_x_reps(x, y, gamma, sigma, B)
    ))
  }
  if (k == 3) {
    set.seed(3)
    p <- 200; p_true <- 15; sigma <- 1.5; n <- 100; ntest <- 30
    z1 <- rnorm(n); z2 <- rnorm(n); z3 <- rnorm(n)
    x1 <- z1 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
    x2 <- z2 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
    x3 <- z3 + matrix(rnorm(n * 5, 0, 0.1), n, 5)
    x4 <- matrix(rnorm(n * (p - 15)), n, p - 15)
    x <- cbind(x1, x2, x3, x4)
    z1test <- rnorm(ntest); z2test <- rnorm(ntest); z3test <- rnorm(ntest)
    x1test <- z1test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
    x2test <- z2test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
    x3test <- z3test + matrix(rnorm(ntest * 5, 0, 0.1), ntest, 5)
    x4test <- matrix(rnorm(ntest * (p - 15)), ntest, p - 15)
    xtest <- cbind(x1test, x2test, x3test, x4test)
    gamma <- c(rep(1.5, p_true), rep(0, p - p_true))
    y <- drop(x %*% gamma) + sigma * rnorm(n)
    return(list(
      shared = list(x = x, xtest = xtest, gamma = gamma, sigma = sigma),
      reps = fixed_x_reps(x, y, gamma, sigma, B)
    ))
  }
  if (k == 4) {
    ## fixed x; each replicate draws new coefficients and a new y
    set.seed(4)
    n <- 300; ntest <- 100; p <- 1000; rho <- 0.5; p_true <- 30; sigma <- 0.1
    Sigma <- matrix(rho, nrow = p, ncol = p)
    diag(Sigma) <- 1
    x <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
    xtest <- mvrnorm(n = ntest, mu = rep(0, p), Sigma = Sigma)
    reps <- lapply(seq_len(B), function(b) {
      gamma <- c(runif(p_true), rep(0, p - p_true))
      y <- drop(x %*% gamma + rnorm(n, 0, sigma))
      list(y = y, gamma = gamma, foldid = draw_foldid(n))
    })
    return(list(shared = list(x = x, xtest = xtest, sigma = sigma),
      reps = reps))
  }
  stop("model_data() covers models 1-4; use rep_data() for model ", k)
}

## Models 5-10: every replicate depends only on its seed b, so rep_data(k, b)
## runs on the worker. It returns x, y, gamma, sigma, xtest, foldid.
## Models 7-10: the binomial sign-conflict models of
## inst/simulation/binomial_sign_conflict (same simulate_binomial() call, so
## the same x, xtest and coefficients for seed b), with y = x gamma + sigma eps
## drawn after the binomial draws and sigma = sd(x'gamma) / sqrt(3).
conflict_defs <- list(
  ## blocks of 10, rho 0.5; (1, 1, 1, -1) in blocks 1-3: conflict slope +0.5
  "7" = list(design = "block", rho = 0.5, block_size = 10L,
    beta = rep(c(1, 1, 1, -1), 3), active = c(1:4, 11:14, 21:24)),
  ## blocks of 10, rho 0.7; (1, 1, 1, -1.5) in blocks 1-3: conflict slope +0.6
  "8" = list(design = "block", rho = 0.7, block_size = 10L,
    beta = rep(c(1, 1, 1, -1.5), 3), active = c(1:4, 11:14, 21:24)),
  ## blocks of 10, rho 0.8; (2, 2, -2) in blocks 1-3: conflict slope +1.2
  "9" = list(design = "block", rho = 0.8, block_size = 10L,
    beta = rep(c(2, 2, -2), 3), active = c(1:3, 11:13, 21:23)),
  ## AR(1) rho 0.8; (1, 1, 1, -1.5) at x1-x4, x101-x104, x201-x204: slope +0.45
  "10" = list(design = "ar1", rho = 0.8, block_size = 10L,
    beta = rep(c(1, 1, 1, -1.5), 3), active = c(1:4, 101:104, 201:204))
)
CONFLICT_SNR <- 3

rep_data <- function(k, b) {
  if (k == 5) {
    ## uniLasso home court: new x, coefficients and y for seed b
    n <- 100; ntest <- 100; p <- 30; rho <- 0.8; sparsity <- 0.2
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    set.seed(b)
    X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
    xtest <- mvrnorm(ntest, mu = rep(0, p), Sigma = Sigma)
    gamma <- rep(0, p)
    nonzero_idx <- sample(1:p, size = floor(sparsity * p))
    gamma[nonzero_idx] <- runif(length(nonzero_idx), min = 0.5, max = 1)
    z_prime <- rnorm(n)
    sigma_prime <- sqrt(var(X %*% gamma))
    y_prime <- X %*% gamma + sigma_prime[1, 1] * z_prime
    gamma_uni <- uniFit(X, y_prime)$beta
    z <- rnorm(n)
    gamma_effective <- gamma_uni * gamma
    sigma <- sqrt(var(X %*% gamma_effective))[1, 1]
    y <- drop(X %*% gamma_effective + sigma * z)
    return(list(x = X, xtest = xtest, y = y, gamma = gamma_effective,
      sigma = sigma, foldid = draw_foldid(n)))
  }
  if (k == 6) {
    ## 3 blocks of 4 equicorrelated (rho 0.5) signals (1, 1, 1, -1), the other
    ## 488 covariates independent: the -1 has univariate slope +0.5.
    ## var(x'gamma) = 12, sigma = 2 (SNR 3).
    n <- 200; ntest <- 100; p <- 500; rho <- 0.5; sigma <- 2
    gen <- function(m) {
      z <- matrix(rnorm(m * p), m)
      for (g in 0:2) {
        idx <- 4 * g + 1:4
        z[, idx] <- sqrt(rho) * rnorm(m) + sqrt(1 - rho) * z[, idx]
      }
      z
    }
    gamma <- c(rep(c(1, 1, 1, -1), 3), rep(0, p - 12))
    set.seed(b)
    x <- gen(n)
    xtest <- gen(ntest)
    y <- drop(x %*% gamma) + sigma * rnorm(n)
    return(list(x = x, xtest = xtest, y = y, gamma = gamma, sigma = sigma,
      foldid = draw_foldid(n)))
  }
  if (as.character(k) %in% names(conflict_defs)) {
    m <- conflict_defs[[as.character(k)]]
    d <- simulate_binomial(n = 400, p = 1000, beta_active = m$beta,
      active = m$active, design = m$design, rho = m$rho,
      block_size = m$block_size, n_test = 100, seed = b)
    sigma <- d$settings$signal_sd / sqrt(CONFLICT_SNR)
    y <- drop(d$train$x %*% d$beta) + sigma * rnorm(400)
    return(list(x = d$train$x, xtest = d$test$x, y = y, gamma = d$beta,
      sigma = sigma, foldid = draw_foldid(400)))
  }
  stop("unknown model ", k)
}

## ---------------------------------------------------------------------------
## One replicate: fit, score, save. Returns NULL, or the error message.
## ---------------------------------------------------------------------------
rep_file <- function(k, b) file.path(rep_dir, sprintf("m%d_b%03d.rds", k, b))

run_rep <- function(task, shared) {
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    d <- if (is.null(task$data)) rep_data(task$k, task$b) else c(task$data, shared)
    fit <- withCallingHandlers(
      fit_methods_tuned(d$x, d$y, d$xtest, d$foldid, s = "min"),
      warning = function(w) {
        message("model ", task$k, " replicate ", task$b, ": ",
          conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    crit <- eval_criteria(fit$coef, fit$pred, d$gamma,
      drop(d$xtest %*% d$gamma), d$sigma)
    saveRDS(list(crit = crit, seconds = proc.time()[["elapsed"]] - t0),
      task$file)
    NULL
  }, error = function(e) {
    sprintf("model %d replicate %d: %s", task$k, task$b, conditionMessage(e))
  })
}

## ---------------------------------------------------------------------------
## Tables of every model saved in out_dir
## ---------------------------------------------------------------------------
source(file.path(sim_dir, "table_simulation.R"))
write_tables <- function() {
  write_sim_tables(read_sim_results(out_dir),
    file = file.path(out_dir, "linear_tuned_tables.pdf"), per_page = 2,
    main = "Linear simulations with the tuned suni2D and uniLasso settings",
    note = paste("Mean (standard error) across replicates.",
      "_std: suni2D with standardize = TRUE; _loo: uniLasso with LOO."),
    captions = c(
      "Model 1" = "Huang, AR(1) rho = 0.5; n = 100, p = 200, 15 nonzero; sigma = 1.5",
      "Model 2" = "Huang, AR(1) rho = 0.95; n = 100, p = 200, 15 nonzero; sigma = 1.5",
      "Model 3" = "Huang, 3 groups of 5 near-collinear signals; n = 100, p = 200; sigma = 1.5",
      "Model 4" = "uniLasso, equicorrelated rho = 0.5; n = 300, p = 1000, 30 nonzero ~ U(0, 1); sigma = 0.1",
      "Model 5" = "uniLasso home court, AR(1) rho = 0.8; n = 100, p = 30, 20% nonzero; SNR = 1",
      "Model 6" = "Sign conflict: 3 blocks of 4 (rho = 0.5), (1, 1, 1, -1); n = 200, p = 500; SNR = 3",
      "Model 7" = "Sign conflict: blocks of 10 (rho = 0.5), (1, 1, 1, -1) in 3 blocks; n = 400, p = 1000; SNR = 3",
      "Model 8" = "Sign conflict: blocks of 10 (rho = 0.7), (1, 1, 1, -1.5) in 3 blocks; n = 400, p = 1000; SNR = 3",
      "Model 9" = "Sign conflict: blocks of 10 (rho = 0.8), (2, 2, -2) in 3 blocks; n = 400, p = 1000; SNR = 3",
      "Model 10" = "Sign conflict: AR(1) rho = 0.8, (1, 1, 1, -1.5) at 3 places; n = 400, p = 1000; SNR = 3"
    )
  )
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------
cl <- NULL
if (ncores > 1) {
  cl <- parallel::makeCluster(ncores,
    outfile = file.path(out_dir, "workers.log"))
  parallel::clusterCall(cl, function(root, code_files) {
    setwd(root)
    suppressPackageStartupMessages({
      library(MASS)
      library(sulnet)
      library(glmnet)
      library(uniLasso)
    })
    for (f in code_files) source(f)
    NULL
  }, getwd(), code_files)
  parallel::clusterExport(cl, c("run_rep", "rep_data", "draw_foldid",
    "conflict_defs", "CONFLICT_SNR"))
}

for (k in models) {
  dat <- if (k <= 4) model_data(k, max(reps)) else list(shared = list())
  todo <- Filter(function(b) !file.exists(rep_file(k, b)), reps)
  message(sprintf("model %d: %d of %d replicates to run", k, length(todo),
    length(reps)))
  tasks <- lapply(todo, function(b) {
    list(k = k, b = b, file = rep_file(k, b),
      data = if (k <= 4) dat$reps[[b]])
  })
  if (length(tasks)) {
    t0 <- Sys.time()
    errs <- if (is.null(cl)) {
      lapply(tasks, run_rep, shared = dat$shared)
    } else {
      parallel::parLapplyLB(cl, tasks, run_rep, shared = dat$shared)
    }
    errs <- unlist(errs)
    message(sprintf("model %d: finished in %.1f min, %d failed", k,
      as.numeric(difftime(Sys.time(), t0, units = "mins")), length(errs)))
    for (e in errs) message("  ", e)
  }

  ## summary of every saved replicate of model k
  files <- list.files(rep_dir, pattern = sprintf("^m%d_b[0-9]+\\.rds$", k),
    full.names = TRUE)
  if (!length(files)) next
  seeds <- as.integer(sub(".*_b([0-9]+)\\.rds$", "\\1", files))
  files <- files[order(seeds)]
  saved <- lapply(files, readRDS)
  raw <- simplify2array(lapply(saved, `[[`, "crit"))  # M x 6 x B
  raw <- aperm(raw, c(3, 1, 2))
  dimnames(raw) <- list(NULL, method_names_tuned, criteria_names)
  n_ok <- apply(!is.na(raw[, , "sse", drop = FALSE]), 2, sum)
  avg <- apply(raw, c(2, 3), mean, na.rm = TRUE)
  se <- apply(raw, c(2, 3), sd, na.rm = TRUE) / sqrt(n_ok)
  table <- noquote(matrix(sprintf("%.4f (%.4f)", avg, se),
    nrow = nrow(avg), dimnames = dimnames(avg)
  ))
  saveRDS(list(mean = avg, se = se, table = table, n = n_ok, raw = raw,
    reps = sort(seeds), seconds = vapply(saved, `[[`, 0, "seconds")),
    file.path(out_dir, sprintf("res_model%d.rds", k)))
  message("model ", k, " (", length(files), " replicates):")
  print(table)
  write_tables()
}
if (!is.null(cl)) parallel::stopCluster(cl)
