## Linear models 1-5 of inst/simulation/"linear simulation.R", rerun with the
## fixed glmnet settings and the honest-CV adaptive lassos of compare_linear.R
## (2026-11-12). See README.md in this folder.
##
## Usage, from the package root (the folder that contains inst/):
##   Rscript inst/simulation/linear_fixed/linear_fixed_simulation.R \
##     [models] [reps] [ncores] [out_dir]
## Defaults: models = 1:5, reps = 1:200, ncores = 1,
##   out_dir = inst/simulation/linear_fixed
## Example, models 1 and 2 with replicates 1-10 on 8 cores:
##   Rscript <this file> 1:2 1:10 8
##
## The data, the replicate responses and the CV folds are those of
## "linear simulation.R": the random draws are made in the same order and with
## the same seeds as there (fit_methods() draws no random numbers), and are made
## here, before any fitting, so that the replicates can run in parallel.
## ncores > 1 uses a PSOCK cluster, which works on Windows.
##
## Each replicate is saved as <out_dir>/reps/m<k>_b<bbb>.rds. A replicate whose
## file already exists is skipped, so an interrupted run resumes where it
## stopped when rerun with the same arguments. At the end, for each model,
## <out_dir>/res_model<k>.rds is rebuilt from every saved replicate (elements
## mean, se, table, n, raw, as returned by run_simulation(), plus seconds), and
## <out_dir>/linear_fixed_tables.pdf from every res_model<k>.rds.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) {
  if (length(args) >= i && nzchar(args[i])) args[i] else default
}
models <- eval(parse(text = arg(1, "1:5")))
reps <- eval(parse(text = arg(2, "1:200")))
ncores <- as.integer(arg(3, "1"))
out_dir <- arg(4, "inst/simulation/linear_fixed")
sim_dir <- "inst/simulation"

suppressPackageStartupMessages({
  library(MASS)
  library(sulnet)
  library(glmnet)
  library(uniLasso)
})
source(file.path(sim_dir, "compare_linear.R"))
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Data. model_data(k, B) returns
##   shared  data common to all replicates of model k
##   reps    one list per replicate b = 1..B with the rest of
##           x, y, gamma, sigma, xtest, foldid
## made with the draws of "linear simulation.R" in their order there.
## ---------------------------------------------------------------------------
draw_foldid <- function(n, nfolds = 10L) {
  sample(rep(seq_len(nfolds), length.out = n))
}

## Models 1-3: fixed x; replicate 1 is the y drawn with x, replicates 2..B
## redraw y = x gamma + sigma eps, as in run_simulation().
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
  if (k == 5) {
    ## uniLasso home court: new x, coefficients and y for seed b
    n <- 100; ntest <- 100; p <- 30; rho <- 0.8; sparsity <- 0.2
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    reps <- lapply(seq_len(B), function(b) {
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
      list(x = X, xtest = xtest, y = y, gamma = gamma_effective,
        sigma = sigma, foldid = draw_foldid(n))
    })
    return(list(shared = list(), reps = reps))
  }
  stop("unknown model ", k)
}

## ---------------------------------------------------------------------------
## One replicate: fit, score, save. Returns NULL, or the error message.
## ---------------------------------------------------------------------------
rep_file <- function(k, b) file.path(rep_dir, sprintf("m%d_b%03d.rds", k, b))

run_rep <- function(task, shared) {
  d <- c(task$data, shared)
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    fit <- withCallingHandlers(
      fit_methods(d$x, d$y, d$xtest, d$foldid, s = "min"),
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
## Run
## ---------------------------------------------------------------------------
cl <- NULL
if (ncores > 1) {
  cl <- parallel::makeCluster(ncores,
    outfile = file.path(out_dir, "workers.log"))
  parallel::clusterCall(cl, function(root, sim_dir) {
    setwd(root)
    suppressPackageStartupMessages({
      library(sulnet)
      library(glmnet)
      library(uniLasso)
    })
    source(file.path(sim_dir, "compare_linear.R"))
    NULL
  }, getwd(), sim_dir)
  parallel::clusterExport(cl, "run_rep")
}

for (k in models) {
  dat <- model_data(k, max(reps))
  todo <- Filter(function(b) !file.exists(rep_file(k, b)), reps)
  message(sprintf("model %d: %d of %d replicates to run", k, length(todo),
    length(reps)))
  tasks <- lapply(todo, function(b) {
    list(k = k, b = b, file = rep_file(k, b), data = dat$reps[[b]])
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
  dimnames(raw) <- list(NULL, method_names, criteria_names)
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
}
if (!is.null(cl)) parallel::stopCluster(cl)

## ---------------------------------------------------------------------------
## Tables of every model saved in out_dir
## ---------------------------------------------------------------------------
source(file.path(sim_dir, "table_simulation.R"))
write_sim_tables(read_sim_results(out_dir),
  file = file.path(out_dir, "linear_fixed_tables.pdf"),
  main = "Linear simulations, fixed glmnet paths and honest-CV adaptive lassos",
  captions = c(
    "Model 1" = "Huang, AR(1) rho = 0.5; n = 100, p = 200, 15 nonzero; sigma = 1.5",
    "Model 2" = "Huang, AR(1) rho = 0.95; n = 100, p = 200, 15 nonzero; sigma = 1.5",
    "Model 3" = "Huang, 3 groups of 5 near-collinear signals; n = 100, p = 200; sigma = 1.5",
    "Model 4" = "uniLasso, equicorrelated rho = 0.5; n = 300, p = 1000, 30 nonzero ~ U(0, 1); sigma = 0.1",
    "Model 5" = "uniLasso home court, AR(1) rho = 0.8; n = 100, p = 30, 20% nonzero; SNR = 1"
  )
)
