## Binomial sign-conflict simulations: models 7-11, continuing models 1-6 of
## inst/simulation/"binomial simulation.R". See README.md in this folder for
## the settings and the reasoning behind them.
##
## Same pipeline as "binomial simulation.R" (run_simulation_binom() with B = 1
## per seed, s = "min", 10 class-stratified folds shared by all methods), with
## three changes:
##   * uniLasso: cv.uniLasso(..., loo = FALSE)
##   * suni2D (all three parameterizations): standardize = FALSE
##   * suni2D zeta_tau: tau = zmax * c(1e3, 10^seq(0, -3, by = -0.25), 0), where
##     zmax is the largest zeta of the default zeta grid; lambda_alpha and
##     sum_rate keep their default grids (which contain alpha = 0.5 / rate = 1)
##
## Usage, from the package root (the folder that contains inst/):
##   Rscript inst/simulation/binomial_sign_conflict/sign_conflict_simulation.R \
##     [models] [reps] [ncores] [sim_dir] [out_dir]
## Defaults: models = 7:11, reps = 1:200, ncores = 1,
##   sim_dir = inst/simulation (holds compare_binom.R,
##             binomial_simulation_setup.R),
##   out_dir = <sim_dir>/binomial_sign_conflict/rerun
## The default out_dir is a subfolder so that a local run never overwrites the
## 200-replicate results res_model7-11.rds shipped in this folder.
## Example, model 9 seeds 1-10 only: Rscript <this file> 9 1:10
## The Argon run used: Rscript scripts/sign_conflict_simulation.R <k> 1:200 \
##   $NSLOTS scripts results/sign_conflict   (see argon/runjob_sign_conflict.txt)
##
## Each replicate is saved as <out_dir>/reps/m<k>_b<bbb>.rds. A replicate whose
## file already exists is skipped, so an interrupted run (for example by an
## Argon maintenance) resumes where it stopped when rerun with the same
## arguments. At the end, <out_dir>/res_model<k>.rds is rebuilt from every saved
## replicate of model k (summarize_sim() format plus raw, n, seeds, ...).
## ncores > 1 uses parallel::mclapply(), which cannot fork on Windows, so
## ncores is set to 1 there. One replicate of an n = 400 model takes about
## 5-10 minutes on one core; the n = 100 model about 1 minute.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
models <- eval(parse(text = arg(1, "7:11")))
reps <- eval(parse(text = arg(2, "1:200")))
ncores <- as.integer(arg(3, "1"))
sim_dir <- arg(4, "inst/simulation")
out_dir <- arg(5, file.path(sim_dir, "binomial_sign_conflict", "rerun"))
if (ncores > 1 && .Platform$OS.type == "windows") {
  message("mclapply cannot fork on Windows; using 1 core")
  ncores <- 1L
}

suppressPackageStartupMessages({
  library(sulnet)
  library(glmnet)
  library(uniLasso)
})
source(file.path(sim_dir, "compare_binom.R"))
source(file.path(sim_dir, "binomial_simulation_setup.R"))
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Models. p = 1000, 100 test points, intercept 0, seed b for replicate b.
## A "conflict" variable has a negative coefficient but a positive univariate
## slope: in an equicorrelated block, cov(x_j, eta) = beta_j + rho * (sum of
## the other coefficients in the block) > 0. uniLasso (tau = Inf) forces such a
## coefficient to be >= 0.
## ---------------------------------------------------------------------------
model_defs <- list(
  ## blocks of 10, rho 0.5; (1, 1, 1, -1) in blocks 1-3: conflict slope +0.5
  "7" = list(n = 400, design = "block", rho = 0.5, block_size = 10L,
    beta = rep(c(1, 1, 1, -1), 3), active = c(1:4, 11:14, 21:24)),
  ## blocks of 10, rho 0.7; (1, 1, 1, -1.5) in blocks 1-3: conflict slope +0.6
  "8" = list(n = 400, design = "block", rho = 0.7, block_size = 10L,
    beta = rep(c(1, 1, 1, -1.5), 3), active = c(1:4, 11:14, 21:24)),
  ## blocks of 10, rho 0.8; (2, 2, -2) in blocks 1-3: conflict slope +1.2
  "9" = list(n = 400, design = "block", rho = 0.8, block_size = 10L,
    beta = rep(c(2, 2, -2), 3), active = c(1:3, 11:13, 21:23)),
  ## AR(1) rho 0.8; (1, 1, 1, -1.5) at x1-x4, x101-x104, x201-x204: slope +0.45
  "10" = list(n = 400, design = "ar1", rho = 0.8, block_size = 10L,
    beta = rep(c(1, 1, 1, -1.5), 3), active = c(1:4, 101:104, 201:204)),
  ## control without conflicts: blocks of 20, rho 0.9; (1, 1) in blocks 1-5
  "11" = list(n = 100, design = "block", rho = 0.9, block_size = 20L,
    beta = rep(c(1, 1), 5), active = c(1:2, 21:22, 41:42, 61:62, 81:82))
)
model_data <- function(k, b) {
  m <- model_defs[[as.character(k)]]
  simulate_binomial(n = m$n, p = 1000, beta_active = m$beta, active = m$active,
    design = m$design, rho = m$rho, block_size = m$block_size,
    n_test = 100, seed = b)
}

## ---------------------------------------------------------------------------
## compare_binom.R's fit_methods_binom() with two changes: uniLasso without
## LOO, and tau_rel, a zeta_tau tau grid relative to zmax. Other arguments in
## ... (here standardize = FALSE) go to cv.suni2D and to the zmax fit.
## ---------------------------------------------------------------------------
TAU_REL <- c(1e3, 10^seq(0, -3, by = -0.25), 0)
fit_methods_binom <- function(x, y, newx, foldid, s = c("min", "1se"),
                              tau_rel = NULL, ...) {
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
    extra <- list()
    if (method == "zeta_tau" && !is.null(tau_rel)) {
      ## zeta[1] is zmax * (1 + 1e-9) in the default zeta grid
      f3 <- suni2D(x, y, family = "binomial", nzeta = 3, ntau = 3, ...)
      extra$tau <- f3$zeta[1] / (1 + 1e-9) * tau_rel
    }
    cvfit <- do.call(cv.suni2D, c(
      list(x = x, y = y, foldid = foldid, method = method, family = "binomial"),
      extra, list(...)
    ))
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
    cv.uniLasso(x, y, family = "binomial", foldid = foldid, loo = FALSE)
  ))
  fits$lasso <- safely("lasso", glmnet_out(
    cv.glmnet(x, y, family = "binomial", foldid = foldid)
  ))
  fits$adalasso_uni <- safely("adalasso_uni", adalasso_out(
    uniFit(x, y, loo = FALSE, family = "binomial")$beta
  ))
  fits$adalasso_lasso <- safely("adalasso_lasso", adalasso_out(fits$lasso$coef))
  fits <- fits[method_names_binom]
  list(
    coef = do.call(cbind, lapply(fits, `[[`, "coef")),
    link = do.call(cbind, lapply(fits, `[[`, "link"))
  )
}

## ---------------------------------------------------------------------------
## One replicate: data from seed b, then run_simulation_binom() draws the folds
## from the RNG state left by simulate_binomial(), as in "binomial simulation.R".
## ---------------------------------------------------------------------------
rep_file <- function(k, b) file.path(rep_dir, sprintf("m%d_b%03d.rds", k, b))
run_one <- function(k, b) {
  f <- rep_file(k, b)
  if (file.exists(f)) return("skip")
  msgs <- character()
  t1 <- proc.time()[["elapsed"]]
  raw <- withCallingHandlers(
    tryCatch({
      dat <- model_data(k, b)
      run_simulation_binom(dat$train$x, dat$train$y, dat$beta, dat$test$x,
        gamma0 = dat$intercept, B = 1,
        standardize = FALSE, tau_rel = TAU_REL
      )$raw[1, , ]
    }, error = function(e) structure(conditionMessage(e), class = "sim_error")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  ok <- !inherits(raw, "sim_error")
  if (!ok) {
    ## errors are not saved, so a rerun tries the replicate again
    message("model ", k, " rep ", b, " error: ", unclass(raw))
    return("error")
  }
  ## write to a temporary name first so that a killed job never leaves a
  ## truncated file that a rerun would skip
  tmp <- paste0(f, ".tmp")
  saveRDS(list(model = k, b = b, raw = raw, messages = msgs,
    time = proc.time()[["elapsed"]] - t1), tmp)
  file.rename(tmp, f)
  message(format(Sys.time()), " model ", k, " rep ", b, ": ok")
  "ok"
}

tasks <- expand.grid(b = reps, model = models)
todo <- !file.exists(rep_file(tasks$model, tasks$b))
message(nrow(tasks), " replicates requested, ", sum(!todo), " already saved, ",
  sum(todo), " to run on ", ncores, " core(s)")
tasks <- tasks[todo, , drop = FALSE]
t0 <- Sys.time()
if (nrow(tasks)) {
  one <- function(i) run_one(tasks$model[i], tasks$b[i])
  out <- if (ncores > 1) {
    parallel::mclapply(seq_len(nrow(tasks)), one, mc.cores = ncores, mc.preschedule = FALSE)
  } else {
    lapply(seq_len(nrow(tasks)), one)
  }
  print(table(unlist(out)))
  message("run time: ", format(Sys.time() - t0))
}

## ---------------------------------------------------------------------------
## Summaries: <out_dir>/res_model<k>.rds from all saved replicates of model k.
## table / mean / se as summarize_sim(); raw is replicates x methods x criteria.
## ---------------------------------------------------------------------------
for (k in models) {
  fs <- list.files(rep_dir, sprintf("^m%d_b\\d+\\.rds$", k), full.names = TRUE)
  if (!length(fs)) next
  rr <- lapply(fs, readRDS)
  seeds <- vapply(rr, `[[`, numeric(1), "b")
  rr <- rr[order(seeds)]
  seeds <- sort(seeds)
  raw_list <- lapply(rr, `[[`, "raw")
  res_k <- summarize_sim(raw_list)
  res_k$raw <- simplify2array(raw_list) # methods x criteria x replicates
  res_k$raw <- aperm(res_k$raw, c(3, 1, 2)) # replicates x methods x criteria
  dimnames(res_k$raw)[[1]] <- paste0("seed", seeds)
  res_k$n <- length(rr)
  res_k$seeds <- seeds
  res_k$settings <- c(model_defs[[as.character(k)]], p = 1000, n_test = 100)
  res_k$messages <- unique(unlist(lapply(rr, `[[`, "messages")))
  res_k$time <- vapply(rr, `[[`, numeric(1), "time")
  saveRDS(res_k, file.path(out_dir, sprintf("res_model%d.rds", k)))
  cat("\n==== model", k, "(", length(rr), "replicates )\n")
  print(res_k$table)
}
