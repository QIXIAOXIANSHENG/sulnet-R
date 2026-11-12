## Models 1-6 of inst/simulation/"binomial simulation.R", rerun with the
## settings chosen in the grid-tuning study (see README.md in this folder and
## inst/simulation/binomial_sign_conflict/README.md).
##
## Same data (same seeds) and pipeline as "binomial simulation.R"
## (run_simulation_binom() with B = 1 per seed, s = "min", 10 class-stratified
## folds shared by all methods), with these changes:
##   * uniLasso: cv.uniLasso(..., loo = FALSE); the adalasso_uni weights also
##     come from uniFit(..., loo = FALSE), as in compare_binom.R
##   * suni2D (all three parameterizations): standardize = FALSE
##   * suni2D zeta_tau: tau = zmax * c(1e3, 10^seq(0, -3, by = -0.25), 0), where
##     zmax is the largest zeta of the default zeta grid; lambda_alpha and
##     sum_rate keep their default grids (which contain alpha = 0.5 / rate = 1)
##   * every model is scored against its true intercept: gamma0 = dat$intercept
##     for models 2-6 (about -2.19 in model 5, 0 up to 1e-9 elsewhere) and
##     gamma0 = 0 for model 1
##
## Usage, from the package root (the folder that contains inst/):
##   Rscript inst/simulation/binomial_tuned/binomial_tuned_simulation.R \
##     [models] [reps] [ncores] [sim_dir] [out_dir]
## Defaults: models = 1:6, reps = 1:200, ncores = 1,
##   sim_dir = inst/simulation (holds compare_binom.R,
##             binomial_simulation_setup.R),
##   out_dir = <sim_dir>/binomial_tuned/rerun
## The default out_dir is a subfolder so that a local run never overwrites the
## 200-replicate results res_model1-6.rds shipped in this folder.
## Example, model 2 seeds 1-10 only: Rscript <this file> 2 1:10
## The Argon run used: Rscript scripts/binomial_tuned_simulation.R <models> \
##   1:200 $NSLOTS scripts results/tuned   (see argon/runjob_tuned.txt)
##
## Each replicate is saved as <out_dir>/reps/m<k>_b<bbb>.rds. A replicate whose
## file already exists is skipped, so an interrupted run resumes where it
## stopped when rerun with the same arguments. At the end,
## <out_dir>/res_model<k>.rds is rebuilt from every saved replicate of model k
## (summarize_sim() format plus raw, n, seeds, ...).
## ncores > 1 uses parallel::mclapply(), which cannot fork on Windows, so
## ncores is set to 1 there. Minutes per replicate on one core (Argon):
## model 1 about 0.7, models 2-5 about 3-4, model 6 about 10.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
models <- eval(parse(text = arg(1, "1:6")))
reps <- eval(parse(text = arg(2, "1:200")))
ncores <- as.integer(arg(3, "1"))
sim_dir <- arg(4, "inst/simulation")
out_dir <- arg(5, file.path(sim_dir, "binomial_tuned", "rerun"))
if (ncores > 1 && .Platform$OS.type == "windows") {
  message("mclapply cannot fork on Windows; using 1 core")
  ncores <- 1L
}

suppressPackageStartupMessages({
  library(sulnet)
  library(glmnet)
  library(uniLasso)
  library(MASS)
})
source(file.path(sim_dir, "compare_binom.R"))
source(file.path(sim_dir, "binomial_simulation_setup.R"))
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## Models: the data of "binomial simulation.R" for seed b.
## ---------------------------------------------------------------------------
model_settings <- list(
  "1" = list(design = "uniLasso case: AR(1) rho = 0.8, class mean shift 0.5 in x1-x20 (training only)",
    n = 100, p = 300, n_test = 30, nonzero = "x1-x30 = 1", intercept = 0),
  "2" = list(design = "independent", n = 200, p = 1000, n_test = 100,
    nonzero = "x1-x10 = (1, -1) x 5", prevalence = 0.5),
  "3" = list(design = "AR(1) rho = 0.5", n = 200, p = 1000, n_test = 100,
    nonzero = "x1-x10 = (1, -1) x 5", prevalence = 0.5),
  "4" = list(design = "blocks of 20, rho = 0.5", n = 200, p = 1000, n_test = 100,
    nonzero = "x1-x10 = (1, -1) x 5", prevalence = 0.5),
  "5" = list(design = "imbalanced: AR(1) rho = 0.5", n = 200, p = 1000, n_test = 100,
    nonzero = "x1-x10 = (1, -1) x 5", prevalence = 0.2),
  "6" = list(design = "Breheny-Huang (2011), independent", n = 400, p = 1000,
    n_test = 100, nonzero = "x1-x4 = (1, 1, -1, -1)", intercept = 0)
)
## Returns x, y, gamma, xtest, gamma0. The RNG is left where
## "binomial simulation.R" leaves it, so run_simulation_binom() draws the same
## folds.
model_data <- function(k, b) {
  if (k == 1) {
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
    return(list(x = X, y = y, gamma = gamma, xtest = xtest, gamma0 = 0))
  }
  dat <- switch(as.character(k),
    "2" = simulation_setting("independent", n = 200, p = 1000, n_test = 100, seed = b),
    "3" = simulation_setting("ar1", n = 200, p = 1000, n_test = 100, seed = b),
    "4" = simulation_setting("block", n = 200, p = 1000, n_test = 100, seed = b),
    "5" = simulation_setting("imbalanced", n = 200, p = 1000, n_test = 100, seed = b),
    "6" = bh2011(n = 400, p = 1000, amplitude = 1, n_test = 100, seed = b),
    stop("unknown model ", k)
  )
  list(x = dat$train$x, y = dat$train$y, gamma = dat$beta, xtest = dat$test$x,
    gamma0 = dat$intercept)
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
## from the RNG state left by the data step, as in "binomial simulation.R".
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
      run_simulation_binom(dat$x, dat$y, dat$gamma, dat$xtest,
        gamma0 = dat$gamma0, B = 1,
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
## slowest model first so the cores stay busy at the end
tasks <- tasks[order(-(tasks$model == 6), tasks$model, tasks$b), ]
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
  res_k$settings <- c(model_settings[[as.character(k)]],
    list(tau_rel = TAU_REL, standardize = FALSE, uniLasso_loo = FALSE,
      pred.loss = "loss", s = "cv.min / lambda.min"))
  res_k$messages <- unique(unlist(lapply(rr, `[[`, "messages")))
  res_k$time <- vapply(rr, `[[`, numeric(1), "time")
  saveRDS(res_k, file.path(out_dir, sprintf("res_model%d.rds", k)))
  cat("\n==== model", k, "(", length(rr), "replicates )\n")
  print(res_k$table)
}
