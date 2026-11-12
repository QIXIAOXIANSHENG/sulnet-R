## 30-replicate run of the paper pipeline ("binomial simulation.R":
## run_simulation_binom, B = 1 per seed, s = "min") with the settings chosen in
## rounds 1-2 (the same as sign_conflict_simulation.R):
##   * cv.uniLasso(loo = FALSE)
##   * suni2D standardize = FALSE (all three parameterizations)
##   * zeta_tau: tau = zmax * c(1e3, 10^seq(0, -3, by = -0.25), 0), where zmax is
##     the largest zeta of the default grid; default zeta grid
##   * lambda_alpha, sum_rate: default grids (they contain alpha = 0.5 / rate = 1)
## Model ids: 1-9 as in sim_lib.R (1-6 = "binomial simulation.R"), 11-22 as in
## round3.R. Usage, from the package root:
##   Rscript inst/simulation/binomial_sign_conflict/grid_tuning/final_30reps.R \
##     [models] [reps] [ncores] [out_dir]
## Defaults: 1:6, 1:30, 1, <TUNE_DIR>/results/final_30reps. Paths: see sim_lib.R.
## Argon runs (Nov 2026): 1:9 and c(11,12,13,15,16,19,22), seeds 1:30, 40 cores;
## the models 1-6 results are in results_30reps/.
tune_dir <- Sys.getenv("TUNE_DIR", file.path(Sys.getenv("SIM_DIR", "inst/simulation"),
  "binomial_sign_conflict", "grid_tuning"))
source(file.path(tune_dir, "sim_lib.R")) # compare_binom.R, setup, gen_data(), arg()
models <- eval(parse(text = arg(1, "1:6")))
reps <- eval(parse(text = arg(2, "1:30")))
ncores <- as.integer(arg(3, "1"))
out_dir <- arg(4, file.path(tune_dir, "results", "final_30reps"))
if (any(models > 9)) { # round-3 designs 11..22; must come before the overrides below
  Sys.setenv(ROUND3_DEFS_ONLY = "1")
  source(file.path(tune_dir, "round3.R"))
}
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
TAU_REL <- c(1e3, 10^seq(0, -3, by = -0.25), 0)

## compare_binom.R's fit_methods_binom() with two changes: uniLasso without
## LOO, and tau_rel, a zeta_tau tau grid relative to zmax. Other arguments in
## ... go to cv.suni2D (and to the zmax fit).
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

tasks <- expand.grid(b = reps, model = models)
tasks <- tasks[order(-(tasks$model %in% c(6, 8, 9)), tasks$model, tasks$b), ]
message(nrow(tasks), " tasks on ", ncores, " cores")
t0 <- Sys.time()
run_one <- function(i) {
  k <- tasks$model[i]
  b <- tasks$b[i]
  f <- file.path(rep_dir, sprintf("m%d_b%03d.rds", k, b))
  if (file.exists(f)) return("skip")
  msgs <- character()
  t1 <- proc.time()[["elapsed"]]
  raw <- withCallingHandlers(
    tryCatch({
      d <- gen_data(k, b, draw_folds = FALSE)
      run_simulation_binom(d$x, d$y, d$gamma, d$xtest,
        gamma0 = d$gamma0, B = 1,
        standardize = FALSE, tau_rel = TAU_REL
      )$raw[1, , ]
    }, error = function(e) structure(conditionMessage(e), class = "sim_error")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  ok <- !inherits(raw, "sim_error")
  saveRDS(list(
    model = k, b = b, raw = if (ok) raw, error = if (!ok) unclass(raw),
    messages = msgs, time = proc.time()[["elapsed"]] - t1
  ), f)
  message(format(Sys.time()), " model ", k, " rep ", b, ": ", if (ok) "ok" else "error")
  if (ok) "ok" else "error"
}
out <- if (ncores > 1) {
  parallel::mclapply(seq_len(nrow(tasks)), run_one, mc.cores = ncores, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(tasks)), run_one)
}
print(table(unlist(out)))
message("total time: ", format(Sys.time() - t0))

## combine into the summarize_sim() format of inst/simulation/binomial/res_model<k>.rds
for (k in models) {
  fs <- list.files(rep_dir, sprintf("^m%d_b\\d+\\.rds$", k), full.names = TRUE)
  rr <- lapply(fs, readRDS)
  bad <- vapply(rr, function(r) is.null(r$raw), logical(1))
  if (any(bad)) message("model ", k, " errors: ", paste(vapply(rr[bad], `[[`, "", "error"), collapse = " | "))
  rr <- rr[!bad]
  if (!length(rr)) next
  res_k <- summarize_sim(lapply(rr, `[[`, "raw"))
  res_k$n <- length(rr)
  res_k$seeds <- vapply(rr, `[[`, numeric(1), "b")
  res_k$messages <- unique(unlist(lapply(rr, `[[`, "messages")))
  res_k$time <- vapply(rr, `[[`, numeric(1), "time")
  saveRDS(res_k, file.path(out_dir, sprintf("res_model%d.rds", k)))
  cat("\n==== model", k, "(", length(rr), "reps )\n")
  print(res_k$table)
}
