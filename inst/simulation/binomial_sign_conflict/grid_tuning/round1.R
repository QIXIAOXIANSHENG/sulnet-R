## Round 1: dense second-coordinate grids, standardize = FALSE, plus the
## default grids (standardize TRUE and FALSE) for reference, models 1-6.
## Usage, from the package root:
##   Rscript inst/simulation/binomial_sign_conflict/grid_tuning/round1.R \
##     [models] [reps] [ncores] [out_dir]
## Defaults: 1:6, 1:30, 1, <TUNE_DIR>/results/round1. Paths: see sim_lib.R.
## Argon run (Nov 2026): models 1:6, reps 1:30, 40 cores (about 2.5 hours).
## Finished tasks are skipped on a rerun (one file per model and replicate).
tune_dir <- Sys.getenv("TUNE_DIR", file.path(Sys.getenv("SIM_DIR", "inst/simulation"),
  "binomial_sign_conflict", "grid_tuning"))
source(file.path(tune_dir, "sim_lib.R"))
models <- eval(parse(text = arg(1, "1:6")))
reps <- eval(parse(text = arg(2, "1:30")))
ncores <- as.integer(arg(3, "1"))
out_dir <- arg(4, file.path(tune_dir, "results", "round1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## First-coordinate grid for lambda_alpha / sum_rate: the default
## 2 * zmax * 0.01^((l - 1) / 99), l = 1..100, extended upward on the same
## log spacing by K values, so that alpha up to 1 - 1/(2 * 0.01^(-K/99))
## (0.9995 for K = 149) still starts from the null model.
ext_pen1 <- function(zmax, K = 149) 2 * zmax * (1 + 1e-9) * 0.01^((-K:99) / 99)

alpha_dense <- sort(unique(round(c(
  seq(0.5, 1 / 1.01, length.out = 20), # default-like linear grid
  1 - 0.5 * 10^seq(0, -3, by = -0.125) # log-spaced in 1 - alpha
), 12)))
rate_dense <- sort(unique(round(c(
  10^seq(0, 3, by = 0.125),
  100^((0:19) / 19) # default-like log grid 1..100
), 12)))

configs <- list(
  list(name = "zt_defT", method = "zeta_tau", opts = list(standardize = TRUE)),
  list(name = "la_defT", method = "lambda_alpha", opts = list(standardize = TRUE)),
  list(name = "sr_defT", method = "sum_rate", opts = list(standardize = TRUE)),
  list(name = "zt_defF", method = "zeta_tau", opts = list(standardize = FALSE)),
  list(name = "la_defF", method = "lambda_alpha", opts = list(standardize = FALSE)),
  list(name = "sr_defF", method = "sum_rate", opts = list(standardize = FALSE)),
  list(
    name = "zt_dense", method = "zeta_tau", opts = list(standardize = FALSE),
    grid = function(d, zmax) list(tau = zmax * c(1e3, 10^seq(1, -3, by = -0.125), 0))
  ),
  list(
    name = "la_dense", method = "lambda_alpha", opts = list(standardize = FALSE),
    grid = function(d, zmax) list(lambda = ext_pen1(zmax), alpha = alpha_dense)
  ),
  list(
    name = "sr_dense", method = "sum_rate", opts = list(standardize = FALSE),
    grid = function(d, zmax) list(sum = ext_pen1(zmax), rate = rate_dense)
  )
)

tasks <- expand.grid(b = reps, model = models)
## heavier models first for load balancing
tasks <- tasks[order(-(tasks$model == 6), -(tasks$model != 1), tasks$model, tasks$b), ]
message(nrow(tasks), " tasks on ", ncores, " cores; ", length(configs), " suni configs")
t0 <- Sys.time()
run_one <- function(i) {
  k <- tasks$model[i]
  b <- tasks$b[i]
  st <- run_task(k, b, configs, file.path(out_dir, sprintf("m%d_b%03d.rds", k, b)))
  message(format(Sys.time()), " model ", k, " rep ", b, ": ", st)
  st
}
if (ncores > 1) {
  out <- parallel::mclapply(seq_len(nrow(tasks)), run_one,
    mc.cores = ncores, mc.preschedule = FALSE
  )
} else {
  out <- lapply(seq_len(nrow(tasks)), run_one)
}
print(table(unlist(out)))
message("total time: ", format(Sys.time() - t0))
