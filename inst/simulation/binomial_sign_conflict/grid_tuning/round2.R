## Round 2. Usage, from the package root:
##   Rscript inst/simulation/binomial_sign_conflict/grid_tuning/round2.R \
##     [set] [models] [reps] [ncores] [out_dir]
##   set "zeta":  first-coordinate grid variants (tau = {1e3, 0} * zmax), models 1-6
##   set "dense": the round-1 configurations, extended toward tau = 0 /
##                alpha = 0.5 / rate = 1, for the sign-conflict models 7-9 of sim_lib.R
## Defaults: zeta, 1:6 (zeta) or 7:9 (dense), 1:30, 1, <TUNE_DIR>/results/round2_<set>.
## Argon runs (Nov 2026): "zeta" 1:6 1:30 on 20 cores (5 min);
## "dense" 7:9 1:30 on 40 cores (1.8 hours). Paths: see sim_lib.R.
tune_dir <- Sys.getenv("TUNE_DIR", file.path(Sys.getenv("SIM_DIR", "inst/simulation"),
  "binomial_sign_conflict", "grid_tuning"))
source(file.path(tune_dir, "sim_lib.R"))
set <- arg(1, "zeta")
models <- eval(parse(text = arg(2, if (set == "zeta") "1:6" else "7:9")))
reps <- eval(parse(text = arg(3, "1:30")))
ncores <- as.integer(arg(4, "1"))
out_dir <- arg(5, file.path(tune_dir, "results", paste0("round2_", set)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ext_pen1 <- function(zmax, K = 149) 2 * zmax * (1 + 1e-9) * 0.01^((-K:99) / 99)
alpha_dense <- sort(unique(round(c(
  seq(0.5, 1 / 1.01, length.out = 20), 1 - 0.5 * 10^seq(0, -3, by = -0.125),
  0.5 + 0.5 * 10^seq(-1, -4, by = -0.5) # just above 0.5 (small tau)
), 12)))
rate_dense <- sort(unique(round(c(
  10^seq(0, 3, by = 0.125), 100^((0:19) / 19),
  1 + 2 * 10^seq(-1, -4, by = -0.5) # just above 1 (small tau)
), 12)))
ends <- function(d, zmax) list(tau = zmax * c(1e3, 0))

configs <- switch(set,
  zeta = list(
    list(name = "zt_ends", method = "zeta_tau", opts = list(standardize = FALSE), grid = ends),
    list(name = "zt_ends_lf001", method = "zeta_tau",
      opts = list(standardize = FALSE, lambda.factor = 0.001, nzeta = 150), grid = ends),
    list(name = "zt_ends_nz200", method = "zeta_tau",
      opts = list(standardize = FALSE, nzeta = 200), grid = ends),
    list(name = "zt_ends_lf03", method = "zeta_tau",
      opts = list(standardize = FALSE, lambda.factor = 0.03, nzeta = 77), grid = ends)
  ),
  dense = list(
    list(name = "zt_defT", method = "zeta_tau", opts = list(standardize = TRUE)),
    list(name = "la_defT", method = "lambda_alpha", opts = list(standardize = TRUE)),
    list(name = "sr_defT", method = "sum_rate", opts = list(standardize = TRUE)),
    list(name = "zt_defF", method = "zeta_tau", opts = list(standardize = FALSE)),
    list(name = "la_defF", method = "lambda_alpha", opts = list(standardize = FALSE)),
    list(name = "sr_defF", method = "sum_rate", opts = list(standardize = FALSE)),
    list(name = "zt_dense", method = "zeta_tau", opts = list(standardize = FALSE),
      grid = function(d, zmax) list(tau = zmax * c(1e3, 10^seq(1, -5, by = -0.125), 0))),
    list(name = "la_dense", method = "lambda_alpha", opts = list(standardize = FALSE),
      grid = function(d, zmax) list(lambda = ext_pen1(zmax), alpha = alpha_dense)),
    list(name = "sr_dense", method = "sum_rate", opts = list(standardize = FALSE),
      grid = function(d, zmax) list(sum = ext_pen1(zmax), rate = rate_dense))
  )
)

tasks <- expand.grid(b = reps, model = models)
tasks <- tasks[order(-(tasks$model %in% c(6, 8, 9)), -(tasks$model != 1), tasks$model, tasks$b), ]
message(nrow(tasks), " tasks on ", ncores, " cores; set ", set, ", ", length(configs), " configs")
t0 <- Sys.time()
run_one <- function(i) {
  k <- tasks$model[i]
  b <- tasks$b[i]
  st <- run_task(k, b, configs, file.path(out_dir, sprintf("m%d_b%03d.rds", k, b)))
  message(format(Sys.time()), " model ", k, " rep ", b, ": ", st)
  st
}
out <- if (ncores > 1) {
  parallel::mclapply(seq_len(nrow(tasks)), run_one, mc.cores = ncores, mc.preschedule = FALSE)
} else {
  lapply(seq_len(nrow(tasks)), run_one)
}
print(table(unlist(out)))
message("total time: ", format(Sys.time() - t0))
