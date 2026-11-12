## Round 3: screen of sign-conflict designs. Usage, from the package root:
##   Rscript inst/simulation/binomial_sign_conflict/grid_tuning/round3.R \
##     [designs] [reps] [ncores] [out_dir]
## Defaults: 11:22, 1:30, 1, <TUNE_DIR>/results/round3. Paths: see sim_lib.R.
## designs are ids 11..22 (R1..R12 below). Each task stores baselines and two
## zeta_tau fits (standardize FALSE / TRUE) on a dense tau grid relative to that
## fit's zmax; per-cell CV (loss, Brier) and test criteria are kept.
## Ids 15, 16 and 19 became models 7, 10 and 11 of sign_conflict_simulation.R.
## Argon run (Nov 2026): 11:22, 1:30, two 40-core jobs (about 35 minutes).
## With ROUND3_DEFS_ONLY set, only the definitions are made (used by final_30reps.R).
defs_only <- nzchar(Sys.getenv("ROUND3_DEFS_ONLY"))
tune_dir <- Sys.getenv("TUNE_DIR", file.path(Sys.getenv("SIM_DIR", "inst/simulation"),
  "binomial_sign_conflict", "grid_tuning"))
source(file.path(tune_dir, "sim_lib.R"))

## active positions: the first length(beta) variables of block k
blk <- function(k, beta, B) (k - 1) * B + seq_along(beta)
rep_blocks <- function(beta, blocks, B) {
  list(active = unlist(lapply(blocks, blk, beta = beta, B = B)), beta = rep(beta, length(blocks)))
}
mix <- function(...) {
  l <- list(...)
  list(active = unlist(lapply(l, `[[`, "active")), beta = unlist(lapply(l, `[[`, "beta")))
}
sb <- function(n, a, design, rho, B = 10L, b) {
  simulate_binomial(n = n, p = 1000, beta_active = a$beta, active = a$active,
    design = design, rho = rho, block_size = B, n_test = 100, seed = b)
}
designs <- list(
  ## R1: model 9 at n = 200
  "11" = function(b) sb(200, rep_blocks(c(2, 2, -2), 1:3, 10), "block", 0.8, 10L, b),
  ## R2: model 9 with six conflict blocks
  "12" = function(b) sb(400, rep_blocks(c(2, 2, -2), 1:6, 10), "block", 0.8, 10L, b),
  ## R3: three conflict blocks + three conflict-free blocks
  "13" = function(b) sb(400, mix(rep_blocks(c(2, 2, -2), 1:3, 10), rep_blocks(c(1, 1, 1), 4:6, 10)),
    "block", 0.8, 10L, b),
  ## R4: high correlation, blocks of 20
  "14" = function(b) sb(200, rep_blocks(c(1.5, 1.5, -1.5), 1:3, 20), "block", 0.9, 20L, b),
  ## R5: moderate correlation, weak conflict
  "15" = function(b) sb(400, rep_blocks(c(1, 1, 1, -1), 1:3, 10), "block", 0.5, 10L, b),
  ## R6: AR(1) 0.8 with (1, 1, 1, -1.5) at three places
  "16" = function(b) sb(400, list(active = c(1:4, 101:104, 201:204), beta = rep(c(1, 1, 1, -1.5), 3)),
    "ar1", 0.8, 10L, b),
  ## R7: suppressor pairs (2, -1), rho 0.9
  "17" = function(b) sb(200, rep_blocks(c(2, -1), 1:5, 2), "block", 0.9, 2L, b),
  ## R8: three suppressor pairs + six positive signals whose pair partner is noise
  "18" = function(b) sb(200, mix(rep_blocks(c(2, -1), 1:3, 2), rep_blocks(1, 4:9, 2)), "block", 0.9, 2L, b),
  ## R9: conflict-free, highly correlated proxies, small n
  "19" = function(b) sb(100, rep_blocks(c(1, 1), 1:5, 20), "block", 0.9, 20L, b),
  ## R10: R9-type blocks + two conflict blocks
  "20" = function(b) sb(150, mix(rep_blocks(c(1, 1), 1:3, 20), rep_blocks(c(1.5, 1.5, -1.5), 4:5, 20)),
    "block", 0.9, 20L, b),
  ## R11: global equicorrelation 0.3, 8 positive and 2 negative (conflict) signals
  "21" = function(b) sb(200, list(active = 1:10, beta = c(rep(1, 8), -1, -1)), "equicor", 0.3, 10L, b),
  ## R12: many weak conflicts (true -1, univariate slope +1.4)
  "22" = function(b) sb(400, rep_blocks(c(1.5, 1.5, -1), 1:6, 10), "block", 0.8, 10L, b)
)

gen_data_orig <- gen_data
gen_data <- function(model, b, nfolds = 10, draw_folds = TRUE) {
  if (model <= 9) return(gen_data_orig(model, b, nfolds, draw_folds))
  dat <- designs[[as.character(model)]](b)
  X <- dat$train$x
  storage.mode(X) <- "double"
  y <- as.numeric(drop(dat$train$y))
  n <- nrow(X)
  foldid <- integer(n)
  if (draw_folds) foldid[order(y, runif(n))] <- rep_len(seq_len(nfolds), n)
  list(
    x = X, y = y, xtest = as.matrix(dat$test$x), gamma = as.numeric(dat$beta),
    gamma0 = dat$intercept, eta_new = dat$intercept + drop(dat$test$x %*% dat$beta),
    foldid = foldid
  )
}

TAU_REL <- c(1e3, 10^seq(1, -4, by = -0.25), 0)
zmax_of <- function(d, std) {
  suni2D(d$x, d$y, family = "binomial", standardize = std, nzeta = 3, ntau = 3)$zeta[1] / (1 + 1e-9)
}
configs <- list(
  list(name = "zt_F", method = "zeta_tau", opts = list(standardize = FALSE),
    grid = function(d, zmax) list(tau = zmax * TAU_REL)),
  list(name = "zt_T", method = "zeta_tau", opts = list(standardize = TRUE),
    grid = function(d, zmax) list(tau = zmax_of(d, TRUE) * TAU_REL))
)

if (!defs_only) {
models <- eval(parse(text = arg(1, "11:22")))
reps <- eval(parse(text = arg(2, "1:30")))
ncores <- as.integer(arg(3, "1"))
out_dir <- arg(4, file.path(tune_dir, "results", "round3"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(b = reps, model = models)
nn <- c("11" = 200, "12" = 400, "13" = 400, "14" = 200, "15" = 400, "16" = 400, "17" = 200,
  "18" = 200, "19" = 100, "20" = 150, "21" = 200, "22" = 400)
tasks <- tasks[order(-nn[as.character(tasks$model)], tasks$model, tasks$b), ]
message(nrow(tasks), " tasks on ", ncores, " cores")
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
}
