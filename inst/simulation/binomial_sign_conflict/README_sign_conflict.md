# Binomial sign-conflict simulations (models 7-11)

This folder holds the binomial simulation models 7-11, which continue models 1-6 of
`inst/simulation/binomial simulation.R`, together with everything used to choose them:

- the scripts;
- the 200-replicate results;
- the Argon job files;
- the grid-tuning study behind the settings.

All runs were done on the Argon cluster (University of Iowa) in November 2026, with sulnet
built from GitHub `QIXIAOXIANSHENG/sulnet-R` at commit `b11feb6` (sulnet 0.4.0). The other
software was R 4.2.2, glmnet 5.2 and uniLasso 2.11.

## Contents

| file | what it is |
|---|---|
| `sign_conflict_simulation.R` | Runs models 7-11 (simulation, fitting, summaries). Resumes interrupted runs. |
| `make_tables.R` | Writes `binomial_sign_conflict_tables.pdf` from `res_model7-11.rds`. |
| `res_model7.rds` ... `res_model11.rds` | 200-replicate results (seeds 1-200), described below. |
| `binomial_sign_conflict_tables.pdf` | The summary tables of models 7-11. |
| `argon/runjob_sign_conflict.txt` | The SGE job used for the 200-replicate runs. |
| `argon/runjob_grid_tuning.txt` | An SGE job for the grid-tuning scripts, listing the commands used. |
| `grid_tuning/` | Scripts of the tuning study (rounds 1-3), runnable locally. |
| `grid_tuning/results_30reps/` | Models 1-6 rerun with the new settings, 30 replicates. Tables are in the pdf there. |

## Common settings (all of models 7-11)

**Data.** Each replicate b = 1, ..., 200 draws one data set with
`simulate_binomial(..., seed = b)` from `binomial_simulation_setup.R`:

- n training observations, as listed in the models table below;
- 100 test observations;
- p = 1000 Gaussian covariates with unit variances;
- intercept 0;
- y ~ Bernoulli(plogis(x'beta)), which gives a prevalence of about 0.5.

**Pipeline.** The pipeline is the same as for models 1-6. `run_simulation_binom()` from
`compare_binom.R` is called with `B = 1` and `s = "min"`. All methods share one assignment of
the 10 class-stratified CV folds. The folds are drawn right after the data, from the same seed.

**Methods.** The 7 methods are the rows of the tables; all are evaluated at the CV minimum:

| method | fit |
|---|---|
| `suni2D_zeta_tau` | `cv.suni2D(method = "zeta_tau", family = "binomial", standardize = FALSE)`; default zeta grid (100 values; lambda.factor 0.01 because n < p); **tau = zmax x {1e3, 10^0, 10^-0.25, ..., 10^-3, 0}** (15 values, zmax = largest zeta of the default grid); CV loss = deviance (the default `pred.loss = "loss"`) |
| `suni2D_lambda_alpha` | same with `method = "lambda_alpha"`, default grid (100 lambda x 20 alpha in [0.5, 0.99]) |
| `suni2D_sum_rate` | same with `method = "sum_rate"`, default grid (100 sum x 20 rate in [1, 100]) |
| `uniLasso` | `cv.uniLasso(family = "binomial", loo = FALSE)` |
| `lasso` | `cv.glmnet(family = "binomial")` |
| `adalasso_uni` | adaptive lasso, weights 1/abs(univariate slopes) from `uniFit(loo = FALSE)`, `standardize = FALSE` |
| `adalasso_lasso` | adaptive lasso, weights 1/abs(lasso coefficients) |

These settings differ from the original `binomial simulation.R` in three places:

1. `cv.uniLasso` uses `loo = FALSE`.
2. suni2D uses `standardize = FALSE`.
3. The zeta_tau tau grid contains 0 and a very large value.

The reasons are given under "Why these settings and models" below. The lambda_alpha and
sum_rate default grids already contain alpha = 0.5 and rate = 1, which correspond to tau = 0.

**Criteria.** The 7 criteria come from `eval_criteria_binom()` in `compare_binom.R`:

| criterion | definition |
|---|---|
| `sse` | squared error of the slope estimates |
| `rpe_kl` | Kullback-Leibler risk of the test probabilities, relative to the entropy (main criterion) |
| `rpe_class` | excess test misclassification, relative to the Bayes error |
| `support` | number of selected variables |
| `sign` | proportion of the p coordinates whose estimated sign equals the true sign (0 counts as a sign) |
| `tpr`, `fpr` | true- and false-positive selection rates |

## Models

A *sign conflict* is a variable whose true coefficient is negative but whose univariate slope
is positive. In an equicorrelated block,
cov(x_j, eta) = beta_j + rho * (sum of the other coefficients in the block). uniLasso
(tau = infinity) forces such a coefficient to be >= 0. suni2D can relax that constraint
through tau.

| model | design | n | active variables and coefficients | conflict slope cov(x_j, eta) | sd(eta) |
|---|---|---|---|---|---|
| 7 | blocks of 10, rho = 0.5 | 400 | (1, 1, 1, -1) at x1-4, x11-14, x21-24 | +0.50 (3 conflicts) | 3.46 |
| 8 | blocks of 10, rho = 0.7 | 400 | (1, 1, 1, -1.5) at x1-4, x11-14, x21-24 | +0.60 (3 conflicts) | 3.07 |
| 9 | blocks of 10, rho = 0.8 | 400 | (2, 2, -2) at x1-3, x11-13, x21-23 | +1.20 (3 conflicts) | 4.10 |
| 10 | AR(1), rho = 0.8 | 400 | (1, 1, 1, -1.5) at x1-4, x101-104, x201-204 | +0.45 (3 conflicts) | 3.41 |
| 11 | blocks of 20, rho = 0.9 (control) | 100 | (1, 1) at x1-2, x21-22, x41-42, x61-62, x81-82 | none | 4.36 |

In every block design, the noise variables inside an active block have nonzero univariate
slopes. These "proxies" number 18 for model 7, 18 for model 8, 21 for model 9 and 90 for
model 11; the AR(1) design of model 10 has 312. Model 11 has no sign conflicts. It checks
that suni2D loses nothing when the univariate signs are right.

## Results (200 replicates)

The full tables (all 7 criteria) are in `binomial_sign_conflict_tables.pdf`. All 200 x 5
replicates finished without errors or convergence messages. Seeds 1-30 are identical to the
earlier 30-replicate runs of the same designs.

rpe_kl, mean (se):

| method | model 7 | model 8 | model 9 | model 10 | model 11 |
|---|---|---|---|---|---|
| suni2D_zeta_tau | **0.294** (0.006) | **0.252** (0.005) | **0.208** (0.005) | **0.342** (0.006) | **0.558** (0.016) |
| suni2D_lambda_alpha | 0.294 (0.006) | 0.252 (0.005) | 0.208 (0.005) | 0.342 (0.006) | 0.558 (0.016) |
| suni2D_sum_rate | 0.294 (0.006) | 0.252 (0.005) | 0.208 (0.005) | 0.342 (0.006) | 0.559 (0.016) |
| uniLasso | 0.442 (0.009) | 0.488 (0.008) | 0.575 (0.010) | 0.489 (0.009) | 0.597 (0.019) |
| lasso | 0.386 (0.006) | 0.425 (0.006) | 0.541 (0.008) | 0.410 (0.006) | 0.631 (0.016) |
| adalasso_uni | 0.333 (0.007) | 0.275 (0.006) | 0.228 (0.006) | 0.371 (0.008) | 0.598 (0.019) |
| adalasso_lasso | 0.429 (0.010) | 0.472 (0.011) | 0.467 (0.014) | 0.432 (0.009) | 1.075 (0.034) |

The next table gives the paired difference in rpe_kl, suni2D_zeta_tau minus the competitor:
the mean (se) over the 200 replicates, and in brackets the share of replicates where suni2D is
lower.

| competitor | model 7 | model 8 | model 9 | model 10 | model 11 |
|---|---|---|---|---|---|
| uniLasso | -0.148 (0.007) [0.97] | -0.236 (0.007) [0.99] | -0.367 (0.008) [1.00] | -0.147 (0.007) [0.97] | -0.039 (0.010) [0.59] |
| lasso | -0.091 (0.005) [0.92] | -0.172 (0.005) [0.98] | -0.333 (0.006) [1.00] | -0.068 (0.006) [0.78] | -0.074 (0.012) [0.76] |
| adalasso_uni | -0.039 (0.005) [0.70] | -0.022 (0.004) [0.64] | -0.021 (0.003) [0.61] | -0.029 (0.005) [0.64] | -0.040 (0.010) [0.60] |
| adalasso_lasso | -0.134 (0.009) [0.81] | -0.220 (0.009) [0.96] | -0.260 (0.013) [0.89] | -0.090 (0.009) [0.76] | -0.517 (0.030) [0.92] |

Other criteria, mean over the 200 replicates:

| | model 7 | model 8 | model 9 | model 10 | model 11 |
|---|---|---|---|---|---|
| support: suni2D / uniLasso / lasso / adalasso_uni | 32 / 42 / 63 / 61 | 40 / 35 / 64 / 59 | 32 / 28 / 53 / 45 | 30 / 37 / 58 / 58 | 18 / 21 / 28 / 21 |
| tpr: suni2D / uniLasso / lasso / adalasso_uni | 0.91 / 0.75 / 0.89 / 0.97 | 0.98 / 0.71 / 0.86 / 0.99 | 1.00 / 0.67 / 0.80 / 1.00 | 0.74 / 0.54 / 0.68 / 0.88 | 0.30 / 0.32 / 0.32 / 0.32 |
| fpr: suni2D / uniLasso / lasso / adalasso_uni | 0.021 / 0.033 / 0.053 / 0.050 | 0.029 / 0.026 / 0.054 / 0.048 | 0.023 / 0.022 / 0.047 / 0.036 | 0.021 / 0.031 / 0.050 / 0.048 | 0.015 / 0.018 / 0.025 / 0.018 |

In short:

- suni2D has the lowest rpe_kl and rpe_class (except a tie in model 11) in every model.
- uniLasso essentially never selects the conflict variables. In model 9 its tpr is 0.665,
  about 6/9, the share of signals whose univariate sign agrees.
- adalasso_uni finds the conflict variables but selects about twice as many variables.
- suni2D's slope `sse` is larger than adalasso_uni's (more shrinkage), but its prediction is
  better.
- The three suni2D parameterizations agree to three decimals.

Median time per replicate on one Argon core: 5-8 minutes for models 7-10, 1.4 minutes for
model 11. Each model ran as one job (40 cores for models 7-10, 16 cores for model 11) of
18-43 minutes.

## Output format

`res_model<k>.rds` is a list with these elements:

- `table`, `mean`, `se`: methods x criteria, in the format of `summarize_sim()`, as in
  `inst/simulation/binomial/res_model1-6.rds`.
- `raw`: replicates x methods x criteria array of the per-replicate values, useful for paired
  comparisons.
- `n` and `seeds`: the number of replicates and their seeds.
- `settings`: the model definition.
- `messages`: distinct warnings, for example solver convergence notes.
- `time`: seconds per replicate.

`table_simulation.R` reads these files directly; `make_tables.R` shows how.

## How to run

From the package root (the folder that contains `inst/`):

```sh
# everything, 1 core (Windows always uses 1 core): about 5-10 min per n = 400 replicate
Rscript inst/simulation/binomial_sign_conflict/sign_conflict_simulation.R
# a subset: models 7 and 9, seeds 1-20
Rscript inst/simulation/binomial_sign_conflict/sign_conflict_simulation.R "c(7,9)" 1:20
# arguments: [models] [reps] [ncores] [sim_dir] [out_dir]
#   defaults: 7:11  1:200  1  inst/simulation  inst/simulation/binomial_sign_conflict/rerun
Rscript inst/simulation/binomial_sign_conflict/make_tables.R   # pdf of res_model7-11.rds
```

Local runs write to `binomial_sign_conflict/rerun/` by default, so the shipped results are
never overwritten. Each replicate is saved as `<out_dir>/reps/m<k>_b<bbb>.rds`, written under a
temporary name and renamed when complete. A rerun with the same arguments skips every saved
replicate, then rebuilds `<out_dir>/res_model<k>.rds` from all the saved replicates. An
interrupted run, for example by an Argon maintenance, therefore continues where it stopped.

**Argon.** The 200-replicate run copied `sign_conflict_simulation.R`, `compare_binom.R` and
`binomial_simulation_setup.R` into `~/sulnet/binomial_simulation/scripts`. One job per model
was submitted from `~/sulnet/binomial_simulation`:

```sh
qsub -N binom_sc_m7 runjob_sign_conflict.txt 7      # same for 8, 9, 10
qsub -N binom_sc_m11 -pe smp 16 runjob_sign_conflict.txt 11
```

Results are in `~/sulnet/binomial_simulation/results/sign_conflict/`: the per-replicate files
in `reps/` and `res_model7-11.rds`. Resubmitting the same job resumes it. A local run
reproduced the Argon values exactly (model 11, seed 1).

## Why these settings and models (grid-tuning study)

The study is in `grid_tuning/`. Each script states its usage and the Argon run it made. All
use 30 replicates (seeds 1-30), the same data and folds as the paper pipeline, and cv.min.
Each task saves the CV curves and test criteria of every grid cell, so `analyze.R` can
evaluate any sub-grid of the second coordinate (tau, alpha or rate) afterwards.

| script | what it ran | main finding |
|---|---|---|
| `round1.R` | models 1-6; default grids with `standardize` TRUE and FALSE; dense tau / alpha / rate grids with `standardize = FALSE` | `standardize = TRUE` makes suni2D about equal to the lasso (models 2-6). Standardizing F_j = beta_j x_j removes abs(beta_j), so tau = 0 becomes the plain lasso. With `FALSE`, tau = 1e3 zmax reproduces uniLasso (loo = FALSE) and tau = 0 reproduces adalasso_uni along the whole path. **tau has no effect in models 1-6**, because no true signal has a sign conflict. alpha > 0.5 or rate > 1 only lowers the effective zeta (model 1: lambda_alpha 2.82 vs zeta_tau 2.63 rpe_kl). Deviance is at least as good as Brier; 1-SE is worse except in models 3 and 5. |
| `round2.R zeta` | models 1-6, zeta grid variants (lambda.factor 0.001 / 0.03, nzeta 200) | The default zeta grid is fine. lambda.factor 0.001 hurts model 1 (+0.18); nothing changes in models 2-6. |
| `round2.R dense` | sign-conflict models 7-9 of `sim_lib.R` (ids 8 and 9 = models 8 and 9 here) | The best achievable rpe_kl falls steadily as tau goes to 0, and CV picks tau = 0 in 30/30 and 29/30 replicates. The default tau grid (lowest value 0.01 zmax, no 0) loses about 0.04 rpe_kl; any grid that contains 0 gives the same result. |
| `round3.R` | 12 designs R1-R12 (ids 11-22): sign conflicts of several strengths and numbers, AR(1), suppressor pairs, global equicorrelation, conflict-free proxies; dense tau grids with `standardize` FALSE and TRUE | **No intermediate tau beat the better endpoint in any design**, even with many correlated noise proxies and small n. `standardize = TRUE` is about equal to the lasso everywhere. Ids 15, 16 and 19 became models 7, 10 and 11. |
| `final_30reps.R` | the paper pipeline with the new settings, seeds 1-30: models 1-6, models 7-9 of `sim_lib.R`, and 7 of the round-3 designs | Models 1-6 are in `results_30reps/`. |
| `analyze.R` | post-hoc sub-grid evaluation of round 1-3 output | |

Round 3, rpe_kl at cv.min, 30 replicates (suni2D = zeta_tau, `standardize = FALSE`, full tau grid):

| id | design | n | suni2D | adalasso_uni | lasso | uniLasso |
|---|---|---|---|---|---|---|
| 11 (R1) | model 9 at n = 200 | 200 | **0.578** | 0.608 | 0.736 | 0.773 |
| 12 (R2) | model 9 with 6 conflict blocks | 400 | 0.615 | **0.592** | 1.108 | 1.064 |
| 13 (R3) | 3 conflict blocks + 3 blocks (1, 1, 1) | 400 | 0.540 | **0.534** | 0.834 | 0.728 |
| 14 (R4) | blocks of 20, rho = 0.9, (1.5, 1.5, -1.5) x 3 | 200 | **0.265** | 0.299 | 0.327 | 0.298 |
| 15 (R5) | = model 7 | 400 | **0.281** | 0.329 | 0.363 | 0.423 |
| 16 (R6) | = model 10 | 400 | **0.325** | 0.347 | 0.404 | 0.464 |
| 17 (R7) | suppressor pairs (2, -1), rho = 0.9, x 5 | 200 | **0.344** | 0.516 | 0.393 | 0.520 |
| 18 (R8) | 3 suppressor pairs + 6 positives with noise partners | 200 | **0.529** | 0.692 | 0.558 | 0.692 |
| 19 (R9) | = model 11 | 100 | **0.545** | 0.577 | 0.661 | 0.576 |
| 20 (R10) | R9 blocks + 2 conflict blocks | 150 | **0.474** | 0.502 | 0.536 | 0.507 |
| 21 (R11) | global equicorrelation 0.3, 8 positives + 2 conflicts | 200 | 0.621 | 0.604 | **0.593** | 0.602 |
| 22 (R12) | 6 weak conflicts (true -1, univariate slope +1.4) | 400 | **0.397** | 0.402 | 0.573 | 0.474 |

Models 7-11 were chosen as the clearest cases: suni2D beats every competitor, including
adalasso_uni, and the design shows the mechanism (models 7-10) or rules out a loss (model 11).

**What drives the gains.** In the sign-conflict designs, CV picks tau = 0. That is the
unconstrained, univariately weighted fit, which uniLasso cannot reach. The remaining edge over
adalasso_uni comes from the CV. cv.suni2D refits the univariate stage inside every fold, while
uniLasso (`loo = FALSE`) and adalasso_uni reuse the full-data univariate fits and select too
many variables. The tau grid works as a data-driven choice between the sign-constrained and the
unconstrained use of the univariate information. Part of the margin over uniLasso exists
because it runs without its leave-one-out step here.

The raw round outputs (per-cell CV and test criteria) and their analyses (`analysis/*.csv`) are
on Argon in `~/sulnet/binomial_simulation/results/round1, round2a (zeta), round2b (dense),
round3, final, final_conflict`. They are not copied here because of their size.

## Notes

- **Argon job logs.** These show `Lmod ... cmake/3.25.0_gcc-9.5.0 ... cannot be loaded` when
  cmake is loaded before `stack/2022.2`. It is harmless (cmake is not needed to run R). The job
  files here load stack first.
- **Convergence messages.** A solver convergence message
  (`Convergence not reached within 100000 passes ...`) means some grid cells of that fit were
  not computed. Those cells are left out of the CV, as in `cv.suni2D`. Such messages are kept
  in `res_model<k>.rds$messages`.
