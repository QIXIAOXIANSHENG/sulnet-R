# Binomial models 1-6 with the tuned settings (200 replicates)

This folder reruns models 1-6 of `inst/simulation/binomial simulation.R`. The data and
seeds are the same as there; the method settings are the ones chosen in the grid-tuning study
in `inst/simulation/binomial_sign_conflict/` (see its README, "Why these settings and
models"). The original results in `inst/simulation/binomial/` are untouched.

The runs were done on the Argon cluster (University of Iowa) on 2026-11-11, with sulnet built
from GitHub `QIXIAOXIANSHENG/sulnet-R` at commit `b11feb6` (sulnet 0.4.0), R 4.2.2,
glmnet 5.2 and uniLasso 2.11.

## Contents

| file | what it is |
|---|---|
| `binomial_tuned_simulation.R` | Runs models 1-6 (simulation, fitting, summaries). Resumes interrupted runs. |
| `make_tables.R` | Writes `binomial_tuned_tables.pdf` from `res_model1-6.rds`. |
| `res_model1.rds` ... `res_model6.rds` | 200-replicate results (seeds 1-200), described below. |
| `binomial_tuned_tables.pdf` | The summary tables of models 1-6. |
| `argon/runjob_tuned.txt` | The SGE job used on Argon, with the submission commands. |
| `diagnostics/model1_cv_diagnostic.R` | Model 1 diagnostic: why suni2D trails uniLasso and adalasso_uni (section below). |
| `diagnostics/model1_cv_analysis.R` | Prints the numbers of that section from `diagnostics/model1_cv_diagnostic.rds`. |
| `diagnostics/model1_cv_diagnostic.rds` | Per-seed output of the diagnostic, seeds 1-200. |

## What changed from `binomial simulation.R`

Everything else (data, seeds, folds, the 7 methods, the 7 criteria) is the same.

1. **suni2D uses `standardize = FALSE`** (all three parameterizations). With the default
   `standardize = TRUE`, the univariate fits F_j = beta_j x_j + c are rescaled to unit
   variance, which removes abs(beta_j) and turns tau = 0 into the plain lasso. In the old
   tables all three suni2D variants therefore matched the lasso in models 2-6.
2. **zeta_tau tau grid** = zmax x {1e3, 10^0, 10^-0.25, ..., 10^-3, 0} (15 values), where
   zmax is the largest zeta of the default zeta grid. tau = 1e3 zmax reproduces uniLasso
   and tau = 0 the univariately weighted adaptive lasso. The default zeta grid is kept, and
   lambda_alpha and sum_rate keep their default grids (which contain alpha = 0.5 and
   rate = 1, the tau = 0 end).
3. **uniLasso uses `loo = FALSE`**: `cv.uniLasso(..., loo = FALSE)`. The adalasso_uni weights
   already came from `uniFit(..., loo = FALSE)`.
4. **True intercept.** Every model is scored against its true intercept:
   `gamma0 = dat$intercept` is passed to `run_simulation_binom()` for models 2-6, and 0 for
   model 1. Only model 5 has a nonzero intercept (-2.19154, which gives prevalence 0.2); the
   others have prevalence 0.5 and intercept 0. The test probabilities used by `rpe_kl` and
   `rpe_class` are plogis(gamma0 + x'gamma).

CV loss is deviance (`pred.loss = "loss"`, the default) and every method is evaluated at the
CV minimum (`s = "min"`), as before.

## Settings

**Pipeline.** For replicate b = 1, ..., 200 the data are drawn from seed b exactly as in
`binomial simulation.R`. `run_simulation_binom()` from `compare_binom.R` is then called with
`B = 1` and `s = "min"`; it draws one assignment of 10 class-stratified CV folds right after
the data, shared by all methods.

**Models.** p Gaussian covariates with unit variances; y ~ Bernoulli(plogis(gamma0 + x'gamma)).

| model | design | n | p | n test | nonzero coefficients | intercept | sd(x'gamma) |
|---|---|---|---|---|---|---|---|
| 1 | uniLasso case: AR(1) rho = 0.8; in training, x1-x20 shifted by +0.5 when a fair coin `ybase` is 1 | 100 | 300 | 30 | x1-x30 = 1 | 0 | 15.2 (+10 shift) |
| 2 | independent | 200 | 1000 | 100 | x1-x10 = (1, -1) x 5 | 0 | 3.16 |
| 3 | AR(1) rho = 0.5 | 200 | 1000 | 100 | x1-x10 = (1, -1) x 5 | 0 | 1.94 |
| 4 | blocks of 20, rho = 0.5 | 200 | 1000 | 100 | x1-x10 = (1, -1) x 5 | 0 | 2.24 |
| 5 | imbalanced: AR(1) rho = 0.5, prevalence 0.2 | 200 | 1000 | 100 | x1-x10 = (1, -1) x 5 | -2.19 | 1.94 |
| 6 | Breheny-Huang (2011), independent | 400 | 1000 | 100 | x1-x4 = (1, 1, -1, -1) | 0 | 2.00 |

Models 2-5 use `simulation_setting()` and model 6 `bh2011()` from
`binomial_simulation_setup.R`. Model 1's test covariates have no mean shift.

**Methods** (the rows of the tables):

| method | fit |
|---|---|
| `suni2D_zeta_tau` | `cv.suni2D(method = "zeta_tau", family = "binomial", standardize = FALSE)`; default zeta grid (100 values); tau grid above |
| `suni2D_lambda_alpha` | same with `method = "lambda_alpha"`, default grid (100 lambda x 20 alpha in [0.5, 0.99]) |
| `suni2D_sum_rate` | same with `method = "sum_rate"`, default grid (100 sum x 20 rate in [1, 100]) |
| `uniLasso` | `cv.uniLasso(family = "binomial", loo = FALSE)` |
| `lasso` | `cv.glmnet(family = "binomial")` |
| `adalasso_uni` | adaptive lasso, weights 1/abs(univariate slopes) from `uniFit(loo = FALSE)`, `standardize = FALSE` |
| `adalasso_lasso` | adaptive lasso, weights 1/abs(lasso coefficients) |

**Criteria** (`eval_criteria_binom()` in `compare_binom.R`): `sse` (squared error of the
slopes), `rpe_kl` (Kullback-Leibler risk of the test probabilities relative to the entropy;
main criterion), `rpe_class` (excess test misclassification relative to the Bayes error),
`support`, `sign` (share of coordinates with the true sign, 0 counting as a sign), `tpr`, `fpr`.

## Results (200 replicates)

The full tables (all 7 criteria) are in `binomial_tuned_tables.pdf`.

- All 1200 replicates finished, and no criterion is NA.
- 5 replicates had a solver convergence message: model 3 seeds 48 and 131, model 4 seed 164, and
  model 5 seeds 24 and 93. In each, some grid cells of one suni2D fit were not computed. Those
  cells are left out of the CV, as `cv.suni2D` does.
- Seeds 1-30 are identical to the earlier 30-replicate run with the same settings
  (`inst/simulation/binomial_sign_conflict/grid_tuning/results_30reps`).

rpe_kl, mean (se); bold = lowest (within 0.0005):

| method | model 1 | model 2 | model 3 | model 4 | model 5 | model 6 |
|---|---|---|---|---|---|---|
| suni2D_zeta_tau | 2.312 (0.111) | 0.536 (0.011) | **0.474** (0.005) | 0.579 (0.007) | 0.497 (0.007) | **0.050** (0.001) |
| suni2D_lambda_alpha | 2.402 (0.127) | **0.535** (0.010) | **0.474** (0.005) | 0.579 (0.007) | 0.497 (0.007) | **0.050** (0.001) |
| suni2D_sum_rate | 2.402 (0.127) | **0.535** (0.010) | **0.474** (0.005) | 0.579 (0.007) | 0.497 (0.007) | **0.050** (0.001) |
| uniLasso | 2.247 (0.107) | 0.707 (0.017) | 0.825 (0.019) | 0.737 (0.014) | 0.978 (0.028) | 0.218 (0.009) |
| lasso | 2.791 (0.124) | 0.558 (0.008) | **0.473** (0.005) | **0.543** (0.007) | **0.496** (0.007) | 0.107 (0.002) |
| adalasso_uni | **2.243** (0.108) | 0.706 (0.017) | 0.824 (0.019) | 0.729 (0.014) | 0.983 (0.028) | 0.217 (0.009) |
| adalasso_lasso | 3.323 (0.146) | 0.932 (0.026) | 0.811 (0.030) | 0.959 (0.029) | 0.909 (0.043) | 0.158 (0.007) |

The next table gives the paired difference in rpe_kl, suni2D_zeta_tau minus the competitor:
the mean (se) over the 200 replicates, and in brackets the share of replicates where suni2D is
lower.

| competitor | model 1 | model 2 | model 3 | model 4 | model 5 | model 6 |
|---|---|---|---|---|---|---|
| uniLasso | +0.065 (0.022) [0.38] | -0.171 (0.015) [0.80] | -0.351 (0.018) [0.99] | -0.157 (0.013) [0.82] | -0.480 (0.026) [0.96] | -0.168 (0.009) [0.96] |
| lasso | -0.479 (0.070) [0.74] | -0.022 (0.005) [0.65] | +0.000 (0.002) [0.46] | +0.036 (0.004) [0.24] | +0.002 (0.003) [0.47] | -0.057 (0.001) [1.00] |
| adalasso_uni | +0.068 (0.021) [0.36] | -0.170 (0.015) [0.80] | -0.350 (0.018) [0.99] | -0.150 (0.012) [0.82] | -0.485 (0.026) [0.96] | -0.167 (0.009) [0.96] |
| adalasso_lasso | -1.012 (0.088) [0.81] | -0.396 (0.023) [0.90] | -0.338 (0.029) [0.82] | -0.380 (0.028) [0.85] | -0.412 (0.041) [0.83] | -0.108 (0.007) [0.91] |

Other criteria, mean over the 200 replicates (suni2D = suni2D_zeta_tau):

| | model 1 | model 2 | model 3 | model 4 | model 5 | model 6 |
|---|---|---|---|---|---|---|
| support: suni2D / uniLasso / lasso / adalasso_uni | 16 / 17 / 29 / 17 | 28 / 60 / 54 / 60 | 5 / 59 / 13 / 59 | 10 / 55 / 30 / 55 | 4 / 50 / 10 / 51 | 10 / 70 / 40 / 70 |
| tpr: suni2D / uniLasso / lasso / adalasso_uni | 0.47 / 0.48 / 0.49 / 0.48 | 0.87 / 0.95 / 0.95 / 0.95 | 0.13 / 0.42 / 0.19 / 0.42 | 0.23 / 0.67 / 0.47 / 0.66 | 0.08 / 0.31 / 0.12 / 0.32 | 1.00 / 1.00 / 1.00 / 1.00 |
| fpr: suni2D / uniLasso / lasso / adalasso_uni | 0.008 / 0.011 / 0.052 / 0.011 | 0.019 / 0.051 / 0.045 / 0.051 | 0.004 / 0.056 / 0.011 / 0.056 | 0.008 / 0.049 / 0.026 / 0.048 | 0.003 / 0.048 / 0.009 / 0.048 | 0.006 / 0.066 / 0.036 / 0.066 |
| rpe_class: suni2D / uniLasso / lasso / adalasso_uni | 2.68 / 2.52 / 3.28 / 2.52 | 0.73 / 0.75 / 0.72 / 0.75 | 1.04 / 0.99 / 1.04 / 0.99 | 1.20 / 0.94 / 1.05 / 0.94 | 0.37 / 0.63 / 0.36 / 0.63 | 0.05 / 0.21 / 0.10 / 0.21 |

In short:

- **Model 6.** suni2D halves the lasso's rpe_kl (0.050 vs 0.107). It is lower in every
  replicate and selects about 10 variables, against 40 for the lasso and 70 for uniLasso.
- **Model 2.** suni2D beats the lasso by 0.022 (se 0.005).
- **Models 3 and 5.** suni2D ties with the lasso. It selects the sparsest models (5 and 4
  variables).
- **Model 4 (blocks).** This is the only clear loss to the lasso: +0.036 (se 0.004).
- **Model 1.** suni2D_zeta_tau is about 0.07 (se 0.02) worse than uniLasso and adalasso_uni,
  and 0.48 better than the lasso. The reason is a difference in how the CV is done, not in
  the fits; see "Model 1: why suni2D trails uniLasso and adalasso_uni" below. lambda_alpha
  and sum_rate are worse here (2.40), because alpha > 0.5 or rate > 1 only lowers the
  effective zeta.
- **uniLasso and adalasso_uni without LOO** select 50-70 variables in models 2-6.
- **No sign conflicts in models 1-6.** As in the grid-tuning study, no true signal has a
  univariate sign conflict, so CV has no reason to pick an intermediate tau. The three suni2D
  parameterizations agree within 0.001 rpe_kl in models 2-6.

Each 40-core job took 12-22 minutes. The median minutes per replicate were 0.6, 3.0, 3.5, 2.2,
2.2 and 4.3 for models 1-6.

## Comparison with the original results (`inst/simulation/binomial/`)

The original files used the default `standardize = TRUE` for suni2D and `loo = TRUE` for
uniLasso. rpe_kl means, original -> tuned:

| model | suni2D_zeta_tau | uniLasso | lasso | adalasso_uni |
|---|---|---|---|---|
| 1 | 2.941 -> 2.312 | 2.389 -> 2.247 | 2.953 -> 2.791 | 2.335 -> 2.243 |
| 2 | 0.558 -> 0.536 | 0.523 -> 0.707 | 0.558 (same) | 0.706 (same) |
| 3 | 0.474 -> 0.474 | 0.487 -> 0.825 | 0.473 (same) | 0.824 (same) |
| 4 | 0.543 -> 0.579 | 0.574 -> 0.737 | 0.543 (same) | 0.729 (same) |
| 5 | 0.496 -> 0.497 | 0.539 -> 0.978 | 0.496 (same) | 0.983 (same) |
| 6 | 0.107 -> 0.050 | 0.052 -> 0.218 | 0.107 (same) | 0.217 (same) |

- **Models 2-6 use the same data and folds as before.** lasso, adalasso_uni and adalasso_lasso
  do not depend on the changed settings, and they reproduce the original means exactly. That
  includes model 5, so the original model-5 file was already scored with the true intercept.
- **`standardize = FALSE`** helps suni2D in models 1, 2 and 6, changes nothing in models 3 and
  5, and costs 0.036 in model 4.
- **uniLasso's leave-one-out step matters a lot.** With `loo = TRUE` it was the best method in
  model 2 (0.523) and close to the tuned suni2D in model 6 (0.052). With `loo = FALSE` its
  rpe_kl rises by 0.14-0.44 in models 2-6. Much of suni2D's margin over uniLasso in these
  tables therefore comes from running uniLasso without LOO.
- **Model 1 uses different random draws here.** Its covariates come from `MASS::mvrnorm()`,
  which builds the draws from `eigen(Sigma)`. The same seed gives different data on Argon
  (R 4.2.2) and on a Windows PC (R 4.5): seed 1 differs, while model 5 seed 1 agrees to 1e-15.
  The cause is most likely a difference in the eigenvectors (their signs) between the LAPACK
  builds. The original model-1 file was made on the PC, so the model-1 rows of the two tables
  come from different, equally distributed samples. For the same reason, a local rerun of
  model 1 will not reproduce these values exactly. Models 2-6 use `rnorm()` and reproduce on
  any machine.

## Model 1: why suni2D trails uniLasso and adalasso_uni

Along the path, tau = 1e3 zmax gives the same fits as uniLasso and tau = 0 the same as
adalasso_uni, and the tau grid contains both. Yet in model 1 suni2D_zeta_tau is about 0.07
worse than both (paired table above). `diagnostics/model1_cv_diagnostic.R` takes the gap
apart seed by seed, with the model 1 data and folds of the pipeline. Its own CV reproduces the
choices of `cv.suni2D` exactly. `diagnostics/model1_cv_analysis.R` prints the numbers below
from `diagnostics/model1_cv_diagnostic.rds` (seeds 1-200).

The diagnostic was run on a Windows PC, so its model 1 draws differ from the Argon run (see
above). The gap is the same size: +0.059 vs uniLasso here, +0.065 on Argon.

**Ruled out.**

- **The paths are the same.** The best test rpe_kl anywhere along the path agrees within
  0.004 (se 0.004) between the tau = 1e3 zmax fit and uniLasso with converged univariate fits
  (step 1 below), and between the tau = 0 fit and adalasso_uni.
- **The choice of tau costs nothing.** Fits with only tau = 1e3 zmax, or only tau = 0, give
  the same result as the full tau grid in 99% of seeds. All the relevant univariate slopes in
  model 1 are positive, so the sign penalty never binds.
- **No path stopped early.** In no seed did a fold path saturate. So the `cv.suni2D` rule that
  leaves out a cell any fold did not compute (`cv.glmnet` instead reuses that fold's last fit)
  never applied, and filling such cells as `cv.glmnet` does changes nothing.

**Where the gap comes from.** rpe_kl, mean over the 200 seeds:

| step | comparison | difference (se) |
|---|---|---|
| total | suni2D_zeta_tau (full tau grid) minus uniLasso | +0.059 (0.023) |
| 1 | uniLasso with converged univariate fits, minus uniLasso | +0.020 (0.008) |
| 2 | the same with honest CV, minus with reused univariate fits | +0.035 (0.021) |
| rest | suni2D (tau = 1e3 zmax only) minus uniLasso with both changes | +0.004 (0.008) |

Against adalasso_uni the total is +0.039 (0.024), which is essentially step 2 alone: with
converged univariate fits, uniLasso is identical to adalasso_uni in 95% of seeds. Likewise,
the tau = 0 fit ties adalasso_uni with honest CV (+0.005, se 0.020).

1. **uniLasso's univariate fits are not fully converged.** `cv.uniLasso` gets them from
   `uniLasso::uniInfo()` with `nit = loob.nit = 2` Newton steps. Their slopes are about 5%
   larger than the converged ones (`sulnet::uniFit()`) that suni2D and adalasso_uni use, and
   in 98% of seeds some slope differs by a factor outside [0.8, 1.25]. In model 1 the
   unconverged version happens to predict slightly better: the best test rpe_kl along its
   path is 2.044, against 2.065 with converged fits.
2. **Honest vs reused univariate fits in the CV.** `cv.suni2D` redoes the univariate fits
   inside every training fold. uniLasso (`loo = FALSE`) and adalasso_uni compute them once on
   all the data and reuse them in every fold, so their CV has already seen the held-out
   observations. Their CV deviance at its minimum is optimistic (0.374, against 0.466 for the
   honest CV), and they choose smaller penalties: a median of 0.0066 against 0.0096 for the
   honest CV, which picks the larger penalty in 69% of seeds. Model 1's signal is very strong
   (sd of x'gamma 15, with 30 correlated true coefficients), so the test-optimal penalty is
   small (median 0.0085), and the optimism happens to help.

**Conclusion.** suni2D is not doing anything wrong in model 1: it ties uniLasso when uniLasso
gets the same honest CV. Refitting the univariate stage inside the folds is the correct CV for
a two-stage method. Usually the reused fits hurt: they are why uniLasso and adalasso_uni
select 50-70 variables in models 2-6 and lose clearly there. uniLasso's leave-one-out step
(`loo = TRUE`) exists to remove this leak, which is also why it did so much better in the
original tables. A like-for-like comparison would give the baselines honest CV (or
`loo = TRUE` for uniLasso); that would also change their results in models 2-6.

To rerun the diagnostic (from the package root; about 3 minutes on 20 cores, and it also runs
on Windows):

```sh
Rscript inst/simulation/binomial_tuned/diagnostics/model1_cv_diagnostic.R 1:200 20
# arguments: [seeds] [ncores] [out_dir]; default out_dir diagnostics/rerun
Rscript inst/simulation/binomial_tuned/diagnostics/model1_cv_analysis.R \
  inst/simulation/binomial_tuned/diagnostics/rerun/model1_cv_diagnostic.rds
```

## Output format

`res_model<k>.rds` is a list with these elements:

- `table`, `mean`, `se`: methods x criteria, in the format of `summarize_sim()`, as in
  `inst/simulation/binomial/res_model1-6.rds`.
- `raw`: replicates x methods x criteria array of the per-replicate values, useful for paired
  comparisons.
- `n` and `seeds`: the number of replicates and their seeds.
- `settings`: the model and method settings.
- `messages`: distinct warnings, for example solver convergence notes.
- `time`: seconds per replicate.

`table_simulation.R` reads these files directly; `make_tables.R` shows how.

## How to run

From the package root (the folder that contains `inst/`):

```sh
# everything, 1 core (Windows always uses 1 core)
Rscript inst/simulation/binomial_tuned/binomial_tuned_simulation.R
# a subset: models 2 and 5, seeds 1-20
Rscript inst/simulation/binomial_tuned/binomial_tuned_simulation.R "c(2,5)" 1:20
# arguments: [models] [reps] [ncores] [sim_dir] [out_dir]
#   defaults: 1:6  1:200  1  inst/simulation  inst/simulation/binomial_tuned/rerun
Rscript inst/simulation/binomial_tuned/make_tables.R   # pdf of res_model1-6.rds
```

Local runs write to `binomial_tuned/rerun/` by default, so the shipped results are never
overwritten. Each replicate is saved as `<out_dir>/reps/m<k>_b<bbb>.rds`, written under a
temporary name and renamed when complete. A rerun with the same arguments skips every saved
replicate, then rebuilds `<out_dir>/res_model<k>.rds` from all the saved replicates, so an
interrupted run continues where it stopped. Minutes per replicate on one core: about 0.6 for
model 1, 2-4 for models 2-5 and 4-10 for model 6. A local run of models 2-6 reproduces the
Argon values (checked for model 5, seed 1); model 1 does not (see above).

**Argon.** `binomial_tuned_simulation.R`, `compare_binom.R` and `binomial_simulation_setup.R`
were copied into `~/sulnet/binomial_simulation/scripts`, and five 40-core jobs were submitted
from `~/sulnet/binomial_simulation` (see `argon/runjob_tuned.txt`):

```sh
qsub -N binom_tuned_m6  runjob_tuned.txt 6
qsub -N binom_tuned_m12 runjob_tuned.txt "c(1,2)"
qsub -N binom_tuned_m3  runjob_tuned.txt 3      # same for 4 and 5
```

Results are in `~/sulnet/binomial_simulation/results/tuned/`: the per-replicate files in
`reps/` and `res_model1-6.rds`. Resubmitting the same job resumes it. Jobs 6746068-6746072
ran on 2026-11-11 and finished without errors.
