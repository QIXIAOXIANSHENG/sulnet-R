# Linear models 1-5 with fixed glmnet paths and honest-CV adaptive lassos (200 replicates)

This folder reruns models 1-5 of `inst/simulation/linear simulation.R` with the corrected
`inst/simulation/compare_linear.R` (2026-11-12). The data, the replicate responses and the CV
folds are the same as in the original run: the suni2D rows reproduce the original per-replicate
results exactly (checked on models 1-3, whose original results keep the per-replicate values).
Only the glmnet-based methods (uniLasso, lasso, adalasso_uni, adalasso_lasso) changed. The
original results in `inst/simulation/linear/` are untouched.

Run locally on 2026-11-12 (Windows 11, 20 parallel workers) with R 4.5.0, sulnet 0.4.0 (repo
commit `f2c915b` plus the change to `compare_linear.R`), glmnet 4.1.10, uniLasso 2.11 and
MASS 7.3.65.

## Contents

| file | what it is |
|---|---|
| `linear_fixed_simulation.R` | Runs models 1-5 (data, fitting, summaries, tables). Resumes interrupted runs. |
| `res_model1.rds` ... `res_model5.rds` | 200-replicate results: `mean`, `se`, `table`, `n`, `raw` (200 x 7 x 6), `reps`, `seconds`. |
| `linear_fixed_tables.pdf` | The summary tables of models 1-5. |
| `reps/m<k>_b<bbb>.rds` | One file per model and replicate (criteria and run time). |

## What changed, and why

The old tables showed adalasso_lasso far behind the other methods (model 2: SSE 43.0, close
to ||gamma||^2 = 43.75), and lasso and uniLasso far behind suni2D in model 4. A diagnosis found
three causes, all in how glmnet was called rather than in the methods:

1. **glmnet stopped its lambda path early.** By default glmnet ends the path once the deviance
   ratio passes 0.999 or improves by less than 1e-5 per step (`glmnet.control(devmax, fdev)`).
   In model 4 (sigma = 0.1, R^2 about 0.9999) this happens before the lambda that CV wants,
   for lasso, uniLasso and both adaptive lassos. `fit_methods()` now sets `fdev = 0,
   devmax = 1` while it runs and restores the previous values on exit.
2. **The default lambda path is too short.** glmnet's default `lambda.min.ratio` is 0.01 when
   n < p. lasso and uniLasso now use `lambda.min.ratio = 1e-4`. The adaptive penalty factors
   1/|b_j| span up to five decades, so with the default path the variables with large weights
   never entered and CV picked the last lambda in every replicate of model 2. The adaptive lassos
   now use 200 lambdas down to 1e-6 lambda_max.
3. **CV of the adaptive lassos reused weights from the full data.** The first stage (the CV-tuned
   lasso, or the univariate slopes) was fitted on all observations, so every held-out fold had
   already influenced the weights and CV underestimated the error of small lambdas. The adaptive
   lassos are now tuned by honest 10-fold CV: in each fold the first stage is refit on the
   training folds (the lasso first stage with the remaining 9 folds as its own CV folds), and
   the final fit uses the first stage on all data. To keep one lambda sequence meaningful across
   folds, the adaptive lasso is fitted as a plain lasso on the rescaled columns x_j |b_j|
   (standardize = FALSE), which is the same problem as `penalty.factor = 1 / |b_j|`.

Everything else (the 7 methods, the 6 criteria, `s = "min"`, 10 folds shared by all methods)
is as before.

## Results compared with the original run

The full run took about 23 minutes (models 1-3 and 5: 11-20 s per replicate on one core,
model 4: about 77 s). All 1000 replicates finished without errors or warnings. The suni2D rows
are identical to the original tables in all five models.

Mean SSE / rpe of the methods that changed, old -> new:

| model | uniLasso | lasso | adalasso_uni | adalasso_lasso |
|---|---|---|---|---|
| 1 | 1.62 / 0.387 -> 1.53 / 0.426 | 1.29 / 0.509 -> 1.34 / 0.525 | 1.45 / 0.440 -> 1.38 / 0.435 | 2.07 / 0.562 -> 1.49 / 0.406 |
| 2 | 6.42 / 0.232 -> 5.89 / 0.188 | 5.98 / 0.238 -> 6.09 / 0.262 | 6.43 / 0.245 -> 5.88 / 0.194 | 43.0 / 1.197 -> 9.81 / 0.317 |
| 3 | 26.7 / 0.188 -> 26.1 / 0.208 | 24.2 / 0.362 -> 24.2 / 0.365 | 24.3 / 0.184 -> 23.9 / 0.195 | 94.4 / 0.485 -> 56.7 / 0.355 |
| 4 | 0.239 / 25.0 -> 0.014 / 1.22 | 0.130 / 25.0 -> 0.016 / 0.83 | 0.262 / 32.5 -> 0.016 / 0.85 | 0.696 / 45.3 -> 0.003 / 0.16 |
| 5 | unchanged | unchanged | 7.94 / 0.139 -> 8.33 / 0.141 | 12.5 / 0.185 -> 10.8 / 0.185 |

* **Model 4** changed the most: every glmnet-based method was stopped early before. Now
  adalasso_lasso has the smallest SSE and rpe of all methods (suni2D lambda_alpha: 0.012 / 0.64).
* **adalasso_lasso** improves in every model. In models 2 and 3 it is still clearly the worst
  method in SSE and TPR. This part is intrinsic: with strongly correlated signals (rho = 0.95, or
  near-collinear groups) the lasso keeps one or two variables of a group and zeroes or nearly
  zeroes the others, and those weights shut out true signals. A diagnostic found the same gap at
  the oracle lambda. This is why Huang, Ma & Zhang (2008) use marginal (univariate) weights.
* **uniLasso in model 1**: rpe rose from 0.387 to 0.426 while SSE and TPR improved. With the old
  settings its CV minimum was at the end of the truncated path in 95% of replicates (checked on
  40 replicates), so the truncation acted as extra shrinkage. CV now chooses inside the path.
* **Model 5** (n > p, SNR = 1): lasso and uniLasso are unchanged, since there the default path
  was already long enough.
