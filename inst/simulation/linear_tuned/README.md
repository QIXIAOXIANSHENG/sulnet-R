# Linear models 1-10 with the tuned suni2D and uniLasso settings (200 replicates)

This folder holds a rerun of the linear simulations after the 2026-11-12 fix of sulnet's Gaussian
solvers. It covers models 1-5 of `inst/simulation/linear simulation.R`. It adds five
sign-conflict models (6-10), the linear counterparts of the binomial models in
`inst/simulation/binomial_sign_conflict`.

The suni2D methods are fitted with both `standardize = FALSE` and `standardize = TRUE`, and
uniLasso with and without LOO. Earlier results are untouched:

- `inst/simulation/linear/`, the original run;
- `inst/simulation/linear_fixed/`, the run with fixed glmnet paths and honest-CV adaptive lassos.

Run locally on 2026-11-12/13 (Windows 11, 20 parallel workers) with R 4.5.0, glmnet 4.1.10,
uniLasso 2.11 and MASS 7.3.65. The sulnet build was sulnet 0.4.0: repo commit `f2c915b` plus the
uncommitted solver fix in `src/suniZTNET.f90`, `src/suniLANET.f90` and `src/suniSRNET.f90`.

## Contents

| file | what it is |
|---|---|
| `linear_tuned_simulation.R` | Runs models 1-10 (data, fitting, summaries, tables). Resumes interrupted runs. |
| `fit_methods_tuned.R` | The 11 methods (`fit_methods_tuned()`), built on the helpers of `compare_linear.R`. |
| `res_model1.rds` ... `res_model10.rds` | 200-replicate results: `mean`, `se`, `table`, `n`, `raw` (200 x 11 x 6), `reps`, `seconds`. |
| `linear_tuned_tables.pdf` | The summary tables of models 1-10. |
| `reps/m<k>_b<bbb>.rds` | One file per model and replicate (criteria and run time). |

To rerun or resume, from the package root:
`Rscript inst/simulation/linear_tuned/linear_tuned_simulation.R 1:10 1:200 <ncores>`.

## Why this rerun

In the earlier tables suni2D_zeta_tau did badly in model 4 (rpe 5.1). There were three
reasons, and all three are dealt with here:

1. **The default zeta grid was too short.** It ends at 0.01 zmax (`lambda.factor = 0.01` when
   n < p). With sigma = 0.1 the penalty that CV wants is far below that. CV picked the last
   zeta in every replicate, and the best cell of the grid was no better. lambda_alpha
   avoided it by choosing alpha near 1, which makes its positive-side threshold
   lambda (1 - alpha) tiny. All suni2D fits here use `lambda.factor = 1e-4`. lasso and
   uniLasso already use `lambda.min.ratio = 1e-4`.
2. **With `standardize = TRUE`, suni2D is the lasso.** Standardizing F_j = beta_j x_j + c
   removes |beta_j|, so the univariate weights are lost; only the sign penalty remains.
   With `standardize = FALSE`, tau = 0 gives the adaptive lasso with univariate weights.
   Both versions are reported here.
3. **With `standardize = FALSE`, the Gaussian solvers stopped too early.** They tested
   convergence on the raw squared change of a coefficient. glmnet, and sulnet's own binomial
   solvers, weight that change by the column's mean square. The columns of F have variance
   cor_j^2 var(y), so with `standardize = FALSE` the default `eps = 1e-8` was 20-60 times
   looser. The solvers also took half coordinate steps (`maj = 2 * var`). The fix takes exact
   coordinate steps and weights the stopping test by the column mean square, which is the mean
   squared change of the fitted values.
   - On a full grid it reaches the same optimum as the old solver at `eps = 1e-15`.
   - In model 4 with `standardize = FALSE`, the median prediction gap at the default `eps` fell
     from 5.4 sigma^2 to 0.01 sigma^2.
   - It is within 0.01 sigma^2 of glmnet on the tau = 0 lasso.

## Methods (rows of the tables)

All methods share the 10 CV folds of a replicate and are evaluated at the CV minimum.

| row | fit |
|---|---|
| `suni2D_zeta_tau`, `suni2D_lambda_alpha`, `suni2D_sum_rate` | `cv.suni2D(method = ..., standardize = FALSE, lambda.factor = 1e-4)`, default grids otherwise (100 x 20), default `eps` |
| `suni2D_zeta_tau_std`, `suni2D_lambda_alpha_std`, `suni2D_sum_rate_std` | the same with `standardize = TRUE` |
| `uniLasso` | `cv.uniLasso(loo = FALSE, lambda.min.ratio = 1e-4)` |
| `uniLasso_loo` | `cv.uniLasso(loo = TRUE, lambda.min.ratio = 1e-4)`, the `uniLasso` row of `linear_fixed` |
| `lasso` | `cv.glmnet(lambda.min.ratio = 1e-4)` |
| `adalasso_uni` | adaptive lasso with weights 1/abs(univariate slopes), honest CV (`compare_linear.R`) |
| `adalasso_lasso` | adaptive lasso with weights 1/abs(CV lasso coefficients), honest CV (`compare_linear.R`) |

glmnet's early path stops are off (`glmnet.control(fdev = 0, devmax = 1)`) during the fits,
as in `compare_linear.R`. With `lambda.factor = 1e-4` the default tau grid runs from zmax down
to 1e-4 zmax, so it covers both a nearly strict sign constraint and nearly no sign penalty.
A grid with an exact 0 gave the same results in a 20-replicate check.

The criteria are those of `compare_linear.R`: `sse`, `rpe` (relative prediction error),
`support`, `sign`, `tpr` and `fpr`.

## Models

| model | design | n, p | coefficients | noise |
|---|---|---|---|---|
| 1 | Huang et al., AR(1) rho = 0.5 within the signals | 100, 200 | 5 x 2.5, 5 x 1.5, 5 x 0.5 | sigma = 1.5 |
| 2 | the same with rho = 0.95 | 100, 200 | as model 1 | sigma = 1.5 |
| 3 | Huang et al., 3 groups of 5 near-collinear signals | 100, 200 | 15 x 1.5 | sigma = 1.5 |
| 4 | uniLasso, equicorrelated rho = 0.5 | 300, 1000 | 30 x U(0, 1) | sigma = 0.1 |
| 5 | uniLasso home court, AR(1) rho = 0.8 | 100, 30 | 20% nonzero | SNR = 1 |
| 6 | sign conflict: 3 blocks of 4 (rho = 0.5), other covariates independent | 200, 500 | (1, 1, 1, -1) in each block | sigma = 2 (SNR 3) |
| 7 | sign conflict: blocks of 10, rho = 0.5 | 400, 1000 | (1, 1, 1, -1) in blocks 1-3 | SNR 3 |
| 8 | sign conflict: blocks of 10, rho = 0.7 | 400, 1000 | (1, 1, 1, -1.5) in blocks 1-3 | SNR 3 |
| 9 | sign conflict: blocks of 10, rho = 0.8 | 400, 1000 | (2, 2, -2) in blocks 1-3 | SNR 3 |
| 10 | sign conflict: AR(1) rho = 0.8 | 400, 1000 | (1, 1, 1, -1.5) at x1-x4, x101-x104, x201-x204 | SNR 3 |

- **Models 1-5:** the data, replicate responses and CV folds are those of the earlier runs.
  For replicate 1 of every model, lasso, uniLasso_loo and both adaptive lassos reproduce
  `linear_fixed` exactly.
- **Sign conflict:** in an equicorrelated block, cov(x_j, x'gamma) = gamma_j + rho * (sum of the
  other coefficients in the block). So the negative coefficient has a positive univariate
  slope (+0.5, +0.6, +1.2; +0.45 for the AR(1) model 10). uniLasso and suni2D with a large tau
  force it to be >= 0.
- **Models 7-10:** these are the binomial models 7-10 with a Gaussian response.
  `simulate_binomial(seed = b)` gives the same x, test x and coefficients as the binomial run.
  y = x'gamma + sigma eps is drawn afterwards, with sigma = sd(x'gamma) / sqrt(3).
- **Model 6:** the design of the 20-replicate check that preceded this run, now drawn with
  seed b.
- **Seeds:** models 5-10 draw replicate b from seed b, with the folds drawn right after the data.

## Results

**The run.**
- All 2000 replicates finished. No fit failed and no fit raised a warning. The only warnings
  in the logs say that glmnet and uniLasso were built under newer R versions. The logs are
  not kept.
- Wall time was 11.3 hours on 20 workers. Models 7-10 took 9.4 of those hours.
- Median time per replicate (all 11 methods), in seconds: 30, 36, 64, 181, 22, 126, 622, 860,
  903 and 944 for models 1-10.
- The machine stopped for about 37 minutes during model 1 (all 20 workers paused at once). So
  20 model-1 replicates show about 2280 s in `seconds`; their actual compute time is about
  30 s.

**rpe** (mean over 200 replicates; the best method of each model is in bold). Standard errors
are 0.002-0.025; see the PDF. Column key:
- `zt`, `la`, `sr`: suni2D zeta_tau, lambda_alpha and sum_rate with `standardize = FALSE`;
- `_std`: the same with `standardize = TRUE`;
- `uni`, `uni_loo`: uniLasso without and with LOO;
- `ada_uni`, `ada_lasso`: the two adaptive lassos.

| model | zt | la | sr | zt_std | la_std | sr_std | uni | uni_loo | lasso | ada_uni | ada_lasso |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.438 | 0.438 | 0.439 | 0.531 | 0.530 | 0.531 | 0.449 | 0.426 | 0.525 | 0.435 | **0.406** |
| 2 | 0.194 | 0.196 | 0.195 | 0.256 | 0.255 | 0.255 | 0.204 | **0.188** | 0.262 | 0.194 | 0.317 |
| 3 | **0.190** | 0.195 | 0.195 | 0.369 | 0.377 | 0.369 | 0.217 | 0.208 | 0.365 | 0.195 | 0.355 |
| 4 | 0.810 | 0.769 | 0.735 | 0.784 | 0.754 | 0.736 | 0.814 | 1.219 | 0.833 | 0.853 | **0.157** |
| 5 | 0.141 | 0.142 | 0.142 | 0.145 | 0.145 | 0.145 | **0.131** | 0.142 | 0.145 | 0.140 | 0.184 |
| 6 | 0.465 | 0.464 | 0.464 | 0.661 | 0.661 | 0.660 | 0.903 | 0.651 | 0.660 | **0.464** | 0.489 |
| 7 | 0.162 | 0.161 | 0.161 | 0.282 | 0.282 | 0.282 | 0.620 | 0.537 | 0.282 | 0.162 | **0.083** |
| 8 | 0.142 | 0.142 | 0.142 | 0.337 | 0.336 | 0.336 | 0.986 | 0.907 | 0.336 | 0.142 | **0.079** |
| 9 | 0.100 | 0.100 | 0.100 | 0.336 | 0.335 | 0.335 | 0.738 | 0.678 | 0.335 | 0.100 | **0.066** |
| 10 | 0.217 | 0.215 | 0.215 | 0.363 | 0.363 | 0.363 | 0.780 | 0.700 | 0.363 | 0.216 | **0.159** |

**True positive rate in the sign-conflict models.** 3 of the 12 nonzero coefficients are
negative (3 of 9 in model 9).

| model | zt | zt_std | uni | uni_loo | lasso | ada_uni | ada_lasso |
|---|---:|---:|---:|---:|---:|---:|---:|
| 6 | 0.940 | 0.902 | 0.749 | 0.746 | 0.907 | 0.942 | 0.850 |
| 7 | 0.995 | 1.000 | 0.750 | 0.750 | 1.000 | 0.995 | 0.999 |
| 8 | 1.000 | 1.000 | 0.747 | 0.744 | 1.000 | 1.000 | 1.000 |
| 9 | 1.000 | 1.000 | 0.667 | 0.667 | 1.000 | 1.000 | 1.000 |
| 10 | 0.977 | 0.950 | 0.530 | 0.514 | 0.951 | 0.978 | 0.907 |

**Findings.**

1. **Model 4 is fixed.** zeta_tau went from 5.1 in `linear_fixed` to 0.81 (`standardize = FALSE`)
   and 0.78 (`TRUE`). sum_rate is the best of the three parameterizations here (0.735).
   adalasso_lasso (0.157) is far ahead of every other method in this design.
2. **`standardize = FALSE` turns suni2D into adalasso_uni; `standardize = TRUE` turns it into the
   lasso.**
   - With `FALSE`, zeta_tau is within 0.004 of adalasso_uni in every model except model 4.
     There it is 0.04 better, and sum_rate is 0.12 better.
   - With `TRUE`, suni2D is within 0.015 of the lasso everywhere except model 4.
   - `FALSE` beats `TRUE` by 0.06-0.24 in models 1-3 and 6-10. Model 5 is a tie (-0.003,
     se 0.003). In model 4 `TRUE` wins by 0.026 (paired differences, zeta_tau).
   - The three parameterizations agree to within 0.01, except in model 4.
3. **The sign penalty helps only in model 4.** suni2D with `standardize = FALSE` matches the
   adaptive lasso with univariate weights almost exactly. This suggests CV keeps the sign
   penalty weak, even in the sign-conflict models. The chosen tau is not saved, so this is an
   inference. In models 6-10 the weak penalty is what lets suni2D keep the negative
   coefficients that uniLasso cannot.
4. **The suni2D `_std` rows are slightly worse than the suni2D rows of `linear_fixed` in models 1-3.**
   - The numbers: 0.531 vs 0.506, 0.256 vs 0.232 and 0.369 vs 0.350.
   - That run used the old solver, whose early stops act as extra shrinkage (see
     reason 3 above). It also used a shorter grid.
5. **uniLasso: LOO vs no LOO** (paired differences, no LOO minus LOO).
   - No LOO is better in model 4 (-0.41) and model 5 (-0.011).
   - It is slightly worse in models 1-2 (+0.02) and about equal in model 3.
   - It is clearly worse in the sign-conflict models (+0.06 to +0.25).
6. **Sign conflict.** Both uniLasso fits miss every negative coefficient in models 6-9. Their
   tpr equals the share of positive coefficients. In model 10 (AR(1)) they also lose about a
   third of the positive ones.
   - The two uniLasso fits have the two worst rpe values in models 7-10, 3-7 times that of
     suni2D.
   - In model 6, uniLasso without LOO is the worst (0.90). uniLasso_loo (0.65) is level with
     the lasso, 1.4 times suni2D.
   - The other methods find 85-100% of the signals. adalasso_lasso is the most accurate in
     models 7-10, and suni2D / adalasso_uni in model 6.
7. **Baselines reproduce.** In models 1-5, the lasso, uniLasso_loo, adalasso_uni and
   adalasso_lasso rows equal those of `linear_fixed` exactly, in every replicate and criterion
   (same data, folds and fits).
