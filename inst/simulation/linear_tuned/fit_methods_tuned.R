## Methods of the tuned linear simulation (inst/simulation/linear_tuned).
##
## fit_methods() of compare_linear.R with two changes:
##   * suni2D is fitted with standardize = FALSE (rows suni2D_<method>) and
##     with standardize = TRUE (rows suni2D_<method>_std), both with
##     lambda.factor = SUNI_LMR, so that the first-coordinate grid reaches
##     1e-4 of its largest value (the default 0.01 for n < p is too short),
##     and with the default eps.
##   * uniLasso is fitted without LOO (row uniLasso) and with LOO
##     (row uniLasso_loo, the uniLasso row of inst/simulation/linear_fixed).
## lasso, adalasso_uni and adalasso_lasso are those of compare_linear.R.
##
## Needs compare_linear.R to be sourced first (ada_cv, GLMNET_LMR) and sulnet
## with the 2026-11-12 Gaussian solver fix (exact coordinate steps and a
## stopping rule weighted by the column mean squares); with the earlier
## solver, fits with standardize = FALSE stop far from convergence.
## Draws no random numbers.

SUNI_LMR <- 1e-4
method_names_tuned <- c(
  "suni2D_zeta_tau", "suni2D_lambda_alpha", "suni2D_sum_rate",
  "suni2D_zeta_tau_std", "suni2D_lambda_alpha_std", "suni2D_sum_rate_std",
  "uniLasso", "uniLasso_loo", "lasso", "adalasso_uni", "adalasso_lasso"
)

fit_methods_tuned <- function(x, y, newx, foldid, s = c("min", "1se")) {
  s <- match.arg(s)
  s_glmnet <- paste0("lambda.", s)
  s_suni <- paste0("cv.", s)
  p <- ncol(x)
  ctrl <- glmnet.control()
  glmnet.control(fdev = 0, devmax = 1)
  on.exit(glmnet.control(fdev = ctrl$fdev, devmax = ctrl$devmax), add = TRUE)
  na_fit <- list(coef = rep(NA_real_, p), pred = rep(NA_real_, nrow(newx)))
  safely <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      warning(name, " failed: ", conditionMessage(e), call. = FALSE)
      na_fit
    })
  }
  suni_out <- function(method, std) {
    cvfit <- cv.suni2D(x, y, foldid = foldid, method = method,
      standardize = std, lambda.factor = SUNI_LMR)
    list(
      coef = as.numeric(coef(cvfit, s = s_suni)[[1]])[-1],
      pred = as.numeric(predict(cvfit, newx, s = s_suni)[[1]])
    )
  }
  glmnet_out <- function(cvfit) {
    list(
      coef = as.numeric(coef(cvfit, s = s_glmnet))[-1],
      pred = as.numeric(predict(cvfit, newx, s = s_glmnet))
    )
  }
  ## first stages of the adaptive lassos, as in compare_linear.R
  first_uni <- function(x, y, foldid) drop(cov(x, y)) / apply(x, 2, var)
  first_lasso <- function(x, y, foldid) {
    as.numeric(coef(
      cv.glmnet(x, y, foldid = foldid, lambda.min.ratio = GLMNET_LMR),
      s = s_glmnet
    ))[-1]
  }

  fits <- list()
  for (std in c(FALSE, TRUE)) {
    for (method in c("zeta_tau", "lambda_alpha", "sum_rate")) {
      name <- paste0("suni2D_", method, if (std) "_std")
      fits[[name]] <- safely(name, suni_out(method, std))
    }
  }
  fits$uniLasso <- safely("uniLasso", glmnet_out(cv.uniLasso(x, y,
    foldid = foldid, loo = FALSE, lambda.min.ratio = GLMNET_LMR)))
  fits$uniLasso_loo <- safely("uniLasso_loo", glmnet_out(cv.uniLasso(x, y,
    foldid = foldid, loo = TRUE, lambda.min.ratio = GLMNET_LMR)))
  fits$lasso <- safely("lasso", glmnet_out(
    cv.glmnet(x, y, foldid = foldid, lambda.min.ratio = GLMNET_LMR)
  ))
  fits$adalasso_uni <- safely("adalasso_uni", ada_cv(
    x, y, newx, foldid, first_uni(x, y, foldid), first_uni, s
  ))
  ## the full-data first stage is the lasso fit above
  fits$adalasso_lasso <- safely("adalasso_lasso", if (anyNA(fits$lasso$coef)) {
    stop("the lasso fit failed")
  } else {
    ada_cv(x, y, newx, foldid, fits$lasso$coef, first_lasso, s)
  })

  fits <- fits[method_names_tuned]
  list(
    coef = do.call(cbind, lapply(fits, `[[`, "coef")),
    pred = do.call(cbind, lapply(fits, `[[`, "pred"))
  )
}
