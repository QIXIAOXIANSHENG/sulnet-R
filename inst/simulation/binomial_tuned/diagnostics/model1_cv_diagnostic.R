## Model 1 diagnostic: why suni2D_zeta_tau trails uniLasso and adalasso_uni in
## model 1 although its tau grid contains both endpoints (tau = 1e3 zmax for
## uniLasso, tau = 0 for adalasso_uni). See "Model 1: why suni2D trails uniLasso
## and adalasso_uni" in ../README_tuned.md.
##
## For each seed, with model 1 data and folds exactly as in the pipeline, it fits:
##   pipe      cv.suni2D with the full tau grid (the pipeline fit)
##   s_inf     cv.suni2D with tau = 1e3 zmax only (the uniLasso end)
##   s_0       cv.suni2D with tau = 0 only (the adalasso_uni end)
##   *_ff      the same fold fits, but a fold whose path stopped early
##             contributes its last computed cell to the later cells, as
##             cv.glmnet does, instead of making the cell NA
##   uni       cv.uniLasso(loo = FALSE): univariate fits from uniLasso::uniInfo
##             (2 Newton steps), computed once on all the data
##   cc        ("uni_conv") cv.uniLasso(loo = FALSE) with fully converged
##             univariate fits (sulnet::uniFit, as suni2D and adalasso_uni use)
##   ada       adalasso_uni as in compare_binom.R
##   *_h       ("honest") uni, cc and ada with the univariate fits redone in
##             each training fold, as cv.suni2D does; full-data lambda grid
##   *_orc     the best test rpe_kl anywhere along the full-data path
##
## Usage, from the package root (the folder that contains inst/):
##   Rscript inst/simulation/binomial_tuned/diagnostics/model1_cv_diagnostic.R \
##     [seeds] [ncores] [out_dir]
## Defaults: 1:200, 1, inst/simulation/binomial_tuned/diagnostics/rerun.
## ncores > 1 uses a PSOCK cluster, which also works on Windows. Each seed is
## saved as <out_dir>/reps/b<bbb>.rds, and saved seeds are skipped on a rerun.
## At the end <out_dir>/model1_cv_diagnostic.rds holds the list of all saved
## seeds; model1_cv_analysis.R summarizes it. The shipped
## model1_cv_diagnostic.rds (seeds 1-200) was made on a Windows PC (R 4.5.x,
## sulnet b11feb6, glmnet 4.1-10, uniLasso 2.11) in about 3 minutes on 20 cores.
## Model 1 draws depend on the platform (see the README), so these seeds are
## not the same data as the Argon run.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
seeds <- eval(parse(text = arg(1, "1:200")))
ncores <- as.integer(arg(2, "1"))
out_dir <- arg(3, "inst/simulation/binomial_tuned/diagnostics/rerun")
repo <- normalizePath(".")
rep_dir <- file.path(out_dir, "reps")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)

worker_init <- function(repo) {
  suppressPackageStartupMessages({
    library(sulnet); library(glmnet); library(uniLasso); library(MASS)
  })
  source(file.path(repo, "inst/simulation/compare_binom.R"))
  invisible(TRUE)
}

one_seed <- function(b, out_dir) {
  f_out <- file.path(out_dir, sprintf("b%03d.rds", b))
  if (file.exists(f_out)) return("skip")
  ## ---- data and folds, as binomial_tuned_simulation.R + run_simulation_binom
  n <- 100; ntest <- 30; p <- 300; rho <- 0.8
  Sigma <- rho^abs(outer(seq_len(p), seq_len(p), "-"))
  gamma <- rep(0, p); gamma[1:30] <- 1
  X <- matrix(0, n, p)
  set.seed(b)
  ybase <- rbinom(n, 1, 0.5)
  for (i in 1:n) {
    mu <- rep(0, p); if (ybase[i] == 1) mu[1:20] <- 0.5
    X[i, ] <- mvrnorm(1, mu = mu, Sigma = Sigma)
  }
  xtest <- mvrnorm(ntest, rep(0, p), Sigma = Sigma)
  y <- as.numeric(rbinom(n, 1, plogis(X %*% gamma)))
  storage.mode(X) <- "double"
  foldid <- integer(n)
  foldid[order(y, runif(n))] <- rep_len(1:10, n)
  eta_new <- drop(xtest %*% gamma)
  crit <- function(coef, link) eval_criteria_binom(as.matrix(coef), as.matrix(link), gamma, eta_new)
  kl <- function(coef, link) crit(coef, link)[, "rpe_kl"]
  ys <- c(-1, 1)[y + 1]
  fmax <- log(1 / 1e-5 - 1)

  ## ---- suni2D: full fit, fold fits, CV with the NA rule and fill-forward
  f3 <- suni2D(X, y, family = "binomial", standardize = FALSE, nzeta = 3, ntau = 3)
  zmax <- f3$zeta[1] / (1 + 1e-9)
  TAU_REL <- c(1e3, 10^seq(0, -3, by = -0.25), 0)
  suni_cv <- function(tau) {
    full <- suni2D(X, y, family = "binomial", standardize = FALSE, tau = tau)
    zeta <- full$zeta
    folds <- lapply(1:10, function(k) {
      w <- foldid == k
      suni2D(X[!w, ], y[!w], family = "binomial", standardize = FALSE,
        zeta = zeta, tau = tau)
    })
    nt <- length(full$tau); nz <- length(zeta)
    ## held-out link predictions: list over tau of n x nz
    pm <- lapply(1:nt, function(a) matrix(NA_real_, n, nz))
    pm_ff <- pm
    last_fold <- matrix(NA_integer_, 10, nt)
    for (k in 1:10) {
      w <- foldid == k
      pr <- predict(folds[[k]], X[w, , drop = FALSE], type = "link")
      for (a in 1:nt) {
        m <- pr[[a]]
        pm[[a]][w, ] <- m
        ok <- which(!is.na(folds[[k]]$gamma0[, a]))
        last <- if (length(ok)) max(ok) else 0L
        last_fold[k, a] <- last
        mf <- m
        if (last > 0 && last < nz) mf[, (last + 1):nz] <- m[, last]
        pm_ff[[a]][w, ] <- mf
      }
    }
    cvm_of <- function(pml) t(sapply(pml, function(m) {
      m <- pmin(pmax(m, -fmax), fmax)
      colMeans(2 * log(1 + exp(-ys * m)))
    })) # nt x nz
    comp <- t(!is.na(full$gamma0)) # nt x nz, computed by the full fit
    cvm <- cvm_of(pm); cvm[!comp] <- NA
    cvm_ff <- cvm_of(pm_ff); cvm_ff[!comp] <- NA
    preds <- predict(full, xtest, type = "link")
    test_kl <- t(sapply(1:nt, function(a) {
      v <- kl(full$gamma[[a]], preds[[a]]); v[!comp[a, ]] <- NA; v
    }))
    if (nt == 1) { test_kl <- matrix(test_kl, 1) }
    pick <- function(cv) {
      if (all(is.na(cv))) return(c(a = NA, z = NA))
      ij <- arrayInd(which.min(cv), dim(cv)); c(a = ij[1], z = ij[2])
    }
    s_na <- pick(cvm); s_ff <- pick(cvm_ff)
    ## also the packaged cv.suni2D, to check the NA-rule replication
    cvfit <- cv.suni2D(X, y, foldid = foldid, family = "binomial",
      standardize = FALSE, tau = tau, pred.loss = "loss")
    pkg_kl <- kl(as.numeric(coef(cvfit, s = "cv.min")[[1]])[-1],
      as.numeric(predict(cvfit, xtest, s = "cv.min", type = "link")[[1]]))
    list(
      kl_na = test_kl[s_na["a"], s_na["z"]], kl_ff = test_kl[s_ff["a"], s_ff["z"]],
      kl_pkg = unname(pkg_kl), kl_orc = min(test_kl, na.rm = TRUE),
      sel_na = s_na, sel_ff = s_ff, zeta = zeta, tau = full$tau,
      n_full = rowSums(comp), n_cv = rowSums(!is.na(cvm)),
      last_fold = last_fold,
      zeta_na = zeta[s_na["z"]], zeta_ff = zeta[s_ff["z"]],
      tau_na = full$tau[s_na["a"]], tau_ff = full$tau[s_ff["a"]],
      test_kl = test_kl, cvm = cvm,
      df_na = full$df[s_na["z"], s_na["a"]]
    )
  }
  s_inf <- suni_cv(1e3 * zmax)
  s_0 <- suni_cv(0)
  pipe <- suni_cv(zmax * TAU_REL)

  ## ---- glmnet-based baselines, with full-data and honest univariate fits
  glm_path_kl <- function(gfit) {
    cf <- as.matrix(gfit$beta)
    lk <- predict(gfit, xtest, type = "link")
    kl(cf, lk)
  }
  uf <- uniFit(X, y, loo = FALSE, family = "binomial")
  ## uniLasso, as cv.uniLasso(loo = FALSE)
  cu <- cv.uniLasso(X, y, family = "binomial", foldid = foldid, loo = FALSE)
  uni_kl <- kl(as.numeric(coef(cu, s = "lambda.min"))[-1],
    as.numeric(predict(cu, xtest, s = "lambda.min", type = "link")))
  uni_path <- glm_path_kl(cu$glmnet.fit)
  ## adalasso_uni
  ca <- cv.glmnet(X, y, family = "binomial", foldid = foldid,
    penalty.factor = 1 / abs(uf$beta), standardize = FALSE)
  ada_kl <- kl(as.numeric(coef(ca, s = "lambda.min"))[-1],
    as.numeric(predict(ca, xtest, s = "lambda.min", type = "link")))
  ada_path <- glm_path_kl(ca$glmnet.fit)
  ## honest versions: univariate fits redone in each training fold; the
  ## full-data lambda grid; folds clamped beyond their path end as cv.glmnet
  honest <- function(type) {
    lam <- switch(type, uni = cu$lambda, uni_conv = cc$lambda, ada = ca$lambda)
    pm <- matrix(NA_real_, n, length(lam))
    for (k in 1:10) {
      w <- foldid == k
      ## uniLasso features come from uniLasso::uniInfo (as in cv.uniLasso);
      ## the adalasso weights from sulnet::uniFit (as in compare_binom.R)
      ufk <- if (type == "uni") uniInfo(X[!w, ], y[!w], "binomial", loo = FALSE) else {
        u <- uniFit(X[!w, ], y[!w], loo = FALSE, family = "binomial")
        list(beta = u$beta, beta0 = u$b0)
      }
      if (type %in% c("uni", "uni_conv")) {
        Ftr <- sweep(sweep(X[!w, ], 2, ufk$beta, "*"), 2, ufk$beta0, "+")
        Fte <- sweep(sweep(X[w, , drop = FALSE], 2, ufk$beta, "*"), 2, ufk$beta0, "+")
        g <- glmnet(Ftr, y[!w], family = "binomial", lower.limits = 0,
          standardize = FALSE, lambda = lam)
        pm[w, ] <- predict(g, Fte, s = lam, type = "link")
      } else {
        g <- glmnet(X[!w, ], y[!w], family = "binomial",
          penalty.factor = 1 / abs(ufk$beta), standardize = FALSE, lambda = lam)
        pm[w, ] <- predict(g, X[w, , drop = FALSE], s = lam, type = "link")
      }
    }
    pm <- pmin(pmax(pm, -fmax), fmax)
    cvm <- colMeans(2 * log(1 + exp(-ys * pm)))
    j <- which.min(cvm)
    path <- switch(type, uni = uni_path, uni_conv = cc_path, ada = ada_path)
    list(kl = unname(path[j]), j = j, lam = lam[j], cvm = cvm)
  }
  ## uniLasso with fully converged univariate fits (sulnet::uniFit), as suni2D
  cc <- cv.uniLasso(X, y, family = "binomial", foldid = foldid, loo = FALSE,
    info = list(beta = uf$beta, beta0 = uf$b0))
  cc_kl <- kl(as.numeric(coef(cc, s = "lambda.min"))[-1],
    as.numeric(predict(cc, xtest, s = "lambda.min", type = "link")))
  cc_path <- glm_path_kl(cc$glmnet.fit)
  ui <- uniInfo(X, y, "binomial", loo = FALSE)
  nz <- uf$beta != 0
  bratio <- ui$beta / uf$beta
  bcmp <- c(ratio_active = median(bratio[1:30]), ratio_noise = median(bratio[31:300]),
    ratio_min = min(bratio), ratio_max = max(bratio),
    cor = cor(ui$beta, uf$beta))
  uni_h <- honest("uni")
  cc_h <- honest("uni_conv")
  ada_h <- honest("ada")

  res <- list(b = b, zmax = zmax, s_inf = s_inf, s_0 = s_0, pipe = pipe,
    uni_kl = unname(uni_kl), uni_orc = min(uni_path), uni_lam = cu$lambda,
    uni_jmin = match(cu$lambda.min, cu$lambda),
    ada_kl = unname(ada_kl), ada_orc = min(ada_path), ada_lam = ca$lambda,
    ada_jmin = match(ca$lambda.min, ca$lambda),
    uni_h = uni_h, ada_h = ada_h, cc_h = cc_h,
    cc_kl = unname(cc_kl), cc_orc = min(cc_path), cc_lam = cc$lambda,
    cc_jmin = match(cc$lambda.min, cc$lambda), bcmp = bcmp,
    uni_path = uni_path, uni_cvm = cu$cvm, cc_path = cc_path, cc_cvm = cc$cvm,
    ada_path = ada_path, ada_cvm = ca$cvm,
    uni_df = cu$nzero[match(cu$lambda.min, cu$lambda)],
    ada_df = ca$nzero[match(ca$lambda.min, ca$lambda)])
  saveRDS(res, f_out)
  "ok"
}

if (ncores > 1) {
  cl <- parallel::makeCluster(ncores)
  parallel::clusterCall(cl, worker_init, repo)
  out <- parallel::parLapplyLB(cl, seeds, one_seed, out_dir = rep_dir)
  parallel::stopCluster(cl)
} else {
  worker_init(repo)
  out <- lapply(seeds, one_seed, out_dir = rep_dir)
}
print(table(unlist(out)))

## all saved seeds in one file
fs <- list.files(rep_dir, "^b[0-9]+[.]rds$", full.names = TRUE)
saveRDS(lapply(fs, readRDS), file.path(out_dir, "model1_cv_diagnostic.rds"))
message("written ", file.path(out_dir, "model1_cv_diagnostic.rds"), " (", length(fs), " seeds)")
