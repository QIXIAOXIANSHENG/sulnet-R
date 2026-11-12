## Summary of the model 1 diagnostic (model1_cv_diagnostic.R); prints the
## numbers quoted in "Model 1: why suni2D trails uniLasso and adalasso_uni" in
## ../README_tuned.md.
## Usage, from the package root:
##   Rscript inst/simulation/binomial_tuned/diagnostics/model1_cv_analysis.R [file]
## Default file: inst/simulation/binomial_tuned/diagnostics/model1_cv_diagnostic.rds
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1 && nzchar(args[1])) args[1] else
  "inst/simulation/binomial_tuned/diagnostics/model1_cv_diagnostic.rds"
R <- readRDS(f)
g <- function(fun) sapply(R, fun)

## test rpe_kl at the CV choice (and the path oracle) per seed
K <- data.frame(
  pipe = g(function(r) r$pipe$kl_na), pipe_pkg = g(function(r) r$pipe$kl_pkg), pipe_ff = g(function(r) r$pipe$kl_ff),
  s_inf = g(function(r) r$s_inf$kl_na), s_inf_pkg = g(function(r) r$s_inf$kl_pkg), s_inf_ff = g(function(r) r$s_inf$kl_ff),
  s_0 = g(function(r) r$s_0$kl_na), s_0_ff = g(function(r) r$s_0$kl_ff),
  uni = g(function(r) r$uni_kl), uni_h = g(function(r) r$uni_h$kl),
  uni_conv = g(function(r) r$cc_kl), uni_conv_h = g(function(r) r$cc_h$kl),
  ada = g(function(r) r$ada_kl), ada_h = g(function(r) r$ada_h$kl),
  s_inf_orc = g(function(r) r$s_inf$kl_orc), s_0_orc = g(function(r) r$s_0$kl_orc), pipe_orc = g(function(r) r$pipe$kl_orc),
  uni_orc = g(function(r) r$uni_orc), uni_conv_orc = g(function(r) r$cc_orc), ada_orc = g(function(r) r$ada_orc)
)
ms <- function(v) sprintf("%.3f (%.3f)", mean(v), sd(v) / sqrt(length(v)))
cat("seeds:", nrow(K), "\nmean rpe_kl (se):\n")
print(noquote(sapply(K, ms)))
pd <- function(a, b) {
  d <- K[[a]] - K[[b]]
  sprintf("%-12s - %-12s %+.3f (%.3f)  share<0 %.2f  share==0 %.2f", a, b, mean(d),
    sd(d) / sqrt(length(d)), mean(d < -1e-10), mean(abs(d) < 1e-10))
}
cat("\npaired differences, mean (se):\n")
for (pr in list(
  c("pipe", "uni"), c("pipe", "ada"),
  c("uni_conv", "uni"), c("uni_conv", "ada"), c("uni_conv_h", "uni_conv"), c("s_inf", "uni_conv_h"),
  c("s_inf", "uni_h"), c("uni_h", "uni"), c("s_0", "ada_h"), c("ada_h", "ada"),
  c("pipe", "s_inf"), c("s_0", "s_inf"),
  c("s_inf_orc", "uni_conv_orc"), c("s_0_orc", "ada_orc"), c("s_inf_orc", "uni_orc"), c("pipe_orc", "s_inf_orc"),
  c("pipe_ff", "pipe"), c("s_inf_ff", "s_inf"), c("pipe_pkg", "pipe"), c("s_inf_pkg", "s_inf")
)) cat(pd(pr[1], pr[2]), "\n")

cat("\nshare of seeds where a fold path stopped before the full path (s_inf):",
  mean(g(function(r) any(r$s_inf$last_fold < r$s_inf$n_full))), "\n")
cat("share where the full path stopped early (s_inf):", mean(g(function(r) r$s_inf$n_full < 100)), "\n")
cat("pipeline: chosen tau / zmax:\n")
print(table(signif(g(function(r) r$pipe$tau_na / r$zmax), 2)))
cat("s_0 and s_inf select the same cell:", mean(g(function(r) all(r$s_0$sel_na == r$s_inf$sel_na))), "\n")

## selected penalties (zeta and lambda are on the same scale)
sel <- data.frame(
  s_inf = g(function(r) r$s_inf$zeta_na),
  uni = g(function(r) r$uni_lam[r$uni_jmin]),
  uni_h = g(function(r) r$uni_h$lam),
  uni_conv = g(function(r) r$cc_lam[r$cc_jmin]),
  uni_conv_h = g(function(r) r$cc_h$lam)
)
cat("\nmedian selected penalty:\n")
print(round(sapply(sel, median), 4))
cat("share of seeds where the honest choice is larger (uni_conv_h vs uni_conv):",
  mean(sel$uni_conv_h > sel$uni_conv), "\n")
orc_lam <- g(function(r) r$cc_lam[which.min(r$cc_path)])
cat("median test-optimal penalty on the uni_conv path:", round(median(orc_lam), 4), "\n")
cat("support at the CV choice: s_inf", mean(g(function(r) r$s_inf$df_na)),
  " uni", mean(g(function(r) r$uni_df)), " ada", mean(g(function(r) r$ada_df)), "\n")
cat("mean CV deviance at its minimum: uni_conv (reused fits)", round(mean(g(function(r) min(r$cc_cvm))), 3),
  " uni_conv_h (honest)", round(mean(g(function(r) min(r$cc_h$cvm))), 3),
  " s_inf", round(mean(g(function(r) min(r$s_inf$cvm, na.rm = TRUE))), 3), "\n")

## first stage: uniLasso::uniInfo (2 Newton steps) vs sulnet::uniFit (converged)
B <- t(g(function(r) r$bcmp))
cat("\nuniInfo / uniFit slope ratio: median over seeds of the per-seed median, active",
  round(median(B[, "ratio_active"]), 3), ", noise", round(median(B[, "ratio_noise"]), 3),
  "; share of seeds with some ratio outside [0.8, 1.25]:",
  mean(B[, "ratio_min"] < 0.8 | B[, "ratio_max"] > 1.25), "\n")
