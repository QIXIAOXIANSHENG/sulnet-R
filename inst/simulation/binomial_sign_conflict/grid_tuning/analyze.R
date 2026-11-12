## Post-hoc evaluation of second-coordinate sub-grids from the per-cell CV and
## test criteria saved by round1.R, round2.R and round3.R. Usage, from the
## package root:
##   Rscript inst/simulation/binomial_sign_conflict/grid_tuning/analyze.R \
##     [round] [res_dir] [out_prefix]
## round names the results folder (round1, round2_zeta, round2_dense, round3);
## defaults: round1, <TUNE_DIR>/results/<round>, <TUNE_DIR>/analysis/<round>.
## Writes <out_prefix>_analysis.rds and _summary/_paired/_oracle/_picks.csv.
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
tune_dir <- Sys.getenv("TUNE_DIR", file.path(Sys.getenv("SIM_DIR", "inst/simulation"),
  "binomial_sign_conflict", "grid_tuning"))
round <- arg(1, "round1")
res_dir <- arg(2, file.path(tune_dir, "results", round))
out_prefix <- arg(3, file.path(tune_dir, "analysis", round))
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(sulnet))
criteria <- c("sse", "rpe_class", "rpe_kl", "support", "sign", "tpr", "fpr")

files <- list.files(res_dir, pattern = "^m\\d+_b\\d+\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
ok <- vapply(res, function(r) is.null(r$error), logical(1))
if (any(!ok)) {
  message("errors in: ", paste(basename(files[!ok]), collapse = ", "))
  for (r in res[!ok]) message("  model ", r$model, " rep ", r$b, ": ", r$error)
}
res <- res[ok]
models <- vapply(res, `[[`, numeric(1), "model")
message(length(res), " results; per model: ", paste(table(models), collapse = " "))
warns <- unique(unlist(lapply(res, `[[`, "warnings")))
if (length(warns)) message("distinct warnings (first 10):\n  ", paste(head(warns, 10), collapse = "\n  "))

near <- function(a, b) abs(a - b) <= 1e-7 * pmax(1, abs(b))
match_tol <- function(v, grid) vapply(v, function(x) {
  i <- which(near(grid, x))
  if (length(i)) i[1] else NA_integer_
}, integer(1))

## criteria of the cell selected by getmin2D on rows/cols of a saved config
select_sub <- function(s, rows, cols, loss, rule) {
  if (is.null(s) || !length(rows) || !length(cols)) return(rep(NA_real_, 7))
  m <- s$cvm[[loss]][rows, cols, drop = FALSE]
  sd <- s$cvsd[[loss]][rows, cols, drop = FALSE]
  if (all(is.na(m))) return(rep(NA_real_, 7))
  lm <- sulnet:::getmin2D(s$pen1[cols], s$pen2[rows], m, sd, s$method)[[rule]]
  i2 <- rows[match(lm[[2]], s$pen2[rows])]
  i1 <- cols[match(lm[[1]], s$pen1[cols])]
  s$crit[i2, i1, ]
}

## ---- candidate sub-grids, defined on relative second-coordinate values ----
## zeta_tau: tau / zmax ; lambda_alpha: alpha ; sum_rate: rate
lg <- function(from, to, by) 10^seq(from, to, by = by)
zt3 <- list(
  all = NULL, inf = 1e3, zero = 0, inf_0 = c(1e3, 0),
  interior = lg(1, -4, -0.25), # no endpoints
  inf_d0_3_0 = c(1e3, lg(0, -3, -0.25), 0), # the final-run grid
  inf_d0_2_0 = c(1e3, lg(0, -2, -0.25), 0)
)
cand <- list(
  zt_F = zt3, zt_T = zt3,
  zt_dense = list(
    all = NULL,
    inf = 1e3, zero = 0, inf_0 = c(1e3, 0),
    d0_2 = lg(0, -2, -0.125), # ~ default range, no endpoints
    inf_d0_2_0 = c(1e3, lg(0, -2, -0.125), 0),
    inf_d0_3_0 = c(1e3, lg(0, -3, -0.125), 0),
    inf_d1_2_0 = c(1e3, lg(1, -2, -0.125), 0),
    inf_d0_1_0 = c(1e3, lg(0, -1, -0.125), 0),
    inf_dm1_3_0 = c(1e3, lg(-1, -3, -0.125), 0),
    inf_c0_3_0 = c(1e3, lg(0, -3, -0.5), 0),
    inf_c0_2_0 = c(1e3, 1, 0.1, 0.01, 0),
    inf_c0_1_0 = c(1e3, 1, lg(-0.5, -1, -0.5), 0),
    inf_1_0 = c(1e3, 1, 0),
    inf_d0_2 = c(1e3, lg(0, -2, -0.125)),
    d0_2_0 = c(lg(0, -2, -0.125), 0),
    inf_d0_4_0 = c(1e3, lg(0, -4, -0.125), 0),
    inf_d0_5_0 = c(1e3, lg(0, -5, -0.125), 0),
    inf_c0_5_0 = c(1e3, lg(0, -5, -0.5), 0),
    inf_c0_5_0_1 = c(1e3, lg(0, -5, -1), 0),
    d0_5_0 = c(lg(0, -5, -0.125), 0),
    inf_dm2_5_0 = c(1e3, lg(-2, -5, -0.125), 0)
  ),
  la_dense = list(
    all = NULL,
    lin = seq(0.5, 1 / 1.01, length.out = 20),
    log3 = 1 - 0.5 * lg(0, -3, -0.125),
    log2 = 1 - 0.5 * lg(0, -2, -0.125),
    log3c = 1 - 0.5 * lg(0, -3, -0.5),
    ends = c(0.5, 1 - 0.5e-3),
    ends2 = c(0.5, 1 - 0.5e-2),
    half = 0.5, top = 1 - 0.5e-3,
    lin_top = c(seq(0.5, 1 / 1.01, length.out = 20), 1 - 0.5e-3),
    half_near = c(0.5, 0.5 + 0.5 * lg(-1, -4, -0.5)),
    lin_near = c(seq(0.5, 1 / 1.01, length.out = 20), 0.5 + 0.5 * lg(-1, -4, -0.5))
  ),
  sr_dense = list(
    all = NULL,
    def = 100^((0:19) / 19),
    log3 = lg(0, 3, 0.125),
    log2 = lg(0, 2, 0.125),
    log3c = lg(0, 3, 0.5),
    ends = c(1, 1e3),
    ends2 = c(1, 1e2),
    one = 1, top = 1e3,
    def_top = c(100^((0:19) / 19), 1e3),
    one_near = c(1, 1 + 2 * lg(-1, -4, -0.5)),
    def_near = c(100^((0:19) / 19), 1 + 2 * lg(-1, -4, -0.5))
  )
)

## zmax of a zeta_tau fit: its largest (default-grid) zeta, so that tau is
## relative to the fit's own scale (standardize TRUE and FALSE differ)
zmax_fit <- function(s, zmax) if (is.finite(s$zeta1)) s$zeta1 / (1 + 1e-9) else zmax
rel2 <- function(s, zmax) if (s$method == "zeta_tau") s$pen2 / zmax_fit(s, zmax) else s$pen2
rows_for <- function(s, rel, zmax) {
  if (is.null(rel)) return(seq_along(s$pen2))
  v <- rel2(s, zmax)
  i <- match_tol(rel, v)
  if (anyNA(i)) return(NULL) # candidate not contained in this grid
  sort(unique(i))
}

## ---- collect per-rep criteria ----
rows_out <- list()
add <- function(model, b, config, sub, cols, loss, rule, v) {
  rows_out[[length(rows_out) + 1L]] <<- c(
    list(model = model, b = b, config = config, sub = sub, cols = cols, loss = loss, rule = rule),
    as.list(setNames(as.numeric(v), criteria))
  )
}
pick_rows <- list() # distribution of selected second coordinate (full dense grid)
for (r in res) {
  for (mname in rownames(r$base)) {
    add(r$model, r$b, mname, "-", "-", "-", "-", r$base[mname, criteria])
  }
  for (cf in names(r$suni)) {
    s <- r$suni[[cf]]
    if (is.null(s)) next
    colsets <- list(all = seq_along(s$pen1))
    if (cf %in% c("la_dense", "sr_dense")) colsets$def <- tail(seq_along(s$pen1), 100)
    subs <- if (cf %in% names(cand)) cand[[cf]] else list(all = NULL)
    for (sn in names(subs)) {
      rws <- rows_for(s, subs[[sn]], r$zmax)
      if (is.null(rws)) next
      for (cn in names(colsets)) {
        for (loss in c("loss", "Brier")) {
          for (rule in c("cv.min", "cv.1se")) {
            add(r$model, r$b, cf, sn, cn, loss, rule,
                select_sub(s, rws, colsets[[cn]], loss, rule))
          }
        }
      }
    }
    for (loss in c("loss", "Brier")) {
      ij <- s$sel[[loss]]$cv.min
      if (!is.null(ij)) {
        pick_rows[[length(pick_rows) + 1L]] <- data.frame(
          model = r$model, b = r$b, config = cf, loss = loss,
          pen2 = rel2(s, r$zmax)[ij[["row"]]]
        )
      }
    }
  }
}
long <- do.call(rbind, lapply(rows_out, as.data.frame))
picks <- do.call(rbind, pick_rows)

## ---- oracle per second-coordinate level (dense configs): mean over reps of
## the best-over-first-coordinate test rpe_kl, and of the cvm-selected cell
## within that level ----
orc <- list()
for (r in res) for (cf in names(r$suni)) {
  s <- r$suni[[cf]]
  if (is.null(s)) next
  v <- rel2(s, r$zmax)
  kl <- s$crit[, , "rpe_kl"]
  best <- apply(kl, 1, function(z) if (all(is.na(z))) NA else min(z, na.rm = TRUE))
  cvsel <- vapply(seq_along(v), function(a) {
    m <- s$cvm$loss[a, ]
    if (all(is.na(m))) NA_real_ else kl[a, which.min(m)]
  }, numeric(1))
  orc[[length(orc) + 1L]] <- data.frame(model = r$model, b = r$b, config = cf,
    level = seq_along(v), pen2 = v, oracle_kl = best, cvlevel_kl = cvsel)
}
orc <- do.call(rbind, orc)

## ---- summaries ----
key <- c("model", "config", "sub", "cols", "loss", "rule")
agg <- function(df) {
  sp <- split(df, df[key], drop = TRUE)
  do.call(rbind, lapply(sp, function(d) {
    out <- d[1, key]
    out$n <- sum(!is.na(d$rpe_kl))
    for (cr in criteria) {
      out[[cr]] <- mean(d[[cr]], na.rm = TRUE)
      out[[paste0(cr, "_se")]] <- sd(d[[cr]], na.rm = TRUE) / sqrt(sum(!is.na(d[[cr]])))
    }
    out
  }))
}
summ <- agg(long)
rownames(summ) <- NULL
## paired differences in rpe_kl against uniLasso and adalasso_uni
base_kl <- long[long$config %in% c("uniLasso", "adalasso_uni", "lasso"), c("model", "b", "config", "rpe_kl")]
bw <- reshape(base_kl, idvar = c("model", "b"), timevar = "config", direction = "wide")
lj <- merge(long, bw, by = c("model", "b"))
lj$d_uni <- lj$rpe_kl - lj$rpe_kl.uniLasso
lj$d_ada <- lj$rpe_kl - lj$rpe_kl.adalasso_uni
lj$d_best <- lj$rpe_kl - pmin(lj$rpe_kl.uniLasso, lj$rpe_kl.adalasso_uni)
pd <- aggregate(cbind(d_uni, d_ada, d_best) ~ model + config + sub + cols + loss + rule, lj,
  function(z) c(m = mean(z), se = sd(z) / sqrt(length(z))))
pd <- do.call(data.frame, pd)

saveRDS(list(long = long, summ = summ, picks = picks, oracle = orc, paired = pd),
  paste0(out_prefix, "_analysis.rds"))
write.csv(summ, paste0(out_prefix, "_summary.csv"), row.names = FALSE)
write.csv(pd, paste0(out_prefix, "_paired.csv"), row.names = FALSE)
osum <- aggregate(cbind(oracle_kl, cvlevel_kl) ~ model + config + level + pen2, orc, mean)
write.csv(osum, paste0(out_prefix, "_oracle.csv"), row.names = FALSE)
pk <- picks
pk$bin <- if (nrow(pk)) as.character(signif(pk$pen2, 3)) else character()
write.csv(as.data.frame(table(model = pk$model, config = pk$config, loss = pk$loss, pen2 = pk$bin)),
  paste0(out_prefix, "_picks.csv"), row.names = FALSE)
message("written ", out_prefix, "_*.csv")
