## Boxplots of the real-data comparisons saved by
## inst/linear_true_data/true_dataset.R and
## inst/binomial_true_data/true_dataset.R.
##
## Usage:
##   source("inst/simulation/plot_real_data.R")
##   res_all <- read_real_results("./inst/binomial_true_data/results")
##   plot_real_results(res_all, file = "boxplots.pdf")
##
## Each saved result has an element raw: a list of per-split M x 2 matrices
## with methods in rows and columns (mse or misclass, support).

## Read every .rds result in dir into a list named by dataset. Other files in
## dir (summary.csv, raw.csv, plots) are ignored.
read_real_results <- function(dir) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  names(files) <- sub("\\.rds$", "", basename(files))
  lapply(files, readRDS)
}

## Short x-axis labels for the seven methods. Methods not listed keep their
## names.
method_labels_real <- c(
  suni2D_zeta_tau = "suni2D zt", suni2D_lambda_alpha = "suni2D la",
  suni2D_sum_rate = "suni2D sr", uniLasso = "uniLasso", lasso = "lasso",
  adalasso_uni = "ada uni", adalasso_lasso = "ada lasso"
)

## Range for one axis. With outline = FALSE the outliers are not drawn, so the
## range only covers the whiskers.
box_range <- function(v, outline) {
  if (!outline) v <- boxplot(as.data.frame(v), plot = FALSE)$stats
  v <- v[is.finite(v)]
  if (length(v) == 0) c(0, 1) else range(v)
}

## One panel: for each method, the error box (left axis) at method - 0.2 and
## the support box (right axis) at method + 0.2. The two axes are scaled
## separately for every dataset.
plot_real_panel <- function(
  err, supp, main, sub = NULL, labels, col, outline = TRUE
) {
  M <- ncol(err)
  xlim <- c(0.5, M + 0.5)
  ink <- "#52514e"
  draw <- function(v, at, colr, side) {
    plot.window(xlim = xlim, ylim = box_range(v, outline))
    boxplot(as.data.frame(v),
      at = at, add = TRUE, axes = FALSE, boxwex = 0.35, staplewex = 0.5,
      col = adjustcolor(colr, alpha.f = 0.3), border = colr, medlwd = 2,
      whisklty = 1, outline = outline, outpch = 16, outcex = 0.4
    )
    axis(side, col = ink, col.ticks = colr, col.axis = colr, lwd = 0,
      lwd.ticks = 1)
  }
  plot.new()
  draw(err, seq_len(M) - 0.2, col[1], 2)
  abline(v = seq_len(M - 1) + 0.5, col = "#e1e0d9")
  ## Methods that never finished have no boxes.
  failed <- colSums(!is.na(err)) == 0
  if (any(failed)) {
    text(which(failed), mean(par("usr")[3:4]), "failed",
      srt = 90, col = ink, cex = 0.8
    )
  }
  par(new = TRUE)
  plot.new()
  draw(supp, seq_len(M) + 0.2, col[2], 4)
  axis(1,
    at = seq_len(M), labels = labels, las = 2, col = ink, col.axis = ink,
    lwd = 0, lwd.ticks = 1
  )
  box(col = ink)
  title(main = main, line = if (is.null(sub)) 0.6 else 1.2, cex.main = 1)
  if (!is.null(sub)) mtext(sub, side = 3, line = 0.3, cex = 0.6, col = ink)
}

## Draw one panel per dataset in a single figure, ncol panels per row.
## res     a list of saved results, as returned by read_real_results()
## file    NULL draws on the current device; a path ending in .pdf or .png
##         opens that device, sized from the number of panels, and closes it
## methods subset and order of the methods; NULL keeps all, in saved order
## labels  named vector of x-axis labels, matched to the method names
## outline FALSE hides the outliers and zooms both axes to the whiskers,
##         which helps when a few splits have a much larger mse
## col     colors of the error boxes and left axis, and of the support boxes
##         and right axis
## The error measure (mse or misclass) is the raw column other than support.
## Returns the file path (or NULL) invisibly.
plot_real_results <- function(
  res, file = NULL, ncol = 5, methods = NULL, labels = method_labels_real,
  outline = TRUE, col = c("#2a78d6", "#eb6834"),
  panel_width = 3.2, panel_height = 3
) {
  stopifnot(length(res) > 0, !is.null(names(res)))
  raw1 <- res[[1]]$raw[[1]]
  err_name <- setdiff(colnames(raw1), "support")[1]
  err_label <- switch(err_name,
    mse = "MSE",
    misclass = "Misclassification rate",
    err_name
  )
  if (is.null(methods)) methods <- rownames(raw1)
  lab <- ifelse(methods %in% names(labels), labels[methods], methods)

  nrow <- ceiling(length(res) / ncol)
  if (!is.null(file)) {
    width <- ncol * panel_width
    height <- nrow * panel_height + 0.5
    if (grepl("\\.pdf$", file, ignore.case = TRUE)) {
      pdf(file, width = width, height = height)
    } else if (grepl("\\.png$", file, ignore.case = TRUE)) {
      png(file, width = width, height = height, units = "in", res = 200)
    } else {
      stop("file must end in .pdf or .png")
    }
    on.exit(dev.off(), add = TRUE)
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE, after = FALSE)
  par(
    mfrow = c(nrow, ncol), oma = c(0, 0, 2, 0), mar = c(5, 3, 2.4, 3),
    mgp = c(2, 0.5, 0), tcl = -0.3
  )
  par(cex = 0.7)

  for (name in names(res)) {
    raw <- res[[name]]$raw
    ## nsplit x M matrices, one row per split
    err <- do.call(rbind, lapply(raw, function(r) r[methods, err_name]))
    supp <- do.call(rbind, lapply(raw, function(r) r[methods, "support"]))
    colnames(err) <- colnames(supp) <- methods
    info <- res[[name]]$info
    sub <- if (!is.null(info)) {
      sprintf("n = %d, p = %d", info[["n"]], info[["p"]])
    }
    plot_real_panel(err, supp, name, sub, lab, col, outline)
  }

  ## Legend across the top outer margin.
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
    new = TRUE)
  plot.new()
  legend("top",
    legend = c(paste(err_label, "(left axis)"), "Support (right axis)"),
    fill = adjustcolor(col, alpha.f = 0.3), border = col, text.col = col,
    horiz = TRUE, bty = "n", cex = 1.1
  )
  invisible(file)
}
