## Summary tables of the simulation results saved by
## inst/simulation/linear simulation.R and
## inst/simulation/binomial simulation.R, written to one pdf per family.
##
## Usage:
##   source("inst/simulation/table_simulation.R")
##   res <- read_sim_results("./inst/simulation/linear")
##   write_sim_tables(res, "tables.pdf", main = "Linear simulations")
##
## Each saved result has an element table: a methods x criteria matrix of
## "mean (se)" strings, as built by run_simulation(), run_simulation_binom()
## or summarize_sim(). Needs the gridExtra and gtable packages.

## Read every res_model<k>.rds in dir into a list named "Model <k>", ordered
## by k. Other files in dir (such as the table pdf) are ignored.
read_sim_results <- function(dir) {
  files <- list.files(dir, pattern = "^res_model[0-9]+\\.rds$",
    full.names = TRUE)
  if (length(files) == 0) stop("no res_model<k>.rds files in ", dir)
  k <- as.integer(gsub("[^0-9]", "", basename(files)))
  files <- files[order(k)]
  names(files) <- paste("Model", sort(k))
  lapply(files, readRDS)
}

## One table as a gtable: a bold title and an optional grey caption above a
## booktabs-style table (rules above and below the header and at the bottom,
## light stripes on every other row, numbers right-aligned).
sim_table_grob <- function(tab, title, caption = NULL, fontsize = 9) {
  tab <- unclass(tab)
  theme <- gridExtra::ttheme_minimal(
    base_size = fontsize,
    padding = grid::unit(c(4, 2.4), "mm"),
    core = list(
      fg_params = list(hjust = 1, x = 0.95),
      bg_params = list(fill = c("white", "grey95"), col = NA)
    ),
    colhead = list(fg_params = list(hjust = 1, x = 0.95, fontface = 2)),
    rowhead = list(fg_params = list(hjust = 0, x = 0.05, fontface = 1))
  )
  g <- gridExtra::tableGrob(tab, theme = theme)
  rule <- function(y, lwd) {
    grid::segmentsGrob(0, y, 1, y, gp = grid::gpar(lwd = lwd))
  }
  g <- gtable::gtable_add_grob(g,
    grobs = list(rule(1, 1.5), rule(0, 0.8), rule(0, 1.5)),
    t = c(1, 1, nrow(g)), l = 1, r = ncol(g), z = Inf, clip = "off",
    name = c("rule_top", "rule_mid", "rule_bottom")
  )
  ## Text is top-aligned in its row, so the space below it clears the
  ## descenders (grobHeight of text leaves them out).
  head <- list(title = grid::textGrob(title,
    x = 0, y = 1, hjust = 0, vjust = 1,
    gp = grid::gpar(fontsize = fontsize + 1, fontface = 2)
  ))
  if (!is.null(caption)) {
    head$caption <- grid::textGrob(caption,
      x = 0, y = 1, hjust = 0, vjust = 1,
      gp = grid::gpar(fontsize = fontsize, col = "grey30")
    )
  }
  ## Add from the bottom up, so the title ends up in the first row.
  for (name in rev(names(head))) {
    g <- gtable::gtable_add_rows(g,
      grid::grobHeight(head[[name]]) + grid::unit(2.5, "mm"),
      pos = 0
    )
    g <- gtable::gtable_add_grob(g, head[[name]],
      t = 1, l = 1, r = ncol(g), clip = "off", name = name
    )
  }
  g
}

## Write the tables of all results to one pdf, per_page tables stacked on each
## page. The page size is fitted to the widest table (or caption).
## res      a list of saved results, as returned by read_sim_results()
## file     path ending in .pdf
## main     title at the top of every page; NULL for none
## captions named vector of one-line descriptions, matched to names(res);
##          results without one get no caption
## note     line under main; NULL for none
## Returns the file path invisibly.
write_sim_tables <- function(
  res, file, main = NULL, captions = NULL, per_page = 3, fontsize = 9,
  note = "Entries are mean (standard error) across replicates."
) {
  stopifnot(
    length(res) > 0, !is.null(names(res)),
    grepl("\\.pdf$", file, ignore.case = TRUE), per_page >= 1
  )
  for (pkg in c("gridExtra", "gtable")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("package ", pkg, " is needed: install.packages(\"", pkg, "\")")
    }
  }
  no_table <- !vapply(res, function(r) is.matrix(r$table), logical(1))
  if (any(no_table)) {
    stop("no $table in: ", paste(names(res)[no_table], collapse = ", "))
  }
  grobs <- Map(function(r, name) {
    cap <- if (name %in% names(captions)) captions[[name]]
    sim_table_grob(r$table, name, cap, fontsize)
  }, res, names(res))

  ## Measure the tables on a null pdf device, which has the same font
  ## metrics as the pdf written below.
  cur <- grDevices::dev.cur()
  grDevices::pdf(NULL)
  inch_w <- function(u) grid::convertWidth(u, "in", valueOnly = TRUE)
  inch_h <- function(u) grid::convertHeight(u, "in", valueOnly = TRUE)
  ## tw is the table itself; w also covers a title or caption wider than it.
  tw <- vapply(grobs, function(g) inch_w(sum(g$widths)), 0)
  w <- vapply(seq_along(grobs), function(i) {
    g <- grobs[[i]]
    head <- g$grobs[g$layout$name %in% c("title", "caption")]
    max(tw[i], vapply(head, function(h) inch_w(grid::grobWidth(h)), 0))
  }, 0)
  h <- vapply(grobs, function(g) inch_h(sum(g$heights)), 0)
  grDevices::dev.off()
  if (cur > 1) grDevices::dev.set(cur)

  margin <- 0.5
  gap <- 0.35
  main_h <- if (is.null(main)) 0 else 0.3
  note_h <- if (is.null(note)) 0 else 0.25
  head_h <- main_h + note_h + if (main_h + note_h > 0) 0.15 else 0
  per_page <- min(per_page, length(grobs))
  width <- max(w) + 2 * margin
  height <- 2 * margin + head_h + per_page * max(h) + (per_page - 1) * gap

  grDevices::pdf(file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  pages <- split(seq_along(grobs), ceiling(seq_along(grobs) / per_page))
  for (idx in pages) {
    grid::grid.newpage()
    y <- height - margin
    if (!is.null(main)) {
      grid::grid.text(main,
        x = grid::unit(margin, "in"), y = grid::unit(y, "in"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = fontsize + 4, fontface = 2)
      )
      y <- y - main_h
    }
    if (!is.null(note)) {
      grid::grid.text(note,
        x = grid::unit(margin, "in"), y = grid::unit(y, "in"),
        just = c("left", "top"),
        gp = grid::gpar(fontsize = fontsize, col = "grey30")
      )
      y <- y - note_h
    }
    if (head_h > 0) y <- y - 0.15
    for (i in idx) {
      ## Each table gets a viewport of exactly its own size, left-aligned.
      grid::pushViewport(grid::viewport(
        x = grid::unit(margin, "in"), y = grid::unit(y, "in"),
        width = grid::unit(tw[i], "in"),
        height = grid::unit(h[i], "in"), just = c("left", "top")
      ))
      grid::grid.draw(grobs[[i]])
      grid::popViewport()
      y <- y - h[i] - gap
    }
  }
  invisible(file)
}
