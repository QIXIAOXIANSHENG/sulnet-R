## Summary tables (pdf) of the tuned rerun of models 1-6.
## Usage, from the package root:
##   Rscript inst/simulation/binomial_tuned/make_tables.R [res_dir] [pdf]
## Defaults: inst/simulation/binomial_tuned and
##   <res_dir>/binomial_tuned_tables.pdf
## Needs the gridExtra and gtable packages (see table_simulation.R).
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
sim_dir <- "inst/simulation"
res_dir <- arg(1, file.path(sim_dir, "binomial_tuned"))
pdf_file <- arg(2, file.path(res_dir, "binomial_tuned_tables.pdf"))
source(file.path(sim_dir, "table_simulation.R"))

res <- read_sim_results(res_dir)
n_rep <- vapply(res, function(r) if (is.null(r$n)) NA_integer_ else as.integer(r$n), integer(1))
write_sim_tables(res,
  file = pdf_file,
  main = sprintf(paste("Binomial simulations, tuned settings (%s replicates):",
    "suni2D standardize = FALSE, tau grid with 0, uniLasso loo = FALSE"),
    paste(unique(n_rep), collapse = "/")),
  captions = c(
    "Model 1" = "uniLasso, AR(1) rho = 0.8, class mean shift in x1-x20; n = 100, p = 300, 30 nonzero",
    "Model 2" = "independent; n = 200, p = 1000, 10 nonzero (+/-1)",
    "Model 3" = "AR(1) rho = 0.5; n = 200, p = 1000, 10 nonzero (+/-1)",
    "Model 4" = "block, size 20, rho = 0.5; n = 200, p = 1000, 10 nonzero (+/-1)",
    "Model 5" = "imbalanced, AR(1) rho = 0.5, prevalence 0.2 (intercept -2.19); n = 200, p = 1000, 10 nonzero (+/-1)",
    "Model 6" = "Breheny-Huang (2011), independent; n = 400, p = 1000, 4 nonzero (+/-1)"
  )[names(res)]
)
message("written ", pdf_file)
