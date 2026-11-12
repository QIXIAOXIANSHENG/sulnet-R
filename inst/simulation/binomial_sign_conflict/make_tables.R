## Summary tables (pdf) of the sign-conflict simulations, models 7-11.
## Usage, from the package root:
##   Rscript inst/simulation/binomial_sign_conflict/make_tables.R [res_dir] [pdf]
## Defaults: inst/simulation/binomial_sign_conflict and
##   <res_dir>/binomial_sign_conflict_tables.pdf
## Needs the gridExtra and gtable packages (see table_simulation.R).
args <- commandArgs(trailingOnly = TRUE)
arg <- function(i, default) if (length(args) >= i && nzchar(args[i])) args[i] else default
sim_dir <- "inst/simulation"
res_dir <- arg(1, file.path(sim_dir, "binomial_sign_conflict"))
pdf_file <- arg(2, file.path(res_dir, "binomial_sign_conflict_tables.pdf"))
source(file.path(sim_dir, "table_simulation.R"))

res <- read_sim_results(res_dir)
n_rep <- vapply(res, function(r) if (is.null(r$n)) NA_integer_ else as.integer(r$n), integer(1))
write_sim_tables(res,
  file = pdf_file,
  main = sprintf("Binomial sign-conflict simulations (%s replicates)",
    paste(unique(n_rep), collapse = "/")),
  captions = c(
    "Model 7" = "blocks of 10, rho = 0.5, beta (1,1,1,-1) in blocks 1-3; n = 400, p = 1000; 3 sign conflicts",
    "Model 8" = "blocks of 10, rho = 0.7, beta (1,1,1,-1.5) in blocks 1-3; n = 400, p = 1000; 3 sign conflicts",
    "Model 9" = "blocks of 10, rho = 0.8, beta (2,2,-2) in blocks 1-3; n = 400, p = 1000; 3 sign conflicts",
    "Model 10" = "AR(1) rho = 0.8, beta (1,1,1,-1.5) at x1-4, x101-104, x201-204; n = 400, p = 1000; 3 sign conflicts",
    "Model 11" = "control without conflicts: blocks of 20, rho = 0.9, beta (1,1) in blocks 1-5; n = 100, p = 1000"
  )[names(res)]
)
message("written ", pdf_file)
