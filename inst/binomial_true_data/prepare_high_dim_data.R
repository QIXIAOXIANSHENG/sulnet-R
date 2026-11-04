## Organize the binary datasets downloaded from
## https://iowabiostat.github.io/data-sets/ into analysis-ready CSVs.
## Run from the package root:
##   source("inst/binomial_true_data/prepare_high_dim_data.R")
##
## Input: the files in inst/binomial_true_data/data/high_dim_data/.
## Output: one CSV per dataset in inst/binomial_true_data/data/csv/, in the
## same layout as prepare_uci_data.R: one column per predictor and the 0/1
## response "target" as the last column. Which class is coded 1 is stated in
## each section below.
##
## No file here contains a separate test set (X.test / y.test). If one did,
## read_xy() would append it to the training rows.

hd_dir <- "./inst/binomial_true_data/data/high_dim_data"
out_dir <- "./inst/binomial_true_data/data/csv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## Read a list holding a predictor matrix (X or x) and a response y. A supplied
## test set (X.test, y.test) is appended to the training rows. [[ ]] is used
## because $ would partially match X to X.test.
read_xy <- function(path) {
  dat <- readRDS(path)
  x <- if (!is.null(dat[["X"]])) dat[["X"]] else dat[["x"]]
  y <- dat[["y"]]
  if (!is.null(dat[["X.test"]])) {
    stopifnot(identical(colnames(dat[["X.test"]]), colnames(x)))
    x <- rbind(x, dat[["X.test"]])
    y <- c(y, dat[["y.test"]])
  }
  stopifnot(nrow(x) == length(y))
  list(x = x, y = y)
}

## Write one dataset to out_dir/<name>.csv: the predictors x (a matrix or data
## frame) and the 0/1 response y as the last column "target". Columns constant
## over all rows are dropped, and predictor names are made syntactic and
## unique (as read.csv would; Golub1999 and Singh2002 repeat gene symbols).
## Returns a one-row summary.
write_binary_csv <- function(x, y, name) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  stopifnot(
    nrow(x) == length(y), all(y %in% c(0, 1)), length(unique(y)) == 2,
    !anyNA(x)
  )
  const <- vapply(x, function(v) all(v == v[1]), logical(1))
  if (any(const)) {
    message(name, ": dropping ", sum(const), " constant column(s)")
    x <- x[!const]
  }
  names(x) <- make.names(names(x), unique = TRUE)
  stopifnot(!"target" %in% names(x))
  x$target <- as.integer(y)
  write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  data.frame(
    dataset = name, n = nrow(x), p = ncol(x) - 1L,
    text_cols = sum(vapply(x, is.character, logical(1))),
    n_target1 = sum(y), n_target0 = sum(y == 0)
  )
}

info <- list()

## ---------------------------------------------------------------------------
## Gode2011: surface sensing in swarming bacteria
## https://iowabiostat.github.io/data-sets/Gode2011/Gode2011.html
## Y is 5138 genes x 8 arrays (log2 expression). expCond gives each array's
## Celltype (swimmer / swarmer) and Environment (plate / liquid): a 2 x 2
## design with 2 arrays per combination. Here the arrays are the observations
## and the genes the predictors (8 x 5138), and Celltype is the response.
## target = 1: swarmer (0: swimmer).
## Environment is the other experimental factor, not a measurement, and is
## not used as a predictor. geneInfo (gene annotations) is not used.
## With 8 observations, a 50-50 split leaves 2 arrays per class for training,
## too few for cross-validated logistic fits.
## ---------------------------------------------------------------------------
gode <- readRDS(file.path(hd_dir, "Gode2011.rds"))
x <- t(gode$Y)
cond <- gode$expCond
## The array names (TP1, TP2, OP1, ...) encode the condition: the 2nd letter
## is the environment (P / L) and the 1st letter the cell type (T = swarmer,
## O = swimmer). Check that the rows of expCond follow the columns of Y.
stopifnot(
  nrow(cond) == nrow(x),
  substr(rownames(x), 2, 2) == toupper(substr(cond$Environment, 1, 1)),
  substr(rownames(x), 1, 1) == ifelse(cond$Celltype == "swarmer", "T", "O")
)
info$gode <- write_binary_csv(
  x, as.numeric(cond$Celltype == "swarmer"), "Gode2011"
)

## ---------------------------------------------------------------------------
## Golub1999: gene expression in leukemia patients
## https://iowabiostat.github.io/data-sets/Golub1999/Golub1999.html
## X: 72 bone marrow samples x 7129 probes (log2, vsn-normalized), named by
## gene symbol, or by the Affymetrix probe when no symbol is known.
## y: ALL (47) / AML (25). target = 1: AML.
## ---------------------------------------------------------------------------
dat <- read_xy(file.path(hd_dir, "Golub1999.rds"))
info$golub <- write_binary_csv(
  dat$x, as.numeric(dat$y == "AML"), "Golub1999"
)

## ---------------------------------------------------------------------------
## Scherzer2007: blood biomarkers for early Parkinson's disease
## https://iowabiostat.github.io/data-sets/Scherzer2007/Scherzer2007.html
## X: 105 blood samples x 22283 probes, named by gene symbol.
## y: logical. The site lists 50 Parkinson's patients (cases) and 55 controls,
## and y has 50 TRUE and 55 FALSE, so TRUE = case.
## target = 1: Parkinson's disease (y = TRUE).
## ---------------------------------------------------------------------------
dat <- read_xy(file.path(hd_dir, "Scherzer2007.rds"))
stopifnot(sum(dat$y) == 50, sum(!dat$y) == 55)
info$scherzer <- write_binary_csv(
  dat$x, as.numeric(dat$y), "Scherzer2007"
)

## ---------------------------------------------------------------------------
## Singh2002: case-control study of prostate cancer
## https://iowabiostat.github.io/data-sets/Singh2002/Singh2002.html
## X: 102 prostatectomy samples x 12600 probes (log2, vsn-normalized), named
## by GenBank ID. y: normal (50) / tumor (52). target = 1: tumor.
## ---------------------------------------------------------------------------
dat <- read_xy(file.path(hd_dir, "Singh2002.rds"))
info$singh <- write_binary_csv(
  dat$x, as.numeric(dat$y == "tumor"), "Singh2002"
)

## ---------------------------------------------------------------------------
## glcamd: genetic association study of glaucoma and macular degeneration
## https://iowabiostat.github.io/data-sets/glcamd/glcamd.html
## x: 800 patients x 699 SNPs. y: AMD (400) / GLC (400); each disease group
## is the control group for the other. target = 1: GLC (glaucoma).
## group (the nearest gene of each SNP) is not used.
## ---------------------------------------------------------------------------
dat <- read_xy(file.path(hd_dir, "glcamd.rds"))
x <- dat$x
## The site describes the genotypes as 0/1/2, but 8.7% of the entries are -1,
## in every patient and in 681 of the 699 SNPs. -1 is read as a missing
## genotype. Following the usual GWAS quality control, SNPs with a call rate
## below 90% (more than 80 of 800 genotypes missing) are dropped. Each
## remaining missing genotype is replaced by its SNP's mean over all 800
## patients; the response is not used.
x[x == -1] <- NA
x <- x[, colMeans(is.na(x)) <= 0.1]
miss <- which(is.na(x), arr.ind = TRUE)
x[miss] <- colMeans(x, na.rm = TRUE)[miss[, 2]]
info$glcamd <- write_binary_csv(
  x, as.numeric(dat$y == "GLC"), "glcamd"
)

## ---------------------------------------------------------------------------
## support: SUPPORT study of seriously ill hospitalized patients
## https://iowabiostat.github.io/data-sets/support/support.html
## 1227 patients. Response: Death (death indicator). target = 1: died.
## Income (4 brackets) and Race (5 groups) are text; the rest are numeric. No
## missing values.
## ---------------------------------------------------------------------------
sup <- read.delim(file.path(hd_dir, "support.txt"),
  check.names = FALSE, stringsAsFactors = FALSE
)
## Removed: Time (days on study). With Death it forms the survival outcome, so
## it leaks the response.
info$support <- write_binary_csv(
  sup[!names(sup) %in% c("Death", "Time")], sup$Death, "support"
)

## ---------------------------------------------------------------------------
## Summary of the written files
## ---------------------------------------------------------------------------
high_dim_summary <- do.call(rbind, info)
rownames(high_dim_summary) <- NULL
print(high_dim_summary)
