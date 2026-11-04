## Organize the UCI binary-classification datasets into analysis-ready CSVs.
## Run from the package root:
##   source("inst/binomial_true_data/prepare_uci_data.R")
##
## Input: the zips in inst/binomial_true_data/data/UCI_data/. The first 12
## were downloaded by hand from the UCI repository. The additional datasets in
## the second part are downloaded into the same folder by download_uci() if
## their zip is missing. Every zip is extracted to its own temporary folder;
## nothing is extracted into the repository.
##
## Output: one CSV per dataset in inst/binomial_true_data/data/csv/, with
##   - one column per predictor. Numeric codes stay numeric, except codes of
##     unordered categories (blood group, chest-pain type, ...), which are
##     replaced by text labels. Text columns are read back by read.csv as
##     character, and true_dataset.R turns them into dummy variables.
##   - the 0/1 response "target" as the last column. Which class is coded 1 is
##     stated in each section below.
## Unless a section says otherwise, rows with a missing value are dropped
## (complete cases). Identifier columns, columns that leak the response, and
## columns constant over all rows are removed.
##
## Two of the downloaded datasets are not binary: Wine (3 classes) and Image
## Segmentation (7 classes). They are made binary as one class against the
## rest, as stated in their sections.

uci_dir <- "./inst/binomial_true_data/data/UCI_data"
out_dir <- "./inst/binomial_true_data/data/csv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(600, getOption("timeout")))

## Extract one zip from uci_dir into a fresh temporary folder and return the
## folder.
unzip_uci <- function(zip_name) {
  exdir <- file.path(tempdir(), sub("\\.zip$", "", zip_name))
  unlink(exdir, recursive = TRUE)
  dir.create(exdir, recursive = TRUE)
  utils::unzip(file.path(uci_dir, zip_name), exdir = exdir)
  exdir
}

## Download a dataset's zip from the UCI file server into uci_dir, unless it is
## already there. Returns the zip's file name.
download_uci <- function(url) {
  zip_name <- basename(url)
  dest <- file.path(uci_dir, zip_name)
  if (!file.exists(dest)) download.file(url, dest, mode = "wb")
  zip_name
}

## Minimal ARFF reader. foreign::read.arff() fails on the Chronic Kidney
## Disease file, whose data rows contain stray tabs and spaces, two trailing
## commas and one doubled comma. Here every field is trimmed and empty fields
## are dropped, after which each row must have one field per @attribute. "?" is
## read as NA. All columns are returned as character, for the caller to convert.
read_arff_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  attr_lines <- grep("^\\s*@attribute", lines, ignore.case = TRUE, value = TRUE)
  vars <- gsub("'", "", sub(
    "^\\s*@attribute\\s+(\\S+).*$", "\\1", attr_lines,
    ignore.case = TRUE
  ))
  start <- grep("^\\s*@data", lines, ignore.case = TRUE)
  body <- lines[-seq_len(start)]
  body <- body[nzchar(trimws(body)) & !grepl("^\\s*%", body)]
  fields <- lapply(strsplit(body, ","), function(f) {
    f <- trimws(f)
    f[nzchar(f)]
  })
  stopifnot(all(lengths(fields) == length(vars)))
  m <- do.call(rbind, fields)
  m[m == "?"] <- NA
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- vars
  df
}

## Write one dataset to out_dir/<name>.csv: the predictors in the data frame x
## and the 0/1 response y as the last column "target". Columns constant over
## all rows are dropped, and predictor names are made syntactic and unique (as
## read.csv would). Returns a one-row summary.
write_binary_csv <- function(x, y, name) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  stopifnot(
    nrow(x) == length(y), all(y %in% c(0, 1)), length(unique(y)) == 2,
    !anyNA(x)
  )
  const <- vapply(x, function(v) all(v == v[1]), logical(1))
  if (any(const)) {
    message(name, ": dropping constant column(s) ",
            paste(names(x)[const], collapse = ", "))
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

## ===========================================================================
## Part 1: the datasets downloaded by hand
## ===========================================================================

## ---------------------------------------------------------------------------
## Bone marrow transplant: children
## https://archive.ics.uci.edu/dataset/565/bone+marrow+transplant+children
## 187 children after unrelated-donor stem cell transplantation.
## target = 1: died (survival_status = 1).
## ---------------------------------------------------------------------------
dir <- unzip_uci("bone+marrow+transplant+children.zip")
bmt <- read_arff_simple(file.path(dir, "bone-marrow.arff"))
y <- as.numeric(bmt$survival_status)
## Removed:
## - survival_time: the follow-up time, which encodes the outcome.
## - ANCrecovery, PLTrecovery, time_to_aGvHD_III_IV: times to post-transplant
##   events, coded 1e6 when the event never happened, so the code is not a
##   time. ANCrecovery = 1e6 occurs only for 5 patients who died, so it
##   leaks the outcome. The 1e6 code in time_to_aGvHD_III_IV repeats the
##   binary aGvHDIIIIV, which is kept.
## - extcGvHD (extensive chronic GvHD): missing for 31 patients, all of whom
##   died within 96 days, before it could be assessed. Its missingness
##   encodes the outcome, and dropping those rows instead would remove 31 of
##   the 85 deaths.
drop <- c(
  "survival_status", "survival_time", "ANCrecovery", "PLTrecovery",
  "time_to_aGvHD_III_IV", "extcGvHD"
)
x <- bmt[setdiff(names(bmt), drop)]
## Unordered categories become text: blood groups (codes from the ARFF
## header: 0 = O, 1 = A, -1 = B, 2 = AB), Disease (already text) and HLAgrI
## (type of HLA difference). Every other column is binary or an ordered code
## and stays numeric.
abo <- c("0" = "O", "1" = "A", "-1" = "B", "2" = "AB")
x$DonorABO <- unname(abo[x$DonorABO])
x$RecipientABO <- unname(abo[x$RecipientABO])
x$HLAgrI <- ifelse(is.na(x$HLAgrI), NA, paste0("type", x$HLAgrI))
num <- setdiff(names(x), c("DonorABO", "RecipientABO", "Disease", "HLAgrI"))
x[num] <- lapply(x[num], as.numeric)
## 23 rows have another missing value (CMVstatus 16, RecipientCMV 14, ...):
## 164 remain (70 died, 94 survived).
ok <- complete.cases(x)
info$bmt <- write_binary_csv(x[ok, ], y[ok], "bone_marrow_transplant")

## ---------------------------------------------------------------------------
## Breast Cancer Wisconsin (Diagnostic)
## https://archive.ics.uci.edu/dataset/17/breast+cancer+wisconsin+diagnostic
## 569 breast masses, 30 features of the cell nuclei (mean, standard error
## and worst value of 10 measurements). No missing values.
## target = 1: malignant (M).
## ---------------------------------------------------------------------------
dir <- unzip_uci("breast+cancer+wisconsin+diagnostic.zip")
wdbc <- read.csv(file.path(dir, "wdbc.data"), header = FALSE)
feat <- c(
  "radius", "texture", "perimeter", "area", "smoothness", "compactness",
  "concavity", "concave_points", "symmetry", "fractal_dimension"
)
names(wdbc) <- c(
  "id", "diagnosis",
  paste0(feat, "_", rep(c("mean", "se", "worst"), each = 10))
)
## Removed: id (patient identifier).
info$wdbc <- write_binary_csv(
  wdbc[-(1:2)], as.numeric(wdbc$diagnosis == "M"),
  "breast_cancer_wisconsin_diagnostic"
)

## ---------------------------------------------------------------------------
## Breast Cancer (Ljubljana)
## https://archive.ics.uci.edu/dataset/14/breast+cancer
## 286 patients. The predictors are text categories (age, tumor-size and
## inv-nodes are interval labels such as "30-39") except deg-malig (1-3).
## target = 1: recurrence-events.
## ---------------------------------------------------------------------------
dir <- unzip_uci("breast+cancer.zip")
bc <- read.csv(file.path(dir, "breast-cancer.data"),
  header = FALSE, na.strings = "?", stringsAsFactors = FALSE
)
names(bc) <- c(
  "class", "age", "menopause", "tumor_size", "inv_nodes", "node_caps",
  "deg_malig", "breast", "breast_quad", "irradiat"
)
## 9 rows have a missing value (node_caps 8, breast_quad 1): 277 remain.
bc <- bc[complete.cases(bc), ]
info$bc <- write_binary_csv(
  bc[-1], as.numeric(bc$class == "recurrence-events"), "breast_cancer"
)

## ---------------------------------------------------------------------------
## Chronic Kidney Disease
## https://archive.ics.uci.edu/dataset/336/chronic+kidney+disease
## 400 patients (250 ckd, 150 notckd), 24 predictors. target = 1: ckd.
## ---------------------------------------------------------------------------
dir <- unzip_uci("chronic+kidney+disease.zip")
## The zip contains a RAR archive, which R cannot read. bsdtar (libarchive)
## can: Windows 10+ ships it as System32/tar.exe (the tar on the PATH may be
## GNU tar, which cannot), and macOS ships it as bsdtar.
tar_exe <- if (.Platform$OS.type == "windows") {
  file.path(Sys.getenv("SystemRoot"), "System32", "tar.exe")
} else {
  "bsdtar"
}
status <- system2(tar_exe, c(
  "-xf", shQuote(file.path(dir, "Chronic_Kidney_Disease.rar")),
  "-C", shQuote(dir)
))
stopifnot(status == 0)
## chronic_kidney_disease.arff and chronic_kidney_disease_full.arff hold the
## same 400 rows; the "full" file is used.
ckd <- read_arff_simple(file.path(
  dir, "Chronic_Kidney_Disease", "chronic_kidney_disease_full.arff"
))
y <- as.numeric(ckd$class == "ckd")
x <- ckd[names(ckd) != "class"]
## Numeric: the 11 numeric attributes, plus sg (specific gravity), al
## (albumin) and su (sugar), which the ARFF header declares nominal but which
## are ordered numeric grades. The other 10 attributes are yes/no-type text.
num <- c(
  "age", "bp", "sg", "al", "su", "bgr", "bu", "sc", "sod", "pot", "hemo",
  "pcv", "wbcc", "rbcc"
)
x[num] <- lapply(x[num], as.numeric)
## Missing values are common (rbc 152, rbcc 131, wbcc 106, ...). Complete
## cases: 158 rows (43 ckd, 115 notckd).
ok <- complete.cases(x)
info$ckd <- write_binary_csv(x[ok, ], y[ok], "chronic_kidney_disease")

## ---------------------------------------------------------------------------
## Credit Approval
## https://archive.ics.uci.edu/dataset/27/credit+approval
## 690 credit card applications, 15 anonymized predictors A1-A15
## (A2, A3, A8, A11, A14, A15 numeric; the rest text). target = 1: "+".
## ---------------------------------------------------------------------------
dir <- unzip_uci("credit+approval.zip")
crx <- read.csv(file.path(dir, "crx.data"),
  header = FALSE, na.strings = "?", stringsAsFactors = FALSE,
  colClasses = c(
    "character", "numeric", "numeric", rep("character", 4), "numeric",
    "character", "character", "numeric", "character", "character",
    "numeric", "numeric", "character"
  )
)
names(crx) <- paste0("A", 1:16)
## 37 rows have a missing value: 653 remain.
crx <- crx[complete.cases(crx), ]
info$crx <- write_binary_csv(
  crx[-16], as.numeric(crx$A16 == "+"), "credit_approval"
)

## ---------------------------------------------------------------------------
## Drug Induced Autoimmunity Prediction
## https://archive.ics.uci.edu/dataset/1104/drug_induced_autoimmunity_prediction
## 196 RDKit molecular descriptors. The zip has a training set (477 drugs) and
## a test set (120 drugs); they are combined (597 rows). No missing values.
## target = 1: Label = 1 (autoimmunity-inducing drug).
## ---------------------------------------------------------------------------
dir <- unzip_uci("drug_induced_autoimmunity_prediction.zip")
dia <- rbind(
  read.csv(file.path(dir, "DIA_trainingset_RDKit_descriptors.csv"),
    check.names = FALSE
  ),
  read.csv(file.path(dir, "DIA_testset_RDKit_descriptors.csv"),
    check.names = FALSE
  )
)
## Removed: SMILES (the molecule's identifier). RDKit_ChemDes.xlsx only
## describes the descriptors and is not used.
info$dia <- write_binary_csv(
  dia[!names(dia) %in% c("Label", "SMILES")], dia$Label,
  "drug_induced_autoimmunity"
)

## ---------------------------------------------------------------------------
## Heart Disease (Cleveland)
## https://archive.ics.uci.edu/dataset/45/heart+disease
## processed.cleveland.data: 303 patients, the 14 attributes that published
## experiments use. The other three sites (Hungary, Switzerland, VA) have
## many missing values and are not used. num is 0 (no disease) to 4; as in
## heart-disease.names, presence (1-4) is distinguished from absence (0).
## target = 1: num > 0.
## ---------------------------------------------------------------------------
dir <- unzip_uci("heart+disease.zip")
heart <- read.csv(file.path(dir, "processed.cleveland.data"),
  header = FALSE, na.strings = "?"
)
names(heart) <- c(
  "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg", "thalach",
  "exang", "oldpeak", "slope", "ca", "thal", "num"
)
## 6 rows have a missing value (ca 4, thal 2): 297 remain.
heart <- heart[complete.cases(heart), ]
## Unordered categories become text, with labels from heart-disease.names.
heart$cp <- c(
  "typical_angina", "atypical_angina", "non_anginal", "asymptomatic"
)[heart$cp]
heart$restecg <- c("normal", "st_t_abnormal", "lv_hypertrophy")[heart$restecg + 1]
heart$slope <- c("up", "flat", "down")[heart$slope]
heart$thal <- c("3" = "normal", "6" = "fixed_defect", "7" = "reversible_defect")[
  as.character(heart$thal)
]
info$heart <- write_binary_csv(
  heart[names(heart) != "num"], as.numeric(heart$num > 0),
  "heart_disease_cleveland"
)

## ---------------------------------------------------------------------------
## Heart Failure Clinical Records
## https://archive.ics.uci.edu/dataset/519/heart+failure+clinical+records
## 299 patients, 12 predictors. No missing values.
## target = 1: DEATH_EVENT = 1.
## ---------------------------------------------------------------------------
dir <- unzip_uci("heart+failure+clinical+records.zip")
hf <- read.csv(file.path(dir, "heart_failure_clinical_records_dataset.csv"))
## Removed: time, the follow-up period. Patients who died have short
## follow-up, so it leaks the outcome.
info$hf <- write_binary_csv(
  hf[!names(hf) %in% c("time", "DEATH_EVENT")], hf$DEATH_EVENT,
  "heart_failure"
)

## ---------------------------------------------------------------------------
## Image Segmentation (NOT binary: 7 classes of 330 regions each)
## https://archive.ics.uci.edu/dataset/50/image+segmentation
## segmentation.data (210 rows) and segmentation.test (2100 rows) are combined
## (2310 rows); 19 numeric features. Made binary as one class against the
## other 6. The classes are equally frequent, so the majority-class rule used
## for Wine below does not pick one; BRICKFACE, the first class in the file,
## is used. target = 1: BRICKFACE.
## ---------------------------------------------------------------------------
dir <- unzip_uci("image+segmentation.zip")
## Each file starts with 2 comment lines, a blank line, the 19 feature names
## and another blank line; the class is the first field of each data row.
read_segmentation <- function(path) {
  lines <- readLines(path)
  d <- read.csv(text = lines[-(1:5)], header = FALSE, stringsAsFactors = FALSE)
  names(d) <- c("class", strsplit(lines[4], ",")[[1]])
  d
}
seg <- rbind(
  read_segmentation(file.path(dir, "segmentation.data")),
  read_segmentation(file.path(dir, "segmentation.test"))
)
## REGION-PIXEL-COUNT is 9 for every region and is dropped as constant.
info$seg <- write_binary_csv(
  seg[-1], as.numeric(seg$class == "BRICKFACE"),
  "image_segmentation_brickface"
)

## ---------------------------------------------------------------------------
## Mushroom
## https://archive.ics.uci.edu/dataset/73/mushroom
## 8124 mushrooms, 22 text attributes (one-letter codes listed in
## agaricus-lepiota.names). target = 1: poisonous (p).
## ---------------------------------------------------------------------------
dir <- unzip_uci("mushroom.zip")
mush <- read.csv(file.path(dir, "agaricus-lepiota.data"),
  header = FALSE, colClasses = "character"
)
names(mush) <- c(
  "class", "cap_shape", "cap_surface", "cap_color", "bruises", "odor",
  "gill_attachment", "gill_spacing", "gill_size", "gill_color",
  "stalk_shape", "stalk_root", "stalk_surface_above_ring",
  "stalk_surface_below_ring", "stalk_color_above_ring",
  "stalk_color_below_ring", "veil_type", "veil_color", "ring_number",
  "ring_type", "spore_print_color", "population", "habitat"
)
## stalk_root is "?" in 2480 rows (30%), and no other attribute is missing.
## Dropping those rows would change the population, so "?" is kept as its
## own level, "missing". veil_type is "p" in every row and is dropped as
## constant.
mush$stalk_root[mush$stalk_root == "?"] <- "missing"
info$mush <- write_binary_csv(
  mush[-1], as.numeric(mush$class == "p"), "mushroom"
)

## ---------------------------------------------------------------------------
## Statlog (German Credit Data)
## https://archive.ics.uci.edu/dataset/144/statlog+german+credit+data
## german.data: 1000 applicants, 20 attributes (7 numeric, 13 qualitative
## with codes such as A11, described in german.doc). german.data-numeric, the
## Strathclyde recoding, is not used. No missing values.
## target = 1: bad credit risk (class 2).
## ---------------------------------------------------------------------------
dir <- unzip_uci("statlog+german+credit+data.zip")
german <- read.table(file.path(dir, "german.data"), stringsAsFactors = FALSE)
names(german) <- c(
  "checking_status", "duration", "credit_history", "purpose",
  "credit_amount", "savings", "employment_since", "installment_rate",
  "personal_status_sex", "other_debtors", "residence_since", "property",
  "age", "other_installment_plans", "housing", "n_existing_credits", "job",
  "n_people_liable", "telephone", "foreign_worker", "class"
)
info$german <- write_binary_csv(
  german[-21], as.numeric(german$class == 2), "german_credit"
)

## ---------------------------------------------------------------------------
## Wine (NOT binary: classes 1, 2, 3 with 59, 71, 48 wines)
## https://archive.ics.uci.edu/dataset/109/wine
## 178 wines, 13 chemical measurements. No missing values. Made binary as the
## most frequent class (2) against the rest, the rule OpenML uses for its
## binarized UCI datasets. target = 1: class 2.
## ---------------------------------------------------------------------------
dir <- unzip_uci("wine.zip")
wine <- read.csv(file.path(dir, "wine.data"), header = FALSE)
names(wine) <- c(
  "class", "alcohol", "malic_acid", "ash", "alcalinity_of_ash", "magnesium",
  "total_phenols", "flavanoids", "nonflavanoid_phenols", "proanthocyanins",
  "color_intensity", "hue", "od280_od315", "proline"
)
info$wine <- write_binary_csv(
  wine[-1], as.numeric(wine$class == 2), "wine_class2"
)

## ===========================================================================
## Part 2: commonly used binary UCI datasets that were not downloaded,
## fetched from the UCI file server
## ===========================================================================

## ---------------------------------------------------------------------------
## Ionosphere
## https://archive.ics.uci.edu/dataset/52/ionosphere
## 351 radar returns, 34 numeric features. No missing values.
## target = 1: good return (g).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/52/ionosphere.zip"
))
iono <- read.csv(file.path(dir, "ionosphere.data"), header = FALSE)
names(iono) <- c(paste0("A", 1:34), "class")
## A2 is 0 in every row and is dropped as constant.
info$iono <- write_binary_csv(
  iono[-35], as.numeric(iono$class == "g"), "ionosphere"
)

## ---------------------------------------------------------------------------
## Connectionist Bench (Sonar, Mines vs. Rocks)
## https://archive.ics.uci.edu/dataset/151/connectionist+bench+sonar+mines+vs+rocks
## 208 sonar returns, energy in 60 frequency bands. No missing values.
## target = 1: mine (M).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/151/connectionist+bench+sonar+mines+vs+rocks.zip"
))
sonar <- read.csv(file.path(dir, "sonar.all-data"), header = FALSE)
names(sonar) <- c(paste0("band", 1:60), "class")
info$sonar <- write_binary_csv(
  sonar[-61], as.numeric(sonar$class == "M"), "sonar"
)

## ---------------------------------------------------------------------------
## Spambase
## https://archive.ics.uci.edu/dataset/94/spambase
## 4601 e-mails, 57 numeric features. No missing values.
## target = 1: spam.
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/94/spambase.zip"
))
spam <- read.csv(file.path(dir, "spambase.data"), header = FALSE)
## Feature names from spambase.names; the six char_freq_<symbol> names are
## spelled out so that they stay distinct after make.names().
spam_names <- sub(":.*$", "", grep(
  "^(word_freq|char_freq|capital_run)",
  readLines(file.path(dir, "spambase.names")),
  value = TRUE
))
stopifnot(length(spam_names) == 57, ncol(spam) == 58)
spam_names[49:54] <- paste0(
  "char_freq_", c("semicolon", "paren", "bracket", "exclam", "dollar", "hash")
)
names(spam) <- c(spam_names, "spam")
info$spam <- write_binary_csv(spam[-58], spam$spam, "spambase")

## ---------------------------------------------------------------------------
## Parkinsons
## https://archive.ics.uci.edu/dataset/174/parkinsons
## 195 voice recordings, 22 numeric features. No missing values. The
## recordings come from 32 people (about 6 each), so a random split puts the
## same person in both halves; the dataset is commonly used this way.
## target = 1: Parkinson's disease (status = 1).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/174/parkinsons.zip"
))
park <- read.csv(file.path(dir, "parkinsons.data"), check.names = FALSE)
## Removed: name (recording identifier). The telemonitoring/ folder is a
## different (regression) dataset and is not used.
info$park <- write_binary_csv(
  park[!names(park) %in% c("name", "status")], park$status, "parkinsons"
)

## ---------------------------------------------------------------------------
## Breast Cancer Wisconsin (Original)
## https://archive.ics.uci.edu/dataset/15/breast+cancer+wisconsin+original
## 699 samples, 9 cytology scores (1-10). target = 1: malignant (class 4).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/15/breast+cancer+wisconsin+original.zip"
))
bcw <- read.csv(file.path(dir, "breast-cancer-wisconsin.data"),
  header = FALSE, na.strings = "?"
)
names(bcw) <- c(
  "id", "clump_thickness", "cell_size_uniformity",
  "cell_shape_uniformity", "marginal_adhesion", "single_epithelial_cell_size",
  "bare_nuclei", "bland_chromatin", "normal_nucleoli", "mitoses", "class"
)
## Removed: id (sample code). bare_nuclei is missing in 16 rows: 683 remain.
## The zip also contains wdbc.data (the Diagnostic dataset above) and
## wpbc.data (a prognostic dataset); neither is used here.
bcw <- bcw[complete.cases(bcw), ]
info$bcw <- write_binary_csv(
  bcw[2:10], as.numeric(bcw$class == 4), "breast_cancer_wisconsin_original"
)

## ---------------------------------------------------------------------------
## ILPD (Indian Liver Patient Dataset)
## https://archive.ics.uci.edu/dataset/225/ilpd+indian+liver+patient+dataset
## 583 patients, 10 predictors (Gender is text). The CSV has no header.
## target = 1: liver patient (Selector = 1).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/225/ilpd+indian+liver+patient+dataset.zip"
))
ilpd <- read.csv(file.path(dir, "Indian Liver Patient Dataset (ILPD).csv"),
  header = FALSE, stringsAsFactors = FALSE
)
names(ilpd) <- c(
  "age", "gender", "tb", "db", "alkphos", "sgpt", "sgot", "tp", "alb",
  "ag_ratio", "selector"
)
## ag_ratio is missing in 4 rows: 579 remain.
ilpd <- ilpd[complete.cases(ilpd), ]
info$ilpd <- write_binary_csv(
  ilpd[-11], as.numeric(ilpd$selector == 1), "ilpd"
)

## ---------------------------------------------------------------------------
## QSAR biodegradation
## https://archive.ics.uci.edu/dataset/254/qsar+biodegradation
## 1055 chemicals, 41 numeric molecular descriptors (named on the UCI page; the
## file has no header, so they are called V1-V41 here). No missing values.
## target = 1: ready biodegradable (RB).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/254/qsar+biodegradation.zip"
))
qsar <- read.csv(file.path(dir, "biodeg.csv"), sep = ";", header = FALSE)
info$qsar <- write_binary_csv(
  qsar[1:41], as.numeric(qsar[[42]] == "RB"), "qsar_biodegradation"
)

## ---------------------------------------------------------------------------
## Congressional Voting Records
## https://archive.ics.uci.edu/dataset/105/congressional+voting+records
## 435 members of the 1984 US House, 16 votes (y / n / ?).
## target = 1: republican.
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/105/congressional+voting+records.zip"
))
votes <- read.csv(file.path(dir, "house-votes-84.data"),
  header = FALSE, colClasses = "character"
)
names(votes) <- c(
  "party", "handicapped_infants", "water_project_cost_sharing",
  "adoption_of_the_budget_resolution", "physician_fee_freeze",
  "el_salvador_aid", "religious_groups_in_schools", "anti_satellite_test_ban",
  "aid_to_nicaraguan_contras", "mx_missile", "immigration",
  "synfuels_corporation_cutback", "education_spending",
  "superfund_right_to_sue", "crime", "duty_free_exports",
  "export_administration_act_south_africa"
)
## As house-votes-84.names notes, "?" does not mean unknown: the member did
## not vote yea or nay. It is kept as a third level, "other".
votes[-1][votes[-1] == "?"] <- "other"
info$votes <- write_binary_csv(
  votes[-1], as.numeric(votes$party == "republican"),
  "congressional_voting"
)

## ---------------------------------------------------------------------------
## SPECTF Heart
## https://archive.ics.uci.edu/dataset/96/spectf+heart
## 44 integer features (rest and stress counts, F1R, F1S, ..., F22S) from
## cardiac SPECT images. SPECTF.train (80 rows) and SPECTF.test (187 rows)
## are combined (267 rows). SPECTFincorrect.test is the superseded test file
## (see DonorNote.txt) and is not used. target = 1: abnormal (diagnosis = 1).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/96/spectf+heart.zip"
))
spectf <- rbind(
  read.csv(file.path(dir, "SPECTF.train"), header = FALSE),
  read.csv(file.path(dir, "SPECTF.test"), header = FALSE)
)
names(spectf) <- c("diagnosis", paste0("F", rep(1:22, each = 2), c("R", "S")))
info$spectf <- write_binary_csv(spectf[-1], spectf$diagnosis, "spectf_heart")

## ---------------------------------------------------------------------------
## Climate Model Simulation Crashes
## https://archive.ics.uci.edu/dataset/252/climate+model+simulation+crashes
## 540 simulations (3 Latin hypercube studies of 180 runs), 18 model
## parameters in [0, 1]. outcome is 0 for a failed run and 1 for a success.
## target = 1: failed run (outcome = 0), 46 runs.
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/252/climate+model+simulation+crashes.zip"
))
climate <- read.table(file.path(dir, "pop_failures.dat"), header = TRUE)
## Removed: Study and Run (study and run numbers of the design).
info$climate <- write_binary_csv(
  climate[!names(climate) %in% c("Study", "Run", "outcome")],
  as.numeric(climate$outcome == 0), "climate_model_crashes"
)

## ---------------------------------------------------------------------------
## LSVT Voice Rehabilitation
## https://archive.ics.uci.edu/dataset/282/lsvt+voice+rehabilitation
## 126 voice recordings, 310 numeric dysphonia measures (sheet "Data"). The
## recordings come from 14 people (9 each), so a random split puts the same
## person in both halves. The "Subject demographics" sheet is not used.
## target = 1: unacceptable phonation (class 2 on sheet "Binary response").
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/282/lsvt+voice+rehabilitation.zip"
))
lsvt_file <- file.path(dir, "LSVT_voice_rehabilitation.xlsx")
lsvt_x <- as.data.frame(readxl::read_excel(lsvt_file, sheet = "Data"))
lsvt_y <- readxl::read_excel(lsvt_file, sheet = "Binary response")[[1]]
info$lsvt <- write_binary_csv(
  lsvt_x, as.numeric(lsvt_y == 2), "lsvt_voice_rehabilitation"
)

## ---------------------------------------------------------------------------
## Arcene (NIPS 2003 feature selection challenge)
## https://archive.ics.uci.edu/dataset/167/arcene
## Mass-spectrometry data, 10000 features (7000 real, 3000 random "probes").
## Only the training (100) and validation (100) sets have labels; they are
## combined (200 rows). The unlabeled test set is not used. Labels are -1/1.
## target = 1: cancer (label 1).
## ---------------------------------------------------------------------------
dir <- unzip_uci(download_uci(
  "https://archive.ics.uci.edu/static/public/167/arcene.zip"
))
arcene_x <- rbind(
  read.table(file.path(dir, "ARCENE", "arcene_train.data")),
  read.table(file.path(dir, "ARCENE", "arcene_valid.data"))
)
names(arcene_x) <- paste0("f", seq_len(ncol(arcene_x)))
arcene_y <- c(
  scan(file.path(dir, "ARCENE", "arcene_train.labels"), quiet = TRUE),
  scan(file.path(dir, "arcene_valid.labels"), quiet = TRUE)
)
## 39 features are constant over the 200 rows and are dropped.
info$arcene <- write_binary_csv(
  arcene_x, as.numeric(arcene_y == 1), "arcene"
)

## ---------------------------------------------------------------------------
## Summary of the written files
## ---------------------------------------------------------------------------
uci_summary <- do.call(rbind, info)
rownames(uci_summary) <- NULL
print(uci_summary)
