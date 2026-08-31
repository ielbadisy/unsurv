# =============================================================================
# data-raw/metabric_survdnn_curves.R
#
# Regenerates inst/extdata/metabric_survdnn_curves.rds, the pre-computed set of
# individualized survival curves used by the package vignette
# (vignettes/unsurv-intro.Rmd).
#
# PROVENANCE
# ----------
# These curves are the exact worked example from the unsurv application note:
#
#   El Badisy I (2026). "unsurv: clustering individualized survival curves."
#   Bioinformatics Advances, 6(1), vbag218. <doi:10.1093/bioadv/vbag218>
#
# Pipeline (identical to the paper's scripts/metabric_survdnn_unsurv.R):
#   1. METABRIC cohort from the 'biostatlab' package (data(metabric)); overall
#      survival in months as the outcome.
#   2. Predictors: clinical + molecular numeric variables with >= 90% complete
#      cases, restricted to the 120 highest-variance columns, complete cases only.
#   3. A single 60/20/20 split (seed 20260615): train / partition / validation.
#   4. A deep AFT survival network (survdnn, loss = "aft") fitted on 'train'
#      only, with the paper's architecture and optimiser settings.
#   5. Predicted survival probabilities on a 100-point time grid spanning the
#      observed follow-up (min positive event time -> 90th percentile of the
#      training survival times), for the partition and validation patients.
#
# The saved object therefore contains NO METABRIC covariates or identifiers,
# only model-predicted survival probabilities plus the (time, status) pair for
# each patient. This keeps the vignette self-contained and lets it reproduce
# the paper figure without importing survdnn / biostatlab / torch at build time.
#
# This script is NOT part of the build (data-raw/ is in .Rbuildignore). Run it
# from the package root when the paper pipeline changes:
#   Rscript data-raw/metabric_survdnn_curves.R
# =============================================================================

suppressPackageStartupMessages({
  library(biostatlab)   # METABRIC data
  library(survival)
  library(survdnn)      # deep AFT survival model + predict()
})

set.seed(20260615)

extract_survmat <- function(pred) {
  if (is.matrix(pred)) return(pred)
  if (is.data.frame(pred)) return(as.matrix(pred))
  if (is.list(pred) && !is.null(pred$survival)) return(pred$survival)
  stop("predict() must return a matrix, data.frame, or list with $survival.")
}

## --- 1. METABRIC cohort --------------------------------------------------------

data(metabric, package = "biostatlab")
raw <- metabric
raw$os_time  <- raw$overall_survival_months
raw$os_event <- raw$overall_survival

exclude <- c(
  "patient_id", "overall_survival_months", "overall_survival",
  "os_time", "os_event", "death_from_cancer"
)

## --- 2. Predictor matrix -----------------------------------------------------

numeric_cols   <- names(raw)[vapply(raw, is.numeric, logical(1))]
numeric_cols   <- setdiff(numeric_cols, exclude)
complete_rate  <- vapply(raw[numeric_cols], function(x) mean(!is.na(x)), numeric(1))
numeric_cols   <- numeric_cols[complete_rate >= 0.90]
vars           <- vapply(raw[numeric_cols], stats::var, numeric(1), na.rm = TRUE)
vars           <- vars[is.finite(vars) & vars > 0]
selected       <- names(sort(vars, decreasing = TRUE))[seq_len(min(120L, length(vars)))]

model_df <- raw[, c("os_time", "os_event", selected), drop = FALSE]
model_df <- model_df[stats::complete.cases(model_df), , drop = FALSE]
names(model_df) <- make.names(names(model_df), unique = TRUE)

## --- 3. 60 / 20 / 20 split: train / partition / validation ------------------

n           <- nrow(model_df)
idx_all     <- sample.int(n)
n_train     <- floor(0.60 * n)
n_partition <- floor(0.20 * n)
idx_train      <- idx_all[seq_len(n_train)]
idx_partition  <- idx_all[n_train + seq_len(n_partition)]
idx_validation <- idx_all[(n_train + n_partition + 1):n]

train      <- model_df[idx_train, , drop = FALSE]
partition  <- model_df[idx_partition, , drop = FALSE]
validation <- model_df[idx_validation, , drop = FALSE]

## --- 4. Deep AFT survival network (fitted on train only) --------------------

fml <- stats::as.formula("Surv(os_time, os_event) ~ .")

mod <- survdnn(
  formula     = fml,
  data        = train,
  hidden      = c(64L, 32L),
  activation  = "relu",
  dropout     = 0.10,
  batch_norm  = TRUE,
  epochs      = 120L,
  lr          = 5e-4,
  loss        = "aft",
  optimizer   = "adamw",
  optim_args  = list(weight_decay = 1e-5),
  .seed       = 20260615,
  verbose     = FALSE
)

## --- 5. Predicted survival curves on a shared time grid --------------------

t_max <- as.numeric(stats::quantile(train$os_time, probs = 0.90, na.rm = TRUE))
times <- seq(
  max(0.1, min(train$os_time[train$os_time > 0], na.rm = TRUE)),
  t_max,
  length.out = 100
)

S_partition  <- extract_survmat(predict(mod, newdata = partition,  type = "survival", times = times))
S_validation <- extract_survmat(predict(mod, newdata = validation, type = "survival", times = times))

storage.mode(S_partition)  <- "double"
storage.mode(S_validation) <- "double"
dimnames(S_partition)  <- NULL
dimnames(S_validation) <- NULL

# Round to 5 decimals: keeps the curves visually and numerically faithful to
# the paper while letting xz compression bring the shipped file well under the
# CRAN data-size guideline.
S_partition  <- round(S_partition,  5)
S_validation <- round(S_validation, 5)

metabric_survdnn_curves <- list(
  times        = times,
  S_partition  = S_partition,
  S_validation = S_validation,
  os_time      = list(partition = partition$os_time,  validation = validation$os_time),
  os_event     = list(partition = partition$os_event, validation = validation$os_event),
  provenance   = paste0(
    "survdnn (AFT) predicted survival curves for the METABRIC cohort, ",
    "reproducing the worked example in El Badisy (2026), Bioinformatics ",
    "Advances 6(1), vbag218, doi:10.1093/bioadv/vbag218. Generated by ",
    "data-raw/metabric_survdnn_curves.R with seed 20260615. Contains only ",
    "model-predicted probabilities and (time, status); no METABRIC ",
    "covariates or patient identifiers."
  )
)

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(
  metabric_survdnn_curves,
  file = "inst/extdata/metabric_survdnn_curves.rds",
  version = 2,
  compress = "xz"
)

cat(sprintf(
  "Wrote inst/extdata/metabric_survdnn_curves.rds\n  partition: %d x %d   validation: %d x %d   grid: [%.1f, %.1f]\n",
  nrow(S_partition), ncol(S_partition),
  nrow(S_validation), ncol(S_validation),
  min(times), max(times)
))
