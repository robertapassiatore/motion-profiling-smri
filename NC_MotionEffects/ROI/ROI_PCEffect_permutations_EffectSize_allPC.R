# ===============================================================
# Motion Profiling: Variance Explained (R²) Analysis by Motion PCs
# ---------------------------------------------------------------
# This script performs linear regression analyses to estimate the
# proportion of variance (R²) in brain ROI measures explained by
# five motion principal components (PC1–PC5) across multiple datasets.
#
# For each dataset (.csv):
#   - Residualizes ROI data for covariates (sex, age, age², TIV)
#   - Standardizes numeric variables
#   - Fits linear models of each ROI ~ PC1 + ... + PC5
#   - Computes adjusted R², p-values, and permutation-based empirical p-values
#   - Performs FDR correction
#
# Outputs:
#   - One CSV per dataset per analysis (cortical / subcortical ROIs)
#     containing regression statistics and permutation results.
#
# Author: Roberta Passiatore
# Updated: 2025-08-28
# ===============================================================

# ---- Packages ----
req <- c("readr","dplyr","parallel","ggplot2","lme4","multivarious")
to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(req, library, character.only = TRUE)

# ---- Config ----
input_dir  <- "/documents/Roberta/MotionProfiling/NC_LinearMotionEffects"
output_dir <- "/documents/Roberta/MotionProfiling/NC_LinearMotionEffects/results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)
n_perm <- 5000
use_parallel <- TRUE
n_cores <- max(1, parallel::detectCores() - 1)

pc_cols <- c("PC1","PC2","PC3","PC4","PC5")
covars  <- c("sex","age","agesq","TIV")
roi_cortical_idx    <- 12:111
roi_subcortical_idx <- 112:119

# ---- Helpers ----
get_dataset_name <- function(file) sub("\\.csv$", "", basename(file))

lmp <- function(model) {
  if (!inherits(model, "lm")) stop("Not an 'lm' object.")
  f <- summary(model)$fstatistic
  as.numeric(pf(f[1], f[2], f[3], lower.tail = FALSE))
}

perm_adjR2 <- function(y, X, n_perm, parallel = TRUE, n_cores = 1) {
  n <- length(y)
  p <- ncol(X) + 1
  fit_adjR2 <- function(y, X) {
    m <- lm.fit(x = cbind(1, as.matrix(X)), y = y)
    rss <- sum(residuals.lm(m)^2)
    tss <- sum((y - mean(y))^2)
    r2 <- 1 - rss / tss
    1 - (1 - r2) * (n - 1) / (n - p)
  }
  worker <- function(i) {
    Xperm <- as.data.frame(lapply(X, sample))
    fit_adjR2(y, Xperm)
  }
  if (!parallel) vapply(seq_len(n_perm), worker, numeric(1))
  else if (.Platform$OS.type == "unix")
    parallel::mclapply(seq_len(n_perm), worker, mc.cores = n_cores) |> unlist()
  else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::parLapply(cl, seq_len(n_perm), worker) |> unlist()
  }
}

run_analysis <- function(df, dataset_name, roi_idx, suffix) {
  out <- vector("list", length(roi_idx))
  Xpcs <- df[, pc_cols, drop = FALSE]
  for (ii in seq_along(roi_idx)) {
    col_id <- roi_idx[ii]
    roi_name <- colnames(df)[col_id]
    y <- df[[col_id]]
    form <- as.formula(paste(roi_name, "~", paste(pc_cols, collapse = " + ")))
    model <- lm(form, data = df)
    smry  <- summary(model)
    p_val <- lmp(model)
    adjR  <- smry$adj.r.squared
    se    <- smry$sigma
    perm_Rs <- perm_adjR2(y, Xpcs, n_perm, parallel = use_parallel, n_cores = n_cores)
    emp_p <- mean(abs(perm_Rs) >= abs(adjR))
    avg_R <- mean(perm_Rs)
    out[[ii]] <- data.frame(
      Dataset = dataset_name,
      ROI = roi_name,
      Analysis = "allcov",
      p_value = p_val,
      R = adjR,
      SE = se,
      Empirical_p_value = emp_p,
      Avg_R = avg_R,
      stringsAsFactors = FALSE
    )
  }
  res <- dplyr::bind_rows(out)
  res$pFDR <- p.adjust(res$p_value, method = "fdr")
  fname <- file.path(output_dir, sprintf("regression_results_allPC_permutations_%s_%s.csv", suffix, dataset_name))
  readr::write_csv(res, fname)
  message("Wrote: ", fname)
  invisible(res)
}

# ---- Main Loop ----
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No CSV files found in: ", input_dir)

for (f in files) {
  message("Processing: ", basename(f))
  data <- readr::read_csv(f, show_col_types = FALSE)
  dataset_name <- get_dataset_name(f)
  for (fac in c("group","sex","ID")) {
    if (fac %in% names(data)) data[[fac]] <- as.factor(data[[fac]])
  }
  data_std <- data
  roi_all_idx <- unique(c(roi_cortical_idx, roi_subcortical_idx))
  roi_all_idx <- roi_all_idx[roi_all_idx <= ncol(data_std)]
  have_covs <- all(covars %in% names(data_std))
  if (have_covs && length(roi_all_idx)) {
    data_std[, roi_all_idx] <- multivarious::residualize(
      stats::as.formula(paste("~", paste(covars, collapse = " + "))),
      data_std[, roi_all_idx, drop = FALSE],
      data_std
    )
  }
  num_cols <- vapply(data_std, is.numeric, logical(1))
  data_std[num_cols] <- lapply(data_std[num_cols], scale)_]()]()_
