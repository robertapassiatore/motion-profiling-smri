# ===============================================================
# GROUP DIFFERENCES ON ROI MEASURES (NC Percentiles)
# ---------------------------------------------------------------
# For each dataset (.csv), this script:
#   1) standardizes numeric columns,
#   2) fits linear models ROI ~ group + covariates,
#   3) (optionally) permutes PCs to get empirical p-values,
#   4) saves per-ROI stats (beta, t, SE, p, Cohen's d, pFDR).
#
# Outputs: results/regression_results_permutations_{corticalROI|subcorticalROI}_{DATASET}.csv
#
# Author: Roberta Passiatore | Updated: 2025-08-29
# ===============================================================

suppressPackageStartupMessages({
  req <- c("readr","dplyr","parallel","ggplot2")
  to_install <- req[!sapply(req, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  lapply(req, library, character.only = TRUE)
})

# -------------------- CONFIG --------------------
input_dir <- "~/data"
output_dir <- file.path(input_dir, "results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

roi_cortical_idx    <- 12:111
roi_subcortical_idx <- 112:119

# Empirical p-value settings (applies only to the PC-corrected model)
use_permutations <- TRUE
n_perm <- 500  # set to 5000 for final runs

pc_cols <- c("PC1","PC2","PC3","PC4","PC5")
covars  <- c("sex","age","agesq","TIV") #Site if needed
group_col <- "group"   # expects levels like "1","2" (reference = first level)

set.seed(1234)

# ------------------ HELPERS --------------------
get_dataset_name <- function(file) sub("\\.csv$", "", basename(file))

# Extract the first coefficient row that corresponds to the group effect
get_group_coef_row <- function(coef_mat, group_col = "group") {
  rn <- rownames(coef_mat)
  hit <- grep(paste0("^", group_col), rn, value = TRUE)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

cohens_d_from_lm <- function(model) {
  beta <- coef(model)
  sres <- sd(residuals(model))
  if (!is.finite(sres) || sres == 0) return(NA_real_)
  # pick the group coefficient (first match)
  cm <- coef(summary(model))
  row <- get_group_coef_row(cm)
  if (is.na(row)) return(NA_real_)
  as.numeric(cm[row, "Estimate"]) / sres
}

fit_and_extract <- function(formula, data, group_col = "group") {
  model <- lm(formula, data = data)
  smry  <- summary(model)
  cm    <- coef(smry)
  row   <- get_group_coef_row(cm, group_col)
  if (is.na(row)) {
    return(list(
      beta = NA_real_, p = NA_real_, t = NA_real_, se = NA_real_,
      d = NA_real_, model = model
    ))
  }
  list(
    beta = as.numeric(cm[row, "Estimate"]),
    p    = as.numeric(cm[row, "Pr(>|t|)"]),
    t    = as.numeric(cm[row, "t value"]),
    se   = as.numeric(cm[row, "Std. Error"]),
    d    = cohens_d_from_lm(model),
    model = model
  )
}

empirical_p_for_group <- function(formula, data, pc_cols, observed_beta, n_perm = 500) {
  if (!all(pc_cols %in% names(data))) return(NA_real_)
  perm_b <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    dsh <- data
    dsh[pc_cols] <- lapply(dsh[pc_cols], sample, replace = FALSE)
    fit <- lm(formula, data = dsh)
    cm  <- coef(summary(fit))
    row <- get_group_coef_row(cm)
    perm_b[i] <- if (!is.na(row)) as.numeric(cm[row, "Estimate"]) else NA_real_
  }
  perm_b <- perm_b[is.finite(perm_b)]
  if (!length(perm_b)) return(NA_real_)
  mean(abs(perm_b) >= abs(observed_beta))
}

run_analysis <- function(df, dataset_name, roi_idx, suffix) {
  res_list <- vector("list", length(roi_idx))

  for (ii in seq_along(roi_idx)) {
    col_id <- roi_idx[ii]
    if (col_id > ncol(df)) next
    roi_name <- colnames(df)[col_id]

    # Two models
    form_std  <- as.formula(paste(roi_name, "~", paste(c(group_col, covars), collapse = " + ")))
    form_corr <- as.formula(paste(roi_name, "~", paste(c(group_col, covars, pc_cols), collapse = " + ")))

    # Standard model
    std <- fit_and_extract(form_std, df, group_col)

    # Corrected model + empirical p via PC permutations
    corr <- fit_and_extract(form_corr, df, group_col)
    emp_p <- NA_real_; avg_d <- NA_real_
    if (use_permutations && is.finite(corr$beta)) {
      emp_p <- empirical_p_for_group(form_corr, df, pc_cols, corr$beta, n_perm)
      # quick average Cohen's d across permutations (optional, mirrors original intent)
      # compute via perm betas divided by residual SD of corrected model (fixed)
      # For stability and speed, skip unless you truly need it.
      avg_d <- NA_real_
    }

    # Keep only the "standard" line in the final table (to match original behavior)
    # If you wish to keep both, add an extra row here for "corrected".
    res_list[[ii]] <- data.frame(
      Dataset = dataset_name,
      ROI = roi_name,
      Analysis = "standard",
      p_value = std$p,
      Beta_Coefficient_Group = std$beta,
      Empirical_p_value = emp_p,       # from corrected model permutations (if computed)
      Cohens_d = std$d,
      T_value = std$t,
      SE = std$se,
      Avg_Cohens_d = avg_d,
      stringsAsFactors = FALSE
    )
  }

  res <- dplyr::bind_rows(res_list)
  res$pFDR <- p.adjust(res$p_value, method = "fdr")
  res$N <- nrow(df)

  out_file <- file.path(output_dir,
                        sprintf("regression_results_permutations_%s_%s.csv", suffix, dataset_name))
  readr::write_csv(res, out_file)
  message("Wrote: ", out_file)
  invisible(res)
}

# ----------------- MAIN LOOP -----------------
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (!length(files)) stop("No CSV files found in: ", input_dir)

for (f in files) {
  message("Processing: ", basename(f))
  data <- readr::read_csv(f, show_col_types = FALSE)

  # Drop first column if it's an index
  if (names(data)[1] %in% c("...1","X1","index","Unnamed: 0")) data <- data[, -1]

  dataset_name <- get_dataset_name(f)

  # Factors
  for (fac in c(group_col, "sex")) {
    if (fac %in% names(data)) data[[fac]] <- as.factor(data[[fac]])
  }

  # Standardize numeric columns (z-score)
  num_cols <- vapply(data, is.numeric, logical(1))
  data_std <- data
  data_std[num_cols] <- lapply(data_std[num_cols], scale)

  # Run analyses
  cort_idx <- roi_cortical_idx[roi_cortical_idx <= ncol(data_std)]
  if (length(cort_idx)) run_analysis(data_std, dataset_name, cort_idx, "corticalROI")

  sub_idx <- roi_subcortical_idx[roi_subcortical_idx <= ncol(data_std)]
  if (length(sub_idx)) run_analysis(data_std, dataset_name, sub_idx, "subcorticalROI")
}

# ----------------- File rename helper -----------------
# Rename any output files containing "SCAP" to "WM" in the results folder
files_to_rename <- list.files(path = output_dir, pattern = "SCAP", full.names = TRUE)
for (old_file in files_to_rename) {
  new_file <- gsub("SCAP", "WM", old_file)
  if (file.rename(old_file, new_file)) {
    cat("Renamed:", old_file, "->", new_file, "\n")
  } else {
    cat("Failed to rename:", old_file, "\n")
  }
}
