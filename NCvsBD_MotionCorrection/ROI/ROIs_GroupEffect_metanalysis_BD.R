# ===============================================================
# META-ANALYSIS OF GROUP EFFECTS ON ROI R² VALUES
# ---------------------------------------------------------------
# This script performs a meta-analysis across multiple datasets
# to combine the effects (R values converted to Fisher’s z-scores)
# of motion principal components (PCs) on cortical ROIs.
#
# Steps:
#   1. Load all per-dataset regression outputs
#   2. Convert R to Fisher’s z-scores
#   3. Organize data by site and session
#   4. Perform a random-effects meta-analysis (via metafor)
#      for each ROI and session
#   5. Export combined meta-analytic results
#
# Author: Roberta Passiatore
# Updated: 2025-08-28
# ===============================================================

# --- Required libraries ---
if (!requireNamespace("metafor", quietly = TRUE)) install.packages("metafor")
library(metafor)

# --- Helper functions ---

# Convert t-statistic to partial r
t2r <- function(t, df) {
  sign(t) * sqrt(t^2 / (t^2 + df))
}

# Fisher’s r-to-z transformation
r2z <- function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# --- Setup ---
res.dir <- "~/Documents/Motion-Profiling/analyses/NC_GroupEffects_NCvsBD/results/"
setwd(res.dir)
datasets <- list.files(pattern = "NCvsBD_permutations_corticalROI.*\\.csv$")

if (length(datasets) == 0)
  stop("No dataset files found in ", res.dir)

# --- Load and combine datasets ---
all.data <- do.call(rbind, lapply(datasets, function(dataset) {
  data <- read.csv(dataset)
  dataset_name <- gsub("^regression_results_permutations_corticalROI|\\.csv$", "", basename(dataset))
  parts <- unlist(strsplit(dataset_name, "_"))
  
  # Convert R to Fisher's z
  data$zscore <- r2z(data$R)
  
  data$Site    <- ifelse(length(parts) >= 1, parts[1], NA)
  data$Session <- ifelse(length(parts) >= 2, parts[2], NA)
  
  data
}))

# --- Initialize meta-analysis structure ---
meta.data.res <- data.frame()

# Get ROI list and sessions
ROIs <- unique(all.data$ROI)
sessions <- unique(all.data$Session)
sessions <- sessions[1:3]  # using first three sessions (no FE)

# --- Meta-analysis loop ---
for (sess in sessions) {
  message("Processing session: ", sess)
  data.session <- subset(all.data, Session == sess)
  
  # Initialize per-session result frame
  meta.session <- data.frame(
    ROI = ROIs,
    zval = NA, ci.lb = NA, ci.ub = NA, pval = NA,
    I2 = NA, H2 = NA, QE = NA, QEp = NA,
    session = sess
  )
  
  for (r in seq_along(ROIs)) {
    data.sub <- subset(data.session, ROI == ROIs[r])
    if (nrow(data.sub) == 0) next
    
    data.test <- data.frame(
      study = data.sub$Site,
      yi = data.sub$zscore,
      vi = data.sub$SE^2
    )
    
    # Perform random-effects meta-analysis
    res <- tryCatch(
      metafor::rma(yi, vi, data = data.test, method = "ML"),
      error = function(e) NULL
    )
    
    if (!is.null(res)) {
      meta.session[r, c("zval", "ci.lb", "ci.ub", "pval", "I2", "H2", "QE", "QEp")] <- 
        c(res$zval, res$ci.lb, res$ci.ub, res$pval, res$I2, res$H2, res$QE, res$QEp)
    }
  }
  
  # Append to combined results
  meta.data.res <- rbind(meta.data.res, meta.session)
}

# --- Save results ---
out_file <- file.path(res.dir, "GroupEffect_NCvsBD_corticalROI_meta_analysis.csv")
write.csv(meta.data.res, out_file, row.names = FALSE)

cat("Meta-analysis completed.\nResults saved to:", out_file, "\n")
