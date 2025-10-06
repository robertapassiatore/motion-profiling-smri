# ===============================================================
# SESSION-LEVEL FISHER P + MEAN R (CORTICAL ROI OUTPUT FORMAT)
# ---------------------------------------------------------------
# Aggregates per-dataset cortical ROI results by session (RS, FI, WM):
#   - Fisher-combined p-values across datasets (per ROI, per session)
#   - Mean R across datasets (per ROI, per session)
# Writes CSV exactly matching expected columns:
#   ROI, meanR, pval, session
#
# Author: Roberta Passiatore
# Updated: 2025-08-28
# ===============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
})

library(readr)
library(dplyr)

# ---- Config ----
res.dir        <- "~/Documents/Motion-Profiling/analyses/NC_LinearMotionEffects/results/"
pattern_kind   <- "corticalROI"                    # match your expected file
sessions_keep  <- c("RS","FI","WM")                # sessions to include
use_empirical  <- FALSE                            # TRUE => combine Empirical_p_value instead of p_value
out_file       <- file.path(res.dir, "res.LinearEffect_allPC_corticalROI_metap.csv")

# ---- Helpers ----
fisher_combine <- function(p) {
  p[!is.finite(p) | is.na(p) | p <= 0] <- .Machine$double.xmin
  p[p > 1] <- 1
  k <- length(p)
  X <- -2 * sum(log(p))
  stats::pchisq(X, df = 2 * k, lower.tail = FALSE)
}

# ---- Load all cortical files ----
files <- list.files(
  res.dir,
  pattern = paste0("^regression_results_allPC_permutations_", pattern_kind, "_.*\\.csv$"),
  full.names = TRUE
)

if (!length(files)) stop("No input files found in: ", res.dir)

all_data <- dplyr::bind_rows(lapply(files, function(fp) {
  dat <- readr::read_csv(fp, show_col_types = FALSE)
  base <- basename(fp)
  stem <- gsub(paste0("^regression_results_allPC_permutations_", pattern_kind, "_|\\.csv$"), "", base)
  parts <- unlist(strsplit(stem, "_"))
  dat$Site    <- ifelse(length(parts) >= 1, parts[1], NA_character_)
  dat$Session <- ifelse(length(parts) >= 2, parts[2], NA_character_)
  # keep only the full model if present
  if ("Analysis" %in% names(dat)) dat <- dplyr::filter(dat, Analysis == "allcov")
  dat
}))

# Choose which p to combine
p_col <- if (use_empirical && "Empirical_p_value" %in% names(all_data)) "Empirical_p_value" else "p_value"
if (!p_col %in% names(all_data)) {
  stop("Column '", p_col, "' not found. Available: ", paste(names(all_data), collapse = ", "))
}

# Keep desired sessions
if (!"Session" %in% names(all_data)) stop("Session column not present in input files.")
all_data <- all_data %>% filter(Session %in% sessions_keep)

# ---- Aggregate: per session x ROI ----
summary_df <- all_data %>%
  group_by(Session, ROI) %>%
  summarise(
    k_datasets = sum(is.finite(.data[[p_col]]) & !is.na(.data[[p_col]])),
    pval       = if (k_datasets > 0) fisher_combine(.data[[p_col]]) else NA_real_,
    meanR      = mean(R, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Reorder & rename to match expected columns exactly ----
summary_df <- summary_df %>%
  select(ROI, meanR, pval, session = Session) %>%
  arrange(session, ROI)

# Add 'Unnamed: 0' index column to match expected output
summary_df <- summary_df %>%
  tibble::rowid_to_column(var = "Unnamed: 0") %>%
  mutate(`Unnamed: 0` = `Unnamed: 0` - 1)  # start at 0 to mimic typical index

# ---- Write ----
readr::write_csv(summary_df, out_file)
cat("Wrote: ", out_file, "\n")
