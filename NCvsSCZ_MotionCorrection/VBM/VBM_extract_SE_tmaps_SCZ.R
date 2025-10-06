# ===============================================================
# IMAGE RESLICING + SE MAP COMPUTATION (NIfTI, R)
# ---------------------------------------------------------------
# This script performs:
#   1) Reslicing: aligns input NIfTI masks/images to a common reference 
#      using ANTsR (e.g., grey-matter template or group mask).
#   2) SE map computation: computes a voxelwise standard error (SE) map
#      from a t-map and its associated contrast (beta) map using:
#         SE = beta / t
#
# Notes:
#   - All inputs must be aligned in space.
#   - Masking and geometry checks ensure consistency.
#
# Author: Roberta Passiatore
# Updated: 2025-08-29
# ===============================================================

suppressPackageStartupMessages({
  need <- c("RNifti", "ANTsR")
  to_install <- need[!sapply(need, requireNamespace, quietly = TRUE)]
  if (length(to_install)) install.packages(to_install, dependencies = TRUE)
  library(RNifti)
  library(ANTsR)
})

# ---------- Helper functions ----------
stop_if_mismatch <- function(a, b, what_a = "image A", what_b = "image B") {
  if (!identical(dim(a), dim(b))) {
    stop("Geometry mismatch between ", what_a, " and ", what_b, ".")
  }
}

# ===============================================================
# 1) RESLICE ALL NIFTI FILES IN A FOLDER TO A REFERENCE (ANTsR)
# ===============================================================
reslice_cfg <- list(
  input_dir      = "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/con",
  output_prefix  = "r",
  reference_path = "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/rGM_cat12_bin.nii.gz",
  pattern        = "nii.gz"
)

RUN_RESLICING <- TRUE

if (RUN_RESLICING) {
  message("== Reslicing images to reference ==")
  reference <- antsImageRead(reslice_cfg$reference_path, dimension = 3)
  mask_files <- list.files(reslice_cfg$input_dir, pattern = reslice_cfg$pattern, full.names = TRUE)

  for (mask_path in mask_files) {
    mask <- antsImageRead(mask_path, dimension = 3)
    mask_resliced <- resampleImageToTarget(mask, reference, interpType = 1) # linear interpolation
    file_name <- basename(mask_path)
    output_path <- file.path(reslice_cfg$input_dir, paste0(reslice_cfg$output_prefix, file_name))
    antsImageWrite(mask_resliced, output_path)
    message("  ✓ Resliced saved -> ", output_path)
  }
}

# ===============================================================
# 2) COMPUTE SE MAP FROM T-MAP AND BETA (CON) MAP: SE = beta / t
# ===============================================================
se_cfg <- list(
  tmap_path   = "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/tmaps/rUNIBA2_FI_PCA_spmT_0004_maskGM.nii.gz", # alterantive noPCA
  beta_path   = "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/con/UNIBA2_FI_PCA_con.nii.gz", # alterantive noPCA
  mask_path   = "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/rGM_cat12_bin.nii.gz",
  out_se_path = "UNIBA2_FI_PCA_SE.nii.gz", # alterantive noPCA
  eps         = 1e-6,  # for safe division
  se_cap_abs  = 1.0    # cap |SE| to avoid extreme/unstable values
)

RUN_SE_MAP <- TRUE

if (RUN_SE_MAP) {
  message("== Computing SE map from t and beta ==")
  i_t   <- readNifti(se_cfg$tmap_path)
  i_b   <- readNifti(se_cfg$beta_path)
  i_msk <- readNifti(se_cfg$mask_path)

  stop_if_mismatch(i_t, i_b, "t-map", "beta-map")
  stop_if_mismatch(i_t, i_msk, "t-map", "mask")

  mask_bin <- i_msk > 0
  t_m <- i_t * mask_bin
  b_m <- i_b * mask_bin

  denom <- t_m
  denom[abs(denom) < se_cfg$eps] <- se_cfg$eps
  se_map <- b_m / denom

  # Clean and constrain SE values
  se_map[!mask_bin] <- NA_real_
  se_map[abs(se_map) > se_cfg$se_cap_abs & !is.na(se_map)] <- se_cfg$se_cap_abs
  se_map <- abs(se_map)

  writeNifti(asNifti(se_map, reference = i_t), se_cfg$out_se_path)
  message("  ✓ SE map written -> ", se_cfg$out_se_path)
}

# ------------------ Reference formulas ------------------
# F = (R²/(p-1)) / ((1-R²)/(n-p))
# R² = F / (F + (n-p)/(p-1))
# Adjusted R² = 1 - [(1-R²) * (n-1) / (n-p-1)]
