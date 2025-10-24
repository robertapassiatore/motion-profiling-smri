# ===============================================================
# VBM_PVALUE_SUMUP.R
# ---------------------------------------------------------------
# This script combines voxel-wise p-value maps from multiple datasets
# using Fisher’s method to obtain a single, whole-brain significance map.
# It also applies FDR (q<0.05) and uncorrected (p<0.001) thresholds to 
# the mean R² map, generating thresholded R² images.
#
# Inputs:
#   - Voxelwise p-value maps from individual datasets
#   - Grey matter binary mask (to define analysis voxels)
#   - Mean R² map
#
# Outputs:
#   - Fisher-combined p-map (whole brain)
#   - R² map thresholded by FDR < 0.05
#   - R² map thresholded by uncorrected p < 0.001
#
# Author: Roberta Passiatore
# Updated: 2025-08-26
# ===============================================================

# --- Setup ---
library(RNifti)

# ---- Inputs ----
p_paths <- c(
  "~/data/UNIBA2/RS/pval_0001.nii",
  "~/data/HCP/RS/pval_0001.nii",
  "~/data/DNS/RS/pval_0001._

# --- Path to binary mask ---
mask_path <- "~/Data/rrGM_cat12_bin.nii"
r2_path = "~/Data/mean_r2_0001_RS.nii"

# --- Load mask and get voxel indices ---
mask <- readNifti(mask_path)
mask_idx <- which(mask > 0)

# --- Load p-value maps ---
p_maps <- lapply(p_paths, readNifti)

# --- Extract p-values inside the mask ---
p_vals_list <- lapply(p_maps, function(p) {
  p_values <- p[mask_idx]
  # Avoid log(0) by replacing zero or negative p-values with smallest non-zero float
  p_values[p_values <= 0] <- .Machine$double.eps
  return(p_values)
})

# --- Apply Fisher's method ---
k <- length(p_vals_list) - 1 
logsum <- Reduce(`+`, lapply(p_vals_list, function(p) -2 * log(p)))
fisher_stat <- logsum
fisher_pvals <- pchisq(fisher_stat, df = 1 * k, lower.tail = FALSE)

# --- Create output array ---
combined_map <- array(NA_real_, dim = dim(mask))
combined_map[mask_idx] <- fisher_pvals

# --- Save output map ---
combined_nifti <- asNifti(combined_map, reference = p_maps[[1]])
writeNifti(combined_nifti, "~/results/fisher_pval_0001_RS.nii.gz")

