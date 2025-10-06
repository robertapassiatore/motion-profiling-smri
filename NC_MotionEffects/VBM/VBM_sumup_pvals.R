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
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/UNIBA2/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/HCP/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/DNS/FI/pval_0001._
