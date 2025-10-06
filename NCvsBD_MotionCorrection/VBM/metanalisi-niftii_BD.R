# ===============================================================
# VOXEL-WISE RANDOM-EFFECTS META-ANALYSIS OF GROUP T-MAPS (NC vs BD)
# ---------------------------------------------------------------
# Purpose:
#   Combine multiple site-level group-effect t-maps (and their SE maps)
#   using a voxel-wise random-effects model (metafor::rma.mv) to produce:
#     - meta-analytic z, p, beta maps
#     - FDR-corrected p maps and a thresholded version
#   Optionally, extract cluster peaks from a thresholded map.
#
# Inputs (example layout):
#   - T-maps:   /.../tmaps_GroupEffect/BD/tmaps/r*_PCA*.nii.gz # alternative file noPCA
#   - SE maps:  /.../tmaps_GroupEffect/BD/se/r*_PCA*.nii.gz # alternative file noPCA
#   - Mask:     /.../tmaps_GroupEffect/rGM_cat12_bin.nii.gz
#   - Sample sizes file: res_BD.csv with columns: study, session, N
#
# Outputs (written to output_dir):
#   meta_z_*.nii.gz, meta_p_*.nii.gz, meta_pFDR_*.nii.gz, meta_pFDR_thr_*.nii.gz, meta_beta_*.nii.gz
#
# Notes:
#   - This is computationally heavy; consider running on a machine with ample RAM/CPU.
#   - All input images must share the same geometry (dim/affine).
# Author: Roberta Passiatore
# Updated: 2025-08-29
# ===============================================================

# ---- Packages (install if needed) ----
need <- c("oro.nifti", "metafor", "dplyr")
to_install <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(oro.nifti)
library(metafor)
library(dplyr)

# ---- Helper functions ----

# Convert t-statistic to partial r
t2r <- function(t, df) {
  sign(t) * sqrt(t^2 / (t^2 + df))
}

# Fisher’s r-to-z transformation
r2z <- function(r) 0.5 * log((1 + r) / (1 - r))

# Safe NIfTI read -> numeric array (replace NA with 0)
read_nifti_sanitize <- function(path, reorient = TRUE) {
  nii <- readNIfTI(path, reorient = reorient)
  arr <- nii@.Data
  arr[is.na(arr)] <- 0
  attr(arr, "header") <- nii
  arr
}

# Check images are the same geometry
check_geom_match <- function(ref_header, arr_list) {
  ref_dim <- dim(ref_header)
  ok <- all(vapply(arr_list, function(a) identical(dim(a), ref_dim), logical(1)))
  if (!ok) stop("Geometry mismatch across input images.")
}

# ---- Paths / IO ----
tmap_dir  <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/tmaps"
se_dir    <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/se"
mask_path <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/rGM_cat12_bin.nii.gz"
sample_csv <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/res_BD.csv"
output_dir <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD"

tmap_files <- list.files(tmap_dir, pattern = glob2rx("r*_PCA*.nii.gz"), full.names = TRUE)
se_files   <- list.files(se_dir,   pattern = glob2rx("r*_noPCA*.nii.gz"), full.names = TRUE)

if (length(tmap_files) == 0) stop("No t-map files found in: ", tmap_dir)
if (length(se_files)   == 0) stop("No SE files found in: ", se_dir)
if (length(tmap_files) != length(se_files)) {
  warning("Different counts for t-maps (", length(tmap_files), ") and SE maps (", length(se_files), "). Proceeding with min length.")
  k <- min(length(tmap_files), length(se_files))
  tmap_files <- tmap_files[seq_len(k)]
  se_files   <- se_files[seq_len(k)]
}

# ---- Load data ----
mask_data <- read_nifti_sanitize(mask_path, reorient = TRUE)
ref_hdr   <- attr(mask_data, "header")

# Load t-maps and SE maps as lists
t_list  <- lapply(tmap_files, read_nifti_sanitize)
se_list <- lapply(se_files,   read_nifti_sanitize)

# Geometry checks
check_geom_match(ref_hdr, t_list)
check_geom_match(ref_hdr, se_list)
if (!identical(dim(mask_data), dim(t_list[[1]]))) stop("Mask and images have different dimensions.")

dims <- dim(mask_data)
dim_x <- dims[1]; dim_y <- dims[2]; dim_z <- dims[3]

# Sample sizes / df per study
db <- read.csv(sample_csv)
if (!all(c("study","session","N") %in% names(db))) {
  stop("Sample CSV must contain columns: study, session, N")
}
# Example df computation; adjust p (parameters) if needed per-study:
# (Your original code used different offsets: 11, 10, etc. Keep that logic here if required.)
# For clarity, create a vector of dfs same length as the number of images:
df_vec <- rep(NA_integer_, length(t_list))
# Fill per your previous logic (edit as appropriate):
df_vec[1] <- db$N[1] - 11
df_vec[2] <- db$N[2] - 11
df_vec[3] <- db$N[3] - 11
df_vec[4] <- db$N[4] - 10
df_vec[5] <- db$N[5] - 10
df_vec[6] <- db$N[6] - 10
# If you have more images and know their effective dfs, set them here:
if (length(df_vec) > 6) {
  # Example default for remaining: N - 10
  for (i in 7:length(df_vec)) df_vec[i] <- db$N[i] - 10
}

# ---- Allocate outputs ----
meta_z    <- array(NA_real_, dim = dims)
meta_p    <- array(NA_real_, dim = dims)
meta_beta <- array(NA_real_, dim = dims)

# Build mask indices once
mask_idx <- which(mask_data > 0, arr.ind = TRUE)

cat("Starting voxel-wise random-effects meta-analysis on", nrow(mask_idx), "voxels...\n")

# ---- Main voxel loop (compute-intensive) ----
for (ii in seq_len(nrow(mask_idx))) {
  x <- mask_idx[ii, 1]; y <- mask_idx[ii, 2]; z <- mask_idx[ii, 3]

  # Gather t and SE across studies at (x,y,z)
  t_vals  <- vapply(t_list,  function(a) a[x, y, z], numeric(1))
  se_vals <- vapply(se_list, function(a) a[x, y, z], numeric(1))

  # Convert t -> partial r -> Fisher z
  # Guard against NA/Inf
  if (!all(is.finite(t_vals)) || !all(is.finite(df_vec))) next
  r_vals <- mapply(t2r, t = t_vals, df = df_vec)
  # avoid r hitting +/-1
  r_vals[r_vals >=  0.999999] <-  0.999999
  r_vals[r_vals <= -0.999999] <- -0.999999
  z_vals <- r2z(r_vals)

  # Build data frame for metafor (random effects with study/session nesting)
  meta_df <- data.frame(
    study     = db$study[seq_along(z_vals)],
    session   = db$session[seq_along(z_vals)],
    yi        = z_vals,
    vi        = se_vals^2
  )

  # Remove any NA variance rows
  meta_df <- meta_df[is.finite(meta_df$yi) & is.finite(meta_df$vi) & meta_df$vi > 0, ]
  if (nrow(meta_df) < 2) next

  # Fit random-effects model
  fit <- try(
    metafor::rma.mv(yi = yi, V = vi,
                    data = meta_df,
                    random = list(~1 | study/session),
                    method = "REML"),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) next

  meta_z[x, y, z]    <- fit$zval
  meta_p[x, y, z]    <- fit$pval
  meta_beta[x, y, z] <- as.numeric(fit$beta)  # intercept in z-units
  if (ii %% 10000 == 0) cat("Processed", ii, "of", nrow(mask_idx), "voxels\n")
}

# ---- Post-processing / FDR ----
# Clean extreme z (occasional numerical instability)
meta_z[meta_z >  50] <- NA
meta_z[meta_z < -50] <- NA

# Apply mask
meta_z[mask_data < 1]    <- NA
meta_p[mask_data < 1]    <- NA
meta_beta[mask_data < 1] <- NA

# FDR across masked voxels
p_vec <- meta_p[mask_data > 0]
p_fdr_vec <- rep(NA_real_, length(p_vec))
if (length(p_vec)) {
  p_fdr_vec <- p.adjust(p_vec, method = "fdr")
}
meta_p_fdr <- array(NA_real_, dim = dims)
meta_p_fdr[mask_data > 0] <- p_fdr_vec

# Thresholded FDR map (q < 0.05)
meta_p_fdr_thr <- meta_p_fdr
meta_p_fdr_thr[meta_p_fdr_thr >= 0.05 | !is.finite(meta_p_fdr_thr)] <- NA

# ---- Write outputs ----
to_nifti_like <- function(arr, template_hdr) {
  nifti(arr,
        dim     = dim(template_hdr),
        datatype= template_hdr@datatype,
        pixdim  = template_hdr@pixdim)
}

z_nii     <- to_nifti_like(meta_z,    ref_hdr)
p_nii     <- to_nifti_like(meta_p,    ref_hdr)
pFDR_nii  <- to_nifti_like(meta_p_fdr,ref_hdr)
pFDRt_nii <- to_nifti_like(meta_p_fdr_thr, ref_hdr)
beta_nii  <- to_nifti_like(meta_beta, ref_hdr)

writeNIfTI(z_nii,     file.path(output_dir, "meta_z_NCvsBD_PCA"))
writeNIfTI(p_nii,     file.path(output_dir, "meta_p_NCvsBD_PCA"))
writeNIfTI(pFDR_nii,  file.path(output_dir, "meta_pFDR_NCvsBD_PCA"))
writeNIfTI(pFDRt_nii, file.path(output_dir, "meta_pFDR_thr_NCvsBD_PCA"))
writeNIfTI(beta_nii,  file.path(output_dir, "meta_beta_NCvsBD_PCA"))

cat("Saved NIfTI outputs to:", output_dir, "\n")

# ---- Peak extraction (optional) ----
# If you want cluster peaks from a thresholded z/p map, you can use neurobase::cluster.
# Example below assumes you have a t-map called 'tmap.nii.gz' and threshold = 2.3:
# library(neurobase)
# tmap <- readNIfTI(file.path(output_dir, "tmap.nii.gz"), reorient = FALSE)
# thr  <- 2.3
# tmap_bin <- array(0, dim = dim(tmap)); tmap_bin[tmap >= thr] <- 1
# cl <- neurobase::cluster(tmap_bin, k = 26)
# lbl <- cl$img
# idx <- which(lbl > 0, arr.ind = TRUE)
# df  <- as.data.frame(idx)
# df$cluster <- lbl[lbl > 0]
# df$tval    <- tmap[lbl > 0]
# library(dplyr)
# cluster_summary <- df %>%
#   group_by(cluster) %>%
#   summarise(
#     size   = n(),
#     max_t  = max(tval),
#     mean_t = mean(tval),
#     peak_x = x[which.max(tval)],
#     peak_y = y[which.max(tval)],
#     peak_z = z[which.max(tval)]
#   ) %>% arrange(desc(max_t))
# write.csv(cluster_summary, file.path(output_dir, "cluster_peaks.csv"), row.names = FALSE)
