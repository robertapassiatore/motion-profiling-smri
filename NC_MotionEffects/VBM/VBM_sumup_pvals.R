library(RNifti)

# --- Define full paths to the p-maps ---
p_paths <- c(
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/UNIBA2/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/HCP/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/DNS/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/LIBD/FI/pval_0001.nii",
  "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/results/UNIBA2/FI/pval_0001.nii"
)

# --- Path to binary mask ---
mask_path <- "/Users/robertapassiatore/Downloads/rrGM_cat12_bin.nii"
r2_path = "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/mean_r2_0001_FI.nii"

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
writeNifti(combined_nifti, "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/fisher_pval_0001_FI.nii.gz")



#### crop R2 map
# --- Load data ---
mask <- readNifti(mask_path)
mask_idx <- which(mask > 0)

r2_map_full <- readNifti(r2_path)
p_maps <- lapply(p_paths, readNifti)

# --- Extract p-values within mask and apply Fisher's method ---
p_vals_list <- lapply(p_maps, function(p) {
  vals <- p[mask_idx]
  vals[vals <= 0] <- .Machine$double.eps
  return(vals)
})

k <- length(p_vals_list) - 1
fisher_stat <- Reduce(`+`, lapply(p_vals_list, function(p) -2 * log(p)))
fisher_pvals <- pchisq(fisher_stat, df = 2 * k, lower.tail = FALSE)

# --- FDR correction ---
fisher_pvals_fdr <- p.adjust(fisher_pvals, method = "fdr")

# --- Apply FDR p < 0.05 threshold to R² map ---
r2_thresh <- array(NA_real_, dim = dim(r2_map_full))
r2_thresh[mask_idx[fisher_pvals_fdr < 0.05]] <- r2_map_full[mask_idx[fisher_pvals_fdr < 0.05]]

# --- Apply p < 0.001 threshold to R² map ---
r2_thresh <- array(NA_real_, dim = dim(r2_map_full))
r2_thresh[mask_idx[fisher_pvals < 0.001]] <- r2_map_full[mask_idx[fisher_pvals < 0.001]]

# --- Save result ---
r2_thresh_nifti <- asNifti(r2_thresh, reference = r2_map_full)
writeNifti(r2_thresh_nifti, "/Volumes/HD2/PROJECTS/RP_MotionProfiling/vbm_nc/mean_r2_0001_FI_uncorr_thresh_0.001.nii.gz")



