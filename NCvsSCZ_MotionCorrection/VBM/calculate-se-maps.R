#F = (R²/(p-1)) / ((1-R²)/(n-p))
#R² = F / (F + (n-p)/(p-1))
#Adjusted R² = 1 - [(1-R²) * (n-1) / (n-p-1)]

library(RNifti)
library(ANTsR)
library(pixmap)

# === CONFIGURATION ===
input_dir <- "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/con"
output_prefix <- "r"
reference_path <- "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/rGM_cat12_bin.nii.gz"

# Load reference image
reference <- antsImageRead(reference_path, dimension = 3)

# List mask files in input_dir
mask_files <- list.files(input_dir, pattern = "nii.gz", full.names = TRUE)

# Loop over each mask file
for (mask_path in mask_files) {
  # Load mask
  mask <- antsImageRead(mask_path, dimension = 3)
  
  # Resample to match reference
  mask_resliced <- resampleImageToTarget(mask, reference, interpType = 1)
  
  # Construct output path
  file_name <- basename(mask_path)
  output_path <- file.path(input_dir, paste0(output_prefix, file_name))
  
  # Save resliced mask
  antsImageWrite(mask_resliced, output_path)
  plot(mask_resliced, axis = 3, nslices = 12)
  
  cat("[✓] Resliced mask saved to", output_path, "\n")
}


#plot(reference, axis = 3, nslices = 9)

library(RNifti)

# Load NIfTI files
i1 <- readNifti("/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/BD/tmaps/rUNIBA2_FI_noPCA_spmT_0004_maskGM.nii.gz")  # t-map
i2 <- readNifti("/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/BD/con/UNIBA2_FI_noPCA_con.nii.gz")  # contrast (beta)
mask <- readNifti("rGM_cat12_bin.nii.gz")

# Ensure the mask is binary
mask_bin <- mask > 0

# Check all images have the same dimensions
stopifnot(all(dim(i1) == dim(i2), dim(i1) == dim(mask_bin)))

# Apply mask
i1_masked <- i1 * mask_bin
i2_masked <- i2 * mask_bin

# Avoid division by zero
epsilon <- 1e-6

# Compute standard error (SE = beta / t)
se_map <- i2_masked / (i1_masked + epsilon)

# Define conditions
out_of_mask <- !mask_bin
too_small_t <- abs(i1_masked) < epsilon
too_large_se <- se_map >= 1
too_small_se <- se_map <= -1

# Apply combined masking rule
se_map[out_of_mask | too_small_t | too_large_se | too_small_se] <- epsilon

# Invert sign if SE is negative (make all values positive)
se_map[se_map < 0] <- -se_map[se_map < 0]

# Save result
writeNifti(se_map, "UNIBA2_FI_noPCA_SE.nii.gz")

#cat("[✓] SE map saved as LIBD_WM_PCA_SE.nii.gz\n")

input_dir <- "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/tmaps"
setwd(input_dir)
tmaps <- list.files(pattern=glob2rx("r*.nii.gz"))
input_dir <- "/Users/roberta.passiatore/Downloads/tmaps_GroupEffect/SCZ/se"
setwd(input_dir)
se <- list.files(pattern=glob2rx("r*.nii.gz"))


# Build a data frame for meta-analysis
data.test <- data.frame(
  study    = data.sub$Site,
  z = data.sub$zscore,
  session = data.sub$Session,
  std.err  = data.sub$SE
)

# Perform meta-analysis using the metafor package
#res <- rma(yi = cohens_d, sei = std.err, data = data.test, method = "ML")
# Two‐level: random intercept for each study
res <- rma.mv(yi = z, V = std.err^2, 
              data = data.test,
              random   = list(~ 1|study/session),
              method = "ML")
summary(res)

#### p value maps

# Load NIfTI files
z_img <- readNifti("/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/szmap_NCvsBD_noPCA.nii")  # t-map
mask_img <- readNifti("/Users/robertapassiatore/Downloads/rGM_cat12_bin.nii")

# --- Convert z-scores to p-values within the mask ---
p_map <- array(NA_real_, dim = dim(z_img))

# Two-tailed p-value calculation
p_map[mask_img > 0] <- 2 * (1 - pnorm(abs(z_img[mask_img > 0])))

# --- Create a new NIfTI image for the p-value map ---
p_nifti <- asNifti(p_map, reference = z_img)

# --- Save the result ---
writeNifti(p_nifti, file = "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/meta_p_NCvsBD_noPCA.nii.gz")



library(RNifti)

# Load images
z_img <- readNifti("/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/szmap_NCvsBD_PCA.nii")  # t-map
diff_img <- readNifti("/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/diff_zmap_NCvsBD.nii")  # t-map

mask_img <- readNifti("/Users/robertapassiatore/Downloads/rGM_cat12_bin.nii")

# Get voxel indices within the mask
mask_indices <- which(mask_img > 0)

# Compute uncorrected two-tailed p-values
p_values <- 2 * (1 - pnorm(abs(z_img[mask_indices])))

# Apply FDR correction (or use "bonferroni")
p_corrected <- p.adjust(p_values, method = "fdr")

# Create full p-value map and fill in corrected p-values
p_map <- array(NA_real_, dim = dim(z_img))
p_map[mask_indices] <- p_corrected


# Save as NIfTI
p_nifti <- asNifti(p_map, reference = z_img)
writeNifti(p_nifti, "meta_pFDR_NCvsSCZ_PCA.nii.gz")



# Step 1: Compute uncorrected two-tailed p-values
p_uncorrected <- 2 * (1 - pnorm(abs(z_img[mask_indices])))

# Step 2: FDR correction
p_corrected <- p.adjust(p_uncorrected, method = "fdr")

# Step 3: Create boolean mask for significant voxels (p < 0.05)
sig_voxels <- logical(length(mask_indices))
sig_voxels[p_corrected < 0.001] <- TRUE

# Step 4: Create output Z-map (default to NA)
z_thresh <- array(NA_real_, dim = dim(z_img))
z_thresh[mask_indices[sig_voxels]] <- z_img[mask_indices[sig_voxels]]

# Step 5: Save the thresholded Z-map
z_thresh_nifti <- asNifti(z_thresh, reference = z_img)
writeNifti(z_thresh_nifti, "zmap_FDR_thresh_0.05.nii.gz")




# Step 1: Compute uncorrected two-tailed p-values
p_uncorrected <- 2 * (1 - pnorm(abs(z_img[mask_indices])))

# Step 2: FDR correction
p_corrected <- p.adjust(p_uncorrected, method = "bonferroni")

# Step 3: Create boolean mask for significant voxels (p < 0.05)
sig_voxels <- logical(length(mask_indices))
sig_voxels[p_corrected < 0.05] <- TRUE

# Step 4: Create output Z-map (default to NA)
z_thresh <- array(NA_real_, dim = dim(z_img))
z_thresh[mask_indices[sig_voxels]] <- z_img[mask_indices[sig_voxels]]

# Step 5: Save the thresholded Z-map
z_thresh_nifti <- asNifti(z_thresh, reference = z_img)
writeNifti(z_thresh_nifti, "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/zmap_FDR_thresh_0.05_PCA.nii.gz")

