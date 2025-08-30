# Install necessary packages if you don't have them
install.packages("neuroim")
install.packages("oro.nifti")
install.packages("metafor")
install.packages("dplyr")
install.packages("AnalyzeFMRI")

# Load the packages
#library(neuroim)
library(oro.nifti)
library(metafor)
library(dplyr)
#library(AnalyzeFMRI)  # optional



## Define helper functions -----------------------------------------

# Convert t-statistic to partial r
t2r = function(t, df) {
  sign(t) * sqrt(t^2 / (t^2 + df))
}

# Fisher’s r-to-z transformation
r2z = function(r) {
  0.5 * log((1 + r) / (1 - r))
}

# Compute standard error from r-squared and n
r2se = function(r, n) {
  sqrt((1 - r^2) / n)
}


input_dir <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/tmaps"
setwd(input_dir)
tmaps <- list.files(pattern=glob2rx("r*_PCA*.nii.gz"))

# Load t-map files
nifti_file1 <- readNIfTI(tmaps[1], reorient = T)
nifti_file1@.Data[is.na(nifti_file1@.Data)] <- 0
tmap_data1 <- nifti_file1@.Data
nifti_file2 <- readNIfTI(tmaps[2], reorient = T)
nifti_file2@.Data[is.na(nifti_file2@.Data)] <- 0
tmap_data2 <- nifti_file2@.Data
nifti_file3 <- readNIfTI(tmaps[3], reorient = T)
nifti_file3@.Data[is.na(nifti_file3@.Data)] <- 0
tmap_data3 <- nifti_file3@.Data
nifti_file4 <- readNIfTI(tmaps[4], reorient = T)
nifti_file4@.Data[is.na(nifti_file4@.Data)] <- 0
tmap_data4 <- nifti_file4@.Data
nifti_file5 <- readNIfTI(tmaps[5], reorient = T)
nifti_file5@.Data[is.na(nifti_file5@.Data)] <- 0
tmap_data5 <- nifti_file5@.Data
nifti_file6 <- readNIfTI(tmaps[6], reorient = T)
nifti_file6@.Data[is.na(nifti_file6@.Data)] <- 0
tmap_data6 <- nifti_file6@.Data
nifti_file7 <- readNIfTI(tmaps[7], reorient = T)
nifti_file7@.Data[is.na(nifti_file7@.Data)] <- 0
tmap_data7 <- nifti_file7@.Data
nifti_file8 <- readNIfTI(tmaps[8], reorient = T)
nifti_file8@.Data[is.na(nifti_file8@.Data)] <- 0
tmap_data8 <- nifti_file8@.Data
nifti_file9 <- readNIfTI(tmaps[9], reorient = T)
nifti_file9@.Data[is.na(nifti_file9@.Data)] <- 0
tmap_data9 <- nifti_file9@.Data
nifti_file10 <- readNIfTI(tmaps[10], reorient = T)
nifti_file10@.Data[is.na(nifti_file10@.Data)] <- 0
tmap_data10 <- nifti_file10@.Data
nifti_file11 <- readNIfTI(tmaps[11], reorient = T)
nifti_file11@.Data[is.na(nifti_file11@.Data)] <- 0
tmap_data11 <- nifti_file11@.Data
nifti_file12 <- readNIfTI(tmaps[12], reorient = T)
nifti_file12@.Data[is.na(nifti_file12@.Data)] <- 0
tmap_data12 <- nifti_file12@.Data

#nifti_file4 <- readNIfTI("/Users/roberta.passiatore/Downloads/metanalisis/wo_removed_data_LIBD/spmT_0001.nii", reorient = T) 
input_dir <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD/se"
setwd(input_dir)
se <- list.files(pattern=glob2rx("r*_noPCA*.nii.gz"))

# Load se files 
nifti_file1 <- readNIfTI(se[1], reorient = T)
nifti_file1@.Data[is.na(nifti_file1@.Data)] <- 0
se_data1 <- nifti_file1@.Data
nifti_file2 <- readNIfTI(se[2], reorient = T)
nifti_file2@.Data[is.na(nifti_file2@.Data)] <- 0
se_data2 <- nifti_file2@.Data
nifti_file3 <- readNIfTI(se[3], reorient = T)
nifti_file3@.Data[is.na(nifti_file3@.Data)] <- 0
se_data3 <- nifti_file3@.Data
nifti_file4 <- readNIfTI(se[4], reorient = T)
nifti_file4@.Data[is.na(nifti_file4@.Data)] <- 0
se_data4 <- nifti_file4@.Data
nifti_file5 <- readNIfTI(se[5], reorient = T)
nifti_file5@.Data[is.na(nifti_file5@.Data)] <- 0
se_data5 <- nifti_file5@.Data
nifti_file6 <- readNIfTI(se[6], reorient = T)
nifti_file6@.Data[is.na(nifti_file6@.Data)] <- 0
se_data6 <- nifti_file6@.Data
nifti_file7 <- readNIfTI(se[7], reorient = T)
nifti_file7@.Data[is.na(nifti_file7@.Data)] <- 0
se_data7 <- nifti_file7@.Data
nifti_file8 <- readNIfTI(se[8], reorient = T)
nifti_file8@.Data[is.na(nifti_file8@.Data)] <- 0
se_data8 <- nifti_file8@.Data
nifti_file9 <- readNIfTI(se[9], reorient = T)
nifti_file9@.Data[is.na(nifti_file9@.Data)] <- 0
se_data9 <- nifti_file9@.Data
nifti_file10 <- readNIfTI(se[10], reorient = T)
nifti_file10@.Data[is.na(nifti_file10@.Data)] <- 0
se_data10 <- nifti_file10@.Data
nifti_file11 <- readNIfTI(se[11], reorient = T)
nifti_file11@.Data[is.na(nifti_file11@.Data)] <- 0
se_data11 <- nifti_file11@.Data
nifti_file12 <- readNIfTI(se[12], reorient = T)
nifti_file12@.Data[is.na(nifti_file12@.Data)] <- 0
se_data12 <- nifti_file12@.Data


# Load a mask file that defines voxels of interest
mask_file <- readNIfTI("/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/rGM_cat12_bin.nii.gz", reorient = TRUE)
mask_data <- mask_file@.Data

# Determine sample sizes for each study
db=read.csv('/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/res_bd.csv')


# Get dimensions from the mask (assumed same as t-maps)
dims <- dim(mask_data)
dim_x <- dims[1]
dim_y <- dims[2]
dim_z <- dims[3]

# Create empty arrays to store meta-analysis results for each voxel
meta_z    <- array(NA, dim = dims)  # meta-analysis z-values
meta_p    <- array(NA, dim = dims)  # meta-analysis p-values
meta_beta <- array(NA, dim = dims)  # meta-analysis effect sizes

# Loop over each voxel in the mask --------------------------------
cat("Starting voxel-wise meta-analysis...\n")
for (x in 1:dim_x) {
  for (y in 1:dim_y) {
    for (z in 1:dim_z) {
      if (mask_data[x, y, z] > 0) {  # Process only voxels within the mask
        
        # Extract t-values from each study for the current voxel
        t1 <- tmap_data1[x, y, z]
        t2 <- tmap_data2[x, y, z]
        t3 <- tmap_data3[x, y, z]
        t4 <- tmap_data4[x, y, z]
        t5 <- tmap_data5[x, y, z]
        t6 <- tmap_data6[x, y, z]
        #t7 <- tmap_data7[x, y, z]
        #t8 <- tmap_data8[x, y, z]
        #t9 <- tmap_data9[x, y, z]
        #t10 <- tmap_data10[x, y, z]
        #t11 <- tmap_data11[x, y, z]
        t#12 <- tmap_data12[x, y, z]

        
        # Degrees of freedom for each study
        df1 <- db$N[1] - 11 # 5 in the model noPCA and 10 in the model PCA + 1 if multiple sites
        df2 <- db$N[2] - 11
        df3 <- db$N[3] - 11
        df4 <- db$N[4] - 10
        df5 <- db$N[5] - 10
        df6 <- db$N[6] - 10
        #df7 <- db$N[7] - 5
        #df8 <- db$N[8] - 5
        #df9 <- db$N[9] - 5
        #df10 <- db$N[10] - 5
        #df11 <- db$N[11] - 5
        #df12 <- db$N[12] - 5
        
   
        
        # Convert t-values to effect sizes (using r2z(t2r(t, df)))
        es1 <- r2z(t2r(t1, df1))
        es2 <- r2z(t2r(t2, df2))
        es3 <- r2z(t2r(t3, df3))
        es4 <- r2z(t2r(t4, df4))
        es5 <- r2z(t2r(t5, df5))
        es6 <- r2z(t2r(t6, df6))
        #es7 <- r2z(t2r(t7, df7))
        #es8 <- r2z(t2r(t8, df8))
        #es9 <- r2z(t2r(t9, df9))
        #es10 <- r2z(t2r(t10, df10))
        #es11 <- r2z(t2r(t11, df11))
        #es12 <- r2z(t2r(t12, df12))

        
        
        # Compute standard errors from effect sizes
        se1 <- se_data1[x, y, z]
        se2 <- se_data2[x, y, z]
        se3 <- se_data3[x, y, z]
        se4 <- se_data4[x, y, z]
        se5 <- se_data5[x, y, z]
        se6 <- se_data6[x, y, z]
        #se7 <- se_data7[x, y, z]
        #se8 <- se_data8[x, y, z]
        #se9 <- se_data9[x, y, z]
        #se10 <- se_data10[x, y, z]
        #se11 <- se_data11[x, y, z]
        #se12 <- se_data12[x, y, z]
 
        
        
        # Build a data frame for meta-analysis
        meta_data <- data.frame(
          study     = db$study,
          session   = db$session,
          z_value   = c(es1, es2, es3, es4, es5, es6),
                        #es7, es8, es9, es10, es11, es12),
          std.error = c(se1, se2, se3, se4, se5, se6)
                        #se7, se8, se9, se10, se11, se12)
        )
        
        
        # Perform meta-analysis using metafor package
        res <- try(
          rma.mv(yi = z_value, V = std.error^2,
                 data = meta_data,
                 random = list(~1 | study/session),
                 method = "REML"),
          silent = TRUE
        )
        
        # Only assign results if model ran successfully
        if (!inherits(res, "try-error")) {
          meta_z[x, y, z]     <- res$zval
          meta_p[x, y, z]     <- res$pval
          meta_beta[x, y, z]  <- res$beta
        }
      }  # end if mask
    }
  }
  cat("Processed slice", x, "of", dim_x, "\n")
}


# Write out the results as NIfTI images -------------------------
# Use the mask header as a template to retain spatial information
# Create array to store FDR-corrected p-values
meta_z[meta_z > 50] = NA
meta_z[meta_z < -50] = NA
meta_z[mask_data < 1] = NA
meta_beta[mask_data < 1] = NA
meta_p[mask_data < 1] = NA
meta_p_fdr = array(NA, dim = dims)

# Step 1: Get logical mask of voxels to include
mask_idx <- which(mask_data > 0, arr.ind = TRUE)

# Step 2: Extract p-values for those voxels
p_vals <- apply(mask_idx, 1, function(idx) {
  x <- idx[1]; y <- idx[2]; z <- idx[3]
  meta_p[x, y, z]
})

# Step 3: Apply FDR correction
p_vals_fdr <- p.adjust(p_vals, method = "fdr")

# Step 4: Write corrected values back to meta_p_fdr
for (i in 1:nrow(mask_idx)) {
  x <- mask_idx[i, 1]
  y <- mask_idx[i, 2]
  z <- mask_idx[i, 3]
  meta_p_fdr[x, y, z] <- p_vals_fdr[i]
}
meta_p_fdr[mask_data < 1] = NA


z_nii <- nifti(meta_z, 
               dim = dims, 
               datatype = mask_file@datatype, 
               pixdim = mask_file@pixdim)
orthographic(z_nii)

p_nii <- nifti(meta_p, 
               dim = dims, 
               datatype = mask_file@datatype, 
               pixdim = mask_file@pixdim)
orthographic(p_nii)

pFDR_nii <- nifti(meta_p_fdr, 
               dim = dims, 
               datatype = mask_file@datatype, 
               pixdim = mask_file@pixdim)
orthographic(pFDR_nii)

# Create a copy to preserve original values
meta_p_fdr_thr <- meta_p_fdr

# Apply threshold: keep only values < 0.05
meta_p_fdr_thr[meta_p_fdr_thr < 0 | is.na(meta_p_fdr_thr)] <- NA
meta_p_fdr_thr[meta_p_fdr_thr > 0.05 | is.na(meta_p_fdr_thr)] <- NA


pFDR_thr_nii <- nifti(meta_p_fdr_thr, 
                  dim = dims, 
                  datatype = mask_file@datatype, 
                  pixdim = mask_file@pixdim)
orthographic(pFDR_thr_nii)

beta_nii <- nifti(meta_beta, 
                  dim = dims, 
                  datatype = mask_file@datatype, 
                  pixdim = mask_file@pixdim)

input_dir <- "/Volumes/HD2/PROJECTS/RP_MotionProfiling/tmaps_GroupEffect/BD"
setwd(input_dir)
writeNIfTI(z_nii, filename = "meta_z_NCvsBD_PCA")
writeNIfTI(p_nii, filename = "meta_p_NCvsBD_PCA")
writeNIfTI(pFDR_nii, filename = "meta_pFDR_NCvsBD_PCA")
writeNIfTI(pFDR_thr_nii, filename = "meta_pFDR_thr_NCvsBD_PCA")
writeNIfTI(beta_nii, filename = "meta_beta_NCvsBD_PCA")

# library(neurobase)
# 
# # Load background (e.g., MNI152) and overlay
# bg <- readNIfTI("MNI152_T1_2mm_brain.nii.gz", reorient = FALSE)
# overlay <- readNIfTI("meta_p_fdr_thresh_lt0.05.nii.gz", reorient = FALSE)
# 
# # Plot overlay
# ortho2(bg, overlay, col.y = hotmetal())

#### extract peak info

tmap <- readNIfTI("tmap.nii.gz", reorient = FALSE)
threshold <- 2.3
tmap_thresh <- tmap
tmap_thresh[tmap < threshold] <- 0
tmap_thresh[tmap >= threshold] <- 1

# Using neurobase (simple and fast)
cluster_map <- cluster(tmap_thresh, connectivity = 26)

# Extract labeled clusters (each voxel has a cluster ID > 0)
cluster_labels <- cluster_map$img

# Get voxel-wise indices and t-values
voxel_info <- which(cluster_labels > 0, arr.ind = TRUE)
voxel_info <- as.data.frame(voxel_info)
voxel_info$cluster <- cluster_labels[cluster_labels > 0]
voxel_info$tval <- tmap[cluster_labels > 0]

# Compute cluster-level summaries
cluster_summary <- voxel_info %>%
  group_by(cluster) %>%
  summarise(
    size = n(),
    max_t = max(tval),
    mean_t = mean(tval),
    peak_x = x[which.max(tval)],
    peak_y = y[which.max(tval)],
    peak_z = z[which.max(tval)]
  ) %>%
  arrange(desc(max_t))






