
# Motion Profiling – R² Estimation and P-value Combination

## Code Overview

This folder contains scripts used to estimate the **variance explained (R²)** by the **five motion principal components (PCs)** on **grey matter volumes**.  
Separate code is provided for **ROI-based** and **VBM (voxel-based morphometry)** analyses.

---

## VBM Analysis

### R² Estimation  
To obtain voxel-wise R² values, contrast images (*t*-maps and *F*-maps) generated from **SPM12** are used as input for:

- `VBM_extract_r2_tmaps.m`  
- `VBM_extract_r2_fmaps.m`

### P-value Estimation  
Similarly, to obtain voxel-wise *p*-values from the same contrast images:

- `VBM_extract_pval_tmaps.m`  
- `VBM_extract_pval_fmaps.m`

### Combining Results Across Datasets  
To combine *p*-values voxel-wise across multiple datasets, the script implementing **Fisher’s method** is provided:

- `VBM_sumup_pvals.R`

**Note:** Example input NIfTI files from the HCP and CNP datasets are included. NIfTI files of the effects shown in *Figure 3* are included in this folder.  
The figure was generated using **MRIcroGL**.

---

## ROI Analysis

### R² Estimation  
To obtain ROI-based R² values, site-specific matrices (subject × ROI + PCs + covariates) are used as input for:

- `ROI_PCEffect_permutations_EffectSize_allPC.R`

### Combining Results Across Datasets  
To combine ROI-based *p*-values across multiple datasets, the script implementing **Fisher’s method** is provided:

- `ROI_LinearEffectPCs_metap.R`

**Note:** Example CSV files from the HCP and CNP datasets are included. CSV files of the effects shown in *Figure 4* are included in this folder.  
The figure was generated using **ggseg**.


---

### Author  
**Roberta Passiatore**  
*Updated: 2025-10-01*  
