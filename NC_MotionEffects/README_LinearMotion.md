
# Motion Profiling – R² Estimation and P-value Combination

## Code Overview

This folder contains scripts used to estimate the **variance explained (R²)** by the **five motion principal components (PCs)** on **grey matter volumes**.  
Separate code is provided for **ROI-based** analyses.

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
*Updated: 2026-03-05*  
