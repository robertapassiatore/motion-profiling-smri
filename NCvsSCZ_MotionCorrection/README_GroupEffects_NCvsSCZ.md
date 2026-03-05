
# Group Effects Analysis – NC vs SCZ

This folder contains scripts used to compute **group effects analysis** on grey matter volumes between **NC (controls)** and **SCZ (patients)**.  
Separate code is provided for **ROI-based** analyses.

---

## ROI Analysis

### R² Estimation  
To obtain ROI-based R² and z values, site-specific matrices (subject × ROI + PCs + covariates) are used as input for:

- `ROIs_GroupEffect_permutations_EffectSize_SCZ.R`

### Mixed-Effects Meta-Analysis  
To obtain meta-analytic statistics across datasets and sessions, a **mixed-effects meta-analysis** is performed, treating *session* and *dataset* as random effects:

- `ROIs_GroupEffect_metanalysis_SCZ.R`


**Note:** CVS files reporting regression results for all cohorts are provided. CSV files of the effects shown in *Figure 6* are included in this folder.  
The figure was generated using **ggseg** and code `plot_GroupDiff_NCvsSCZvsBD.R`

---

### Author  
**Roberta Passiatore**  
*Updated: 2025-03-05*  
