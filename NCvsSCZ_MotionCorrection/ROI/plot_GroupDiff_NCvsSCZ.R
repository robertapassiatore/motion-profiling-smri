# Install required packages if missing
if (!requireNamespace("ggseg", quietly = TRUE)) install.packages("ggseg")
if (!requireNamespace("ggsegSchaefer", quietly = TRUE)) install.packages("ggsegSchaefer")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("grid", quietly = TRUE)) install.packages("grid")

# Load libraries
library(ggseg)
library(ggsegSchaefer)
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)


results_file = 'res.GroupDiff_NCvsBD_corticalROI_REML_standard_metanalysis.csv'
roi_file = "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/roi-names.csv"
atlas_name = "schaefer17_100"
output_pdf = 'Metanalysis-results_corticalROI_NCvsBD_ME_pFDR.pdf'
analysis = 'standard'

# Read ROI names and regression results with stringsAsFactors = FALSE
rois <- read.csv(file = roi_file, stringsAsFactors = FALSE)
results <- read.csv(file = results_file, stringsAsFactors = FALSE)

# Create y-axis annotation and wrap as a patchwork element
# y_label <- textGrob("    FI                   WM                    RS",
#                     rot = 90,  # Rotate text vertically
#                     gp = gpar(fontsize = 14, fontface = "bold"))
# y_annotation <- wrap_elements(full = y_label)

y_label <- textGrob("         Corrected        Standard",
                    rot = 90,  # Rotate text vertically
                    gp = gpar(fontsize = 14, fontface = "bold"))
y_annotation <- wrap_elements(full = y_label)


# Add atlas labels to results by matching on ROI names
results$label <- rois$Atlas[match(results$ROI, rois$ROI)]

# Load the specified atlas from ggsegSchaefer
data(list = atlas_name, package = "ggsegSchaefer")
atlas_obj <- get(atlas_name)

# Define analysis types to loop over
analysis_types <- unique(results$session)

# Create an empty list to store plots
plot_list <- list()

#for(analysis in analysis_types) {
  # Filter results by analysis type
  #results_subset <- results[results$session == analysis, ]
  #results_subset$pFDR = p.adjust(results_subset$pval, method = 'fdr')
  results$pFDR = p.adjust(results$pval, method = 'fdr')
  
  
  # Subset to relevant columns: label, pFDR, and Cohens_d
  #pvalue_df <- results_subset[, c("label", "pFDR", "zval")]
  pvalue_df <- results[, c("label", "pFDR", "zval")]
  
  
  # Create a copy of the atlas object for merging
  atlas_temp <- atlas_obj
  atlas_temp$data <- left_join(atlas_temp$data, pvalue_df, by = "label", 
                               relationship = "many-to-many")
  
  # Set pFDR > 0.05 to NA and mask Cohen's d accordingly
  atlas_temp$data$pFDR[atlas_temp$data$pFDR > 0.03] <- NA
  atlas_temp$data$zval[is.na(atlas_temp$data$pFDR)] <- NA
  
  # Create pFDR plot
  p_pFDR <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = pFDR),
                  position = "dispersed") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = "grey", limits = c(0, 0.05)) +
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                       axis.text.x=element_text(size=10,  family="sans"))
  p_pFDR
  # Optionally add ggtitle() if desired
  
  # Create Cohen's d plot
  p_zval <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = zval),
                  position = "dispersed") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey" , limits = c(-15, 15)) + #
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                          axis.text.x=element_text(size=10,  family="sans"))
  p_zval 
  # Optionally add ggtitle() if desired
  
  # Store the plots in the list with keys "standard_pFDR", etc.
  plot_list[[paste0(analysis, "_pval")]] <- p_pFDR
  plot_list[[paste0(analysis, "_zval")]] <- p_zval
#}



# Combine the plots using patchwork:
# For pFDR plots: y_annotation on the left and standard above corrected.
combined_pFDR <- y_annotation + (plot_list$RS_pval / plot_list$WM_pval / plot_list$FI_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_pFDR
# For Cohen's d plots: same layout.
combined_zval <- y_annotation + (plot_list$RS_zval / plot_list$WM_zval / plot_list$FI_zval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_zval

# Combine the plots using patchwork:
# For pFDR plots: y_annotation on the left and standard above corrected.
combined_pFDR <- y_annotation + (plot_list$standard_pval / plot_list$corrected_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_pFDR
# For Cohen's d plots: same layout.
combined_zval <- y_annotation + (plot_list$standard_zval / plot_list$corrected_zval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_zval

# Open a PDF device and print both combined plots (each on a separate page)
pdf(file = output_pdf, width = 8, height = 6)
print(combined_pFDR)
print(combined_zval)
dev.off()

#Plot Delta
results_file = 'res.GroupDiff_NCvsBD_corticalROI_REML_delta_metanalysis.csv'
output_pdf = 'Metanalysis-results_corticalROI_NCvsBD_ME_Delta_pFDR.pdf'

# Read ROI names and regression results with stringsAsFactors = FALSE
rois <- read.csv(file = roi_file, stringsAsFactors = FALSE)
results <- read.csv(file = results_file, stringsAsFactors = FALSE)

# Create y-axis annotation and wrap as a patchwork element
y_label <- textGrob("    Delta",
                    rot = 90,  # Rotate text vertically
                    gp = gpar(fontsize = 14, fontface = "bold"))
y_annotation <- wrap_elements(full = y_label)

# Add atlas labels to results by matching on ROI names
results$label <- rois$Atlas[match(results$ROI, rois$ROI)]

# Load the specified atlas from ggsegSchaefer
data(list = atlas_name, package = "ggsegSchaefer")
atlas_obj <- get(atlas_name)

# Define analysis types to loop over
analysis_types <- unique(results$session)

# Create an empty list to store plots
plot_list <- list()

for(analysis in analysis_types) {
  # Filter results by analysis type
  results_subset <- results#[results$session == analysis, ]
  results_subset$QEpFDR = p.adjust(results_subset$QEp, method = 'fdr')
  results_subset$pFDR = p.adjust(results_subset$pval, method = 'fdr')
  
  # Subset to relevant columns: label, pFDR, and Cohens_d
  pvalue_df <- results_subset[, c("label", "QEpFDR",  "QE", "zval", 'pFDR')]
  
  # Create a copy of the atlas object for merging
  atlas_temp <- atlas_obj
  atlas_temp$data <- left_join(atlas_temp$data, pvalue_df, by = "label", 
                               relationship = "many-to-many")
  
  # Set pFDR > 0.05 to NA and mask Cohen's d accordingly
  atlas_temp$data$pFDR[atlas_temp$data$pFDR > 0.05] <- NA
  atlas_temp$data$QE[is.na(atlas_temp$data$pFDR)] <- NA
  atlas_temp$data$QEpFDR[is.na(atlas_temp$data$pFDR)] <- NA
  atlas_temp$data$zval[is.na(atlas_temp$data$pFDR)] <- NA
  
  # Create pFDR plot
  p_pFDR <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = QE),
                  position = "dispersed") +
    scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey", limits = c(0, 20)) +
    theme_minimal()
  # Optionally add ggtitle() if desired
  
  # Create Cohen's d plot
  p_delta <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = zval),
                  position = "dispersed") +
    scale_fill_gradient(low='white', high = "darkgreen", na.value = "grey" , limits = c(0, 6)) + #
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                          axis.text.x=element_text(size=10,  family="sans")) +
    labs(fill=expression(delta))
  # Optionally add ggtitle() if desired
  
  # Store the plots in the list with keys "standard_pFDR", etc.
  plot_list[[paste0(analysis, "_pval")]] <- p_pFDR
  plot_list[[paste0(analysis, "_delta")]] <- p_delta
}



# Combine the plots using patchwork:
# For pFDR plots: y_annotation on the left and standard above corrected.
combined_pFDR <- y_annotation + (plot_list$RS_pval / plot_list$WM_pval / plot_list$FI_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_pFDR
# For Cohen's d plots: same layout.
combined_zval <- y_annotation + (plot_list$RS_delta / plot_list$WM_delta / plot_list$FI_delta) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_zval

# Open a PDF device and print both combined plots (each on a separate page)
pdf(file = output_pdf, width = 8, height = 6)
print(combined_pFDR)
print(combined_zval)
dev.off()

# Create y-axis annotation and wrap as a patchwork element
pdf(file = output_pdf, width = 8, height = 6)
plot_list$delta_delta
dev.off()



########### SUBCORTICAL
results_file =  'res.GroupDiff_NCvsBD_subcorticalROI_REML_standard_metanalysis.csv'
roi_file = "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/roi-names-sub.csv"
atlas_name = "aseg"
output_pdf = 'Metanalysis-results_subcorticalROI_NCvsBD_ME_pFDR.pdf'
analysis = 'standard'

# Read ROI names and regression results with stringsAsFactors = FALSE
rois <- read.csv(file = roi_file, stringsAsFactors = FALSE)
results <- read.csv(file = results_file, stringsAsFactors = FALSE)

# Create y-axis annotation and wrap as a patchwork element
y_label <- textGrob("   Corrected        Standard",
                    rot = 90,  # Rotate text vertically
                    gp = gpar(fontsize = 15, fontface = "bold"))
y_annotation <- wrap_elements(full = y_label)

# Add atlas labels to results by matching on ROI names
results$label <- rois$Atlas[match(results$ROI, rois$ROI)]
l.results = results
l.results$label = gsub(l.results$label,pattern = 'Right',replacement = 'Left')
results = rbind(results, l.results)


# Load the specified atlas from ggsegSchaefer
#data(list = atlas_name, package = "ggsegSchaefer")
atlas_obj <- get(atlas_name)

# Define analysis types to loop over
analysis_types <- unique(results$session)

# Create an empty list to store plots
plot_list <- list()

for(analysis in analysis_types) {
  # Filter results by analysis type
  results_subset <- results[results$session == analysis, ]
  results_subset$pFDR = p.adjust(results_subset$pval, method = 'fdr')
  
  # Subset to relevant columns: label, p_value, and Cohens_d
  pvalue_df <- results_subset[, c("label", "pFDR", "zval")]
  
  
  # Create a copy of the atlas object for merging
  ### 3d
  atlas_temp <- atlas_obj
  
  # Merge the atlas data with the p-values
  atlas_temp$data =  left_join(atlas_temp$data, 
                               pvalue_df, by = "label", 
                               relationship = "many-to-many")
  
  
  # Set p_value > 0.05 to NA and mask Cohen's d accordingly
  atlas_temp$data$pFDR[atlas_temp$data$pFDR > 0.05] <- NA
  atlas_temp$data$zval[is.na(atlas_temp$data$pFDR)] <- NA
  
  # Create p_value plot
  p_pval <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = pFDR),
                  position = "dispersed") +
    scale_fill_gradient(low = "yellow", high = "red", na.value = "grey", limits = c(0, 0.05)) +
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                          axis.text.x=element_text(size=10,  family="sans"))
  # Optionally add ggtitle() if desired
  
  
  # Create Cohen's d plot
  p_zval <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = zval),
                  position = "dispersed") +
    scale_fill_gradient(low = "blue", high = "red", na.value = "grey", limits = c(-10, 10)) +
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                          axis.text.x=element_text(size=10,  family="sans"))
  # Optionally add ggtitle() if
  
  # Store the plots in the list with keys "standard_p_value", etc.
  plot_list[[paste0(analysis, "_pval")]] <- p_pval
  plot_list[[paste0(analysis, "_zval")]] <- p_zval
}

# Combine the plots using patchwork:
# For p_value plots: y_annotation on the left and standard above corrected.
combined_p_value <- y_annotation + (plot_list$RS_pval / plot_list$WM_pval / plot_list$FI_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_p_value
# For Cohen's d plots: same layout.
combined_z_value <- y_annotation + (plot_list$RS_zval / plot_list$WM_zval / plot_list$FI_zval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_z_value

# Combine the plots using patchwork:
# For p_value plots: y_annotation on the left and standard above corrected.
combined_p_value <- y_annotation + (plot_list$standard_pval / plot_list$corrected_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_p_value
# For Cohen's d plots: same layout.
combined_z_value <- y_annotation + (plot_list$standard_zval / plot_list$corrected_zval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_z_value
# Open a PDF device and print both combined plots (each on a separate page)
pdf(file = output_pdf, width = 4, height = 4)
print(combined_p_value)
print(combined_z_value)
dev.off()

#Plot Delta
results_file = 'res.GroupDiff_NCvsBD_subcorticalROI_REML_delta_metanalysis.csv'
output_pdf = 'Metanalysis-results_subcorticalROI_NCvsBD_heterogeneity_pFDR.pdf'

# Read ROI names and regression results with stringsAsFactors = FALSE
rois <- read.csv(file = roi_file, stringsAsFactors = FALSE)
results <- read.csv(file = results_file, stringsAsFactors = FALSE)

# Create y-axis annotation and wrap as a patchwork element
y_label <- textGrob("    FI                   WM                    RS",
                    rot = 90,  # Rotate text vertically
                    gp = gpar(fontsize = 14, fontface = "bold"))
y_annotation <- wrap_elements(full = y_label)

# Add atlas labels to results by matching on ROI names
results$label <- rois$Atlas[match(results$ROI, rois$ROI)]
l.results = results
l.results$label = gsub(l.results$label,pattern = 'Right',replacement = 'Left')
results = rbind(results, l.results)

# Load the specified atlas from ggsegSchaefer
atlas_obj <- get(atlas_name)

# Define analysis types to loop over
analysis_types <- unique(results$session)

# Create an empty list to store plots
plot_list <- list()

for(analysis in analysis_types) {
  # Filter results by analysis type
  results_subset <- results[results$session == analysis, ]
  results_subset$QEpFDR = p.adjust(results_subset$QEp, method = 'fdr')
  results_subset$pFDR = p.adjust(results_subset$pval, method = 'fdr')
  
  # Subset to relevant columns: label, pFDR, and Cohens_d
  pvalue_df <- results_subset[, c("label", "QEpFDR",  "QE", "zval", 'pFDR')]
  
  # Create a copy of the atlas object for merging
  atlas_temp <- atlas_obj
  atlas_temp$data <- left_join(atlas_temp$data, pvalue_df, by = "label", 
                               relationship = "many-to-many")
  
  # Set pFDR > 0.05 to NA and mask Cohen's d accordingly
  atlas_temp$data$pFDR[atlas_temp$data$pFDR > 0.05] <- NA
  atlas_temp$data$QE[is.na(atlas_temp$data$pFDR)] <- NA
  atlas_temp$data$QEpFDR[is.na(atlas_temp$data$pFDR)] <- NA
  atlas_temp$data$zval[is.na(atlas_temp$data$pFDR)] <- NA
  
  # Create pFDR plot
  p_pFDR <- ggseg(atlas = atlas_temp,
                  colour = "grey37", size = 0.3,
                  mapping = aes(fill = QE),
                  position = "dispersed") +
    scale_fill_gradient(low = "white", high = "darkblue", na.value = "grey", limits = c(0, 10)) +
    theme_minimal()
  # Optionally add ggtitle() if desired
  
  # Create Cohen's d plot
  p_delta <- ggseg(atlas = atlas_temp,
                   colour = "grey37", size = 0.3,
                   mapping = aes(fill = zval),
                   position = "dispersed") +
    scale_fill_gradient(low='white', high = "darkgreen", na.value = "grey" , limits = c(0, 6)) + #
    theme_brain() + theme(text=element_text(size=16,  family="sans"), 
                          axis.text.x=element_text(size=10,  family="sans"))
  # Optionally add ggtitle() if desired
  
  # Store the plots in the list with keys "standard_pFDR", etc.
  plot_list[[paste0(analysis, "_pval")]] <- p_pFDR
  plot_list[[paste0(analysis, "_delta")]] <- p_delta
}



# Combine the plots using patchwork:
# For pFDR plots: y_annotation on the left and standard above corrected.
combined_pFDR <- y_annotation + (plot_list$RS_pval / plot_list$WM_pval / plot_list$FI_pval) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_pFDR
# For Cohen's d plots: same layout.
combined_zval <- y_annotation + (plot_list$RS_delta / plot_list$WM_delta / plot_list$FI_delta) +
  plot_layout(guides = "collect", widths = c(1, 10)) + plot_annotation(title = '')
combined_zval

# Open a PDF device and print both combined plots (each on a separate page)
pdf(file = output_pdf, width = 8, height = 6)
print(combined_pFDR)
print(combined_zval)
dev.off()



#  association plot Delta ~ QE
results_file = 'res.NCvsBD_GroupDiff_NC_corticalROI_delta_metanalysis.csv'
results <- read.csv(file = results_file, stringsAsFactors = FALSE)
results_file = 'res.NCvsBD_GroupDiff_NC_subcorticalROI_delta_metanalysis.csv'
results <- rbind(results, read.csv(file = results_file, stringsAsFactors = FALSE))

results$pFDR = NA
results$pFDR[which(results$session=='RS')] = p.adjust(results$pval[which(results$session=='RS')], method = 'fdr')
results$pFDR[which(results$session=='WM')] = p.adjust(results$pval[which(results$session=='WM')], method = 'fdr')
results$pFDR[which(results$session=='FI')] = p.adjust(results$pval[which(results$session=='FI')], method = 'fdr')

ggplot(data = results[which(results$Delta>0),], aes(x = Delta, y = QE)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "red", linetype = "dashed") +
  facet_wrap(~ session) +
  stat_cor(method = "pearson",
           #aes(label = paste(after_stat(rr.label), ", p = ", after_stat(p.label))),
           label.x.npc = "left",  # positions label relative to each panel (left)
           label.y.npc = "top",   # positions label relative to each panel (top)
           size = 4,
           parse = FALSE) +       # disable parsing to prevent unexpected comma errors
  theme_pubclean(base_size = 14) +
  labs(title = "Delta by QE (varying ROIs)",
       x = "Delta",
       y = "QE")


df.bd  = results
results_file = 'res.GroupDiff_NCvsSCZ_corticalROI_delta_metanalysis.csv'
results <- read.csv(file = results_file, stringsAsFactors = FALSE)
results_file = 'res.GroupDiff_NCvsSCZ_subcorticalROI_delta_metanalysis.csv'
results <- rbind(results, read.csv(file = results_file, stringsAsFactors = FALSE))
results$pFDR = NA
results$pFDR[which(results$session=='RS')] = p.adjust(results$pval[which(results$session=='RS')], method = 'fdr')
results$pFDR[which(results$session=='WM')] = p.adjust(results$pval[which(results$session=='WM')], method = 'fdr')
results$pFDR[which(results$session=='FI')] = p.adjust(results$pval[which(results$session=='FI')], method = 'fdr')

df.scz = results

df.bd$Diagnosis = 'BD'
names(df.bd) = paste0(names(df.bd),'_BD')
df.scz$Diagnosis = 'SCZ'
names(df.scz) = paste0(names(df.scz),'_SCZ')
df = cbind(df.scz, df.bd)

summary(lm(df$Delta_SCZ~df$Delta_BD))

a = ggplot(data = df, aes(x = Delta_SCZ, y = Delta_BD)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "darkgreen") +
  stat_cor(method = "pearson",
           #aes(label = paste(after_stat(rr.label), ", p = ", after_stat(p.label))),
           label.x.npc = "left",  # positions label relative to each panel (left)
           label.y.npc = "top",   # positions label relative to each panel (top)
           size = 4,
           parse = FALSE) +       # disable parsing to prevent unexpected comma errors
  theme_pubclean(base_size = 14) + ylim(-10,10) + xlim(-10,10) +
  labs(title = "SCZ ~ BD result comparisons across ROIs",
       x = "Delta zscore - SCZ",
       y = "Delta zscore - BD")
a

b = ggplot(data = df, aes(x = zval_SCZ, y = zval_BD)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "darkgreen") +
  stat_cor(method = "pearson",
           #aes(label = paste(after_stat(rr.label), ", p = ", after_stat(p.label))),
           label.x.npc = "left",  # positions label relative to each panel (left)
           label.y.npc = "top",   # positions label relative to each panel (top)
           size = 4,
           parse = FALSE) +       # disable parsing to prevent unexpected comma errors
  theme_pubclean(base_size = 14) + ylim(-10,10) + xlim(-10,10) +
  labs(
       x = "Meta zscore - SCZ",
       y = "Meta zscore - BD")
b

a + b
