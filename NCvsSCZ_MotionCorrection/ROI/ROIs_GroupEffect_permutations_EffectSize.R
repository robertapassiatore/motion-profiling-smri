# Required libraries
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(readr)
library(dplyr)
library(parallel)
library(ggplot2)

# List of datasets
setwd("~/Documents/Motion-Profiling/analyses/NC_Percentiles_GroupDifferences")
datasets <- list.files(pattern = ".csv")
res.dir = "/Users/roberta.passiatore/Documents/Motion-Profiling/analyses/NC_Percentiles_GroupDifferences/results/"

# Function to extract dataset name without extension
get_dataset_name <- function(file) {
  sub("\\.csv$", "", basename(file))
}

# Data frame to store correlations
correlation_results <- data.frame(Dataset = character(),
                                  AnalysisType = character(),
                                  Correlation = numeric(),
                                  stringsAsFactors = FALSE)

# Analysis function
run_analysis <- function(data, dataset_name, start_col, end_col, suffix) {
  results <- data.frame(Dataset = character(),
                        ROI = character(),
                        Analysis = character(),
                        p_value = numeric(),
                        Beta_Coefficient_Group = numeric(),
                        Empirical_p_value = numeric(),
                        Cohens_d = numeric(),
                        T_value = numeric(),
                        SE=numeric(),
                        Avg_Cohens_d = numeric(),
                        stringsAsFactors = FALSE)
  
  model_specs <- list(
    standard = ~ group + sex + age + agesq + TIV,
    corrected = ~ group + sex + age + agesq + TIV + PC1 + PC2 + PC3 + PC4 + PC5
  )
  
  standard_cohens_d <- numeric()
  avg_cohens_d_corrected <- numeric()
  
  # Analysis loop
  for (i in start_col:end_col) {
    roi_name <- colnames(data)[i]
    
    for (model_type in names(model_specs)) {
      formula <- as.formula(paste(roi_name, deparse(model_specs[[model_type]])))
      model <- lm(formula, data = data)
      model_summary <- summary(model)
      beta_coefficients <- coef(model_summary)["group2", "Estimate"]
      p_value_original <- coef(model_summary)["group2", "Pr(>|t|)"]
      t_value_original = coef(model_summary)["group2", "t value"]
      se_original = coef(model_summary)["group2", "Std. Error"]
      
      # Calculate standard deviation of residuals for Cohen's d
      std_dev <- sd(residuals(model))
      cohens_d <- beta_coefficients / std_dev
      
      
      if (model_type == "corrected") {
        permuted_coefficients <- numeric(50) #5000
        cohen_d_values <- numeric(50) #5000
        
        for (j in 1:50) {
          shuffled_data <- data
          shuffled_data[c("PC1", "PC2", "PC3", "PC4", "PC5")] <- lapply(shuffled_data[c("PC1", "PC2", "PC3", "PC4", "PC5")], sample)
          
          perm_model <- lm(formula, data = shuffled_data)
          perm_coefficients <- coef(summary(perm_model))["group2", "Estimate"]
          
          permuted_coefficients[j] <- perm_coefficients
          cohen_d_values[j] <- perm_coefficients / std_dev
        }
        
        empirical_p_value <- mean(abs(permuted_coefficients) >= abs(beta_coefficients))
        avg_cohens_d <- mean(cohen_d_values)
      } else {
        empirical_p_value <- NA
        avg_cohens_d <- NA
      }
      
      results <- rbind(results, data.frame(Dataset = dataset_name,
                                           ROI = roi_name,
                                           Analysis = model_type,
                                           p_value = p_value_original,
                                           Beta_Coefficient_Group = beta_coefficients,
                                           Empirical_p_value = empirical_p_value,
                                           Cohens_d = cohens_d,
                                           T_value = t_value_original,
                                           SE=se_original,
                                           Avg_Cohens_d = avg_cohens_d))
      
      if (model_type == "standard") {
        standard_cohens_d <- c(standard_cohens_d, cohens_d)
      }
      
      if (model_type == "corrected") {
        avg_cohens_d_corrected <- c(avg_cohens_d_corrected, avg_cohens_d)
      }
    }
  }
  
  results$pFDR <- p.adjust(results$p_value, method = "fdr")
  results$N = length(data$sex)
  
  
  results = results[-which(results$Analysis=='corrected'),]
  
  # Write results to CSV file in append mode
  write.table(results, paste0(res.dir,"regression_results_permutations_", suffix, dataset_name, ".csv"), 
              row.names = FALSE, col.names = !file.exists(paste("regression_results_permutations_", suffix, ".csv", sep="")), 
              sep = ",", append = TRUE)
  
  # Correlation analysis between Cohen's d from standard model and Avg Cohen's d from corrected model
  #cor_results <- cor(standard_cohens_d, avg_cohens_d_corrected, use = "complete.obs")
  #print(paste("Correlation between Cohen's d from standard model and average Cohen's d from permutations for", dataset_name, suffix, ":", cor_results))
  
  #correlation_results <- rbind(correlation_results, data.frame(Dataset = dataset_name, AnalysisType = suffix, Correlation = cor_results))
  
  # Plot distribution of Empirical p-values
  #empirical_p_values_filtered <- results$Empirical_p_value[results$Analysis == "corrected"]
  #hist(empirical_p_values_filtered, main = paste("Distribution of Empirical p-values \n for", dataset_name, suffix), xlab = "Empirical p-value", breaks = 50, col = "blue")
}

# Loop through each dataset
for (dataset in datasets) {
  data <- read_csv(dataset)
  data = data[,-1]
  dataset_name <- get_dataset_name(dataset)
  data$group = as.factor(data$group)
  data$sex = as.factor(data$sex)
  
  # Standardize numeric columns
  numeric_columns <- sapply(data, is.numeric)
  data_standardized <- data
  data_standardized[numeric_columns] <- scale(data[numeric_columns])
  
  
  # Run analysis for cortical ROIs
  run_analysis(data_standardized, dataset_name, 12, 111, "corticalROI")
  
  # Run analysis for subcortical ROIs
  run_analysis(data_standardized, dataset_name, 112, 119, "subcorticalROI")
}

#data, dataset_name, start_col, end_col, suffix
# Save the correlation results to a CSV file
#write.csv(correlation_results, paste0("correlation_results_",dataset_name,".csv"), row.names = FALSE)


# # Create the boxplot with facet wrapping
# p <- ggplot(data, aes(x = IQR_rating, y = PC5, fill = group)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.7) +         # Hide outlier points; adjust transparency
#   geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +       # Add jittered data points for clarity
#   facet_wrap(~ factor(group, levels = c("1", "2"), labels = c("BD", "NC")), nrow = 1) +                         # Wrap by the 'measure' variable
#   #scale_fill_brewer(palette = "Set1") +                     # Use a ColorBrewer palette
#   #labs(title = "Boxplot by Group", x = "Group", y = "Value") +
#   theme_minimal(base_size = 14) +                           # Use a clean classic theme
#   theme(
#     legend.position = "none",                             # Remove legend if redundant
#     strip.background = element_rect(fill = "white", color = "black"),
#     strip.text = element_text(face = "bold")
#   )
# 
# # Print the plot
# print(p)

# List files in the directory that contain "FACES" in the name
files_to_rename <- list.files(path = res.dir, pattern = "SCAP", full.names = TRUE)

# Loop through each file and rename it
for (old_file in files_to_rename) {
  # Create the new file name by replacing "FACES" with "FI"
  new_file <- gsub("SCAP", "WM", old_file)
  
  # Rename the file
  success <- file.rename(old_file, new_file)
  
  if (success) {
    cat("Renamed:", old_file, "to", new_file, "\n")
  } else {
    cat("Failed to rename:", old_file, "\n")
  }
}
