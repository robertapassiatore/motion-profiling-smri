# Required libraries
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("multivarious", quietly = TRUE)) install.packages("multivarious")

library(readr)
library(dplyr)
library(parallel)
library(ggplot2)
library(lme4)
#library(lmerTest)
library(reshape2)
library(multivarious)

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# List of datasets
setwd("/documents/Roberta/MotionProfiling/NC_LinearMotionEffects")
datasets <- list.files(pattern = ".csv")
res.dir = "/dcouments/Roberta/MotionProfiling/NC_LinearMotionEffects/results/"

# Function to extract dataset name without extension
get_dataset_name <- function(file) {
  sub("\\.csv$", "", basename(file))
}


# Analysis function
run_analysis <- function(data, dataset_name, start_col, end_col, suffix) {
  results <- data.frame(Dataset = character(),
                        ROI = character(),
                        Analysis = character(),
                        p_value = numeric(),
                        R = numeric(),
                        SE = numeric(),
                        Empirical_p_value = numeric(),
                        #Cohens_d = numeric(),
                        Avg_R = numeric(),
                        stringsAsFactors = FALSE)
  
  model_specs <- list(
    allcov = ~ PC1 + PC2 + PC3 + PC4 + PC5, #sex + age + agesq + TIV 
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
      p_value_original = lmp(model)
      SE = model_summary$sigma
      R = model_summary$adj.r.squared
      
      if (model_type == "allcov") {
        permuted_Rs <- numeric(5000)
        #cohen_d_values <- numeric(5000)
        
        for (j in 1:5000) {
          shuffled_data <- data
          shuffled_data[c("PC1",'PC2','PC3','PC4',"PC5")] <- lapply(shuffled_data[c("PC1",'PC2','PC3','PC4',"PC5")], sample)
          
          perm_model <- lm(formula, data = shuffled_data)
          perm_R <- summary(perm_model)$adj.r.squared
          
          permuted_Rs[j] <- perm_R
        }
        
        empirical_p_value <- mean(abs(permuted_Rs) >= abs(R))
        avg_Rs <- mean(permuted_Rs)
      } else {
        empirical_p_value <- NA
        avg_Rs <- NA
      }
      
      results <- rbind(results, data.frame(Dataset = dataset_name,
                                           ROI = roi_name,
                                           Analysis = model_type,
                                           p_value = p_value_original,
                                           R = R,
                                           SE = SE,
                                           Empirical_p_value = empirical_p_value,
                                           #Cohens_d = cohens_d,
                                           Avg_R = avg_Rs))
      
    }
  }
  
  results$pFDR <- p.adjust(results$p_value, method = "fdr")
  results = results[-which(results$Analysis=='onecov'),]
  
  # Write results to CSV file in append mode
  write.table(results, paste0(res.dir, "regression_results_allPC_permutations_", suffix, dataset_name, ".csv"), 
              row.names = FALSE, col.names = !file.exists(paste("regression_results_allPC_permutations_", suffix, ".csv", sep="")), 
              sep = ",", append = TRUE)
  
}

# Loop through each dataset
for (dataset in datasets) {
  #dataset = datasets[1]
  data <- read_csv(dataset)
  dataset_name <- get_dataset_name(dataset)
  data$group = as.factor(data$group)
  data$sex = as.factor(data$sex)
  data$ID = as.factor(data$ID)
  
  # Standardize numeric columns
  numeric_columns <- sapply(data, is.numeric)
  data_standardized <- data
  data_standardized[,c(12:119)] = residualize(~sex + age + agesq + TIV, data_standardized[,c(12:119)], data_standardized) 
  data_standardized[numeric_columns] <- scale(data[numeric_columns])
  #data_melt = melt(data_standardized, measure.vars = paste0('PC',1:5), variable.name='PC.ID' , value.name =  'PC')
  #data_melt$Subject = rep(1:length(data_standardized$ID),5)
  #data_standardized$PC = rowMeans(data_standardized[,c(7:11)])
  
  
  # Run analysis for cortical ROIs
  run_analysis(data_standardized, dataset_name, 12, 111, "corticalROI")
  
  # Run analysis for subcortical ROIs
  run_analysis(data_standardized, dataset_name, 112, 119, "subcorticalROI")
}

