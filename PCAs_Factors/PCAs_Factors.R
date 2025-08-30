# Check for and install required packages
packages <- c("ggplot2", "gridExtra", "patchwork", "reshape2", "Hmisc", "dplyr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load necessary libraries
library(ggplot2)
library(gridExtra)
library(patchwork)
library(reshape2)
library(Hmisc)
library(dplyr)


####################################################################################
####################################################################################
####################################################################################
### Scree plots
####################################################################################
####################################################################################
####################################################################################
create_double_axis_plot <- function(data, title) {
  max_eigen <- max(data$Eigenvalue, na.rm = TRUE)
  max_cumvar <- max(data$CumulativeVariance, na.rm = TRUE)
  scale_factor <- max_cumvar / max_eigen
  
  p <- ggplot(data, aes(x = PrincipalComponent)) + 
    geom_line(aes(y = Eigenvalue, group = Scan, linetype = Scan)) +
    geom_line(aes(y = CumulativeVariance / scale_factor, group = Scan, linetype = Scan), color = "black") +
    scale_linetype_manual(values=c("solid", "dotted", "dotdash")) +
    scale_color_manual(values=c("black", "black", "black")) +
    scale_x_continuous(breaks=c(0, 4, 8, 12)) +
    scale_y_continuous("Eigenvalue", limits = c(0, max_eigen),
                       sec.axis = sec_axis(~ . * scale_factor, name = "Cumulative Variance")) +
    labs(title = title, x = "Principal Component", y = "Eigenvalue") +
    theme_bw() +
    theme(plot.title = element_text(size = 20,face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.text.y.right = element_text(size = 18),
          axis.title.y.right = element_text(size = 20),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = "right")
  return(p)
}


# Read data for ScreePlots
data1 <- read.csv("UNIBA1_Eigenvalues.Variance.txt")
data2 <- read.csv("UNIBA2_Eigenvalues.Variance.txt")
data3 <- read.csv("LIBD_Eigenvalues.Variance.txt")

# Create individual Scree plots
p1 <- create_double_axis_plot(data1, "UNIBA1")
p2 <- create_double_axis_plot(data2, "UNIBA2")
p3 <- create_double_axis_plot(data3, "LIBD")

# Combine the Scree plots into one figure
combined_scree_plots <- grid.arrange(p1, p2, p3, ncol = 1)
#ggsave("Figure1_ScreePlots.png", combined_scree_plots, width = 12, height = 18, units = "in")

####################################################################################
####################################################################################
####################################################################################
### Loadings for all the PCAs
####################################################################################
####################################################################################
####################################################################################
pcs_data_uniba1 <- read.csv("UNIBA1_PCs.txt")
pcs_data_uniba1$Scan <- factor(pcs_data_uniba1$Scan, levels = c("NB", "FE", "FI"))

pcs_data_uniba2 <- read.csv("UNIBA2_PCs.txt")
pcs_data_uniba2$Scan <- factor(pcs_data_uniba2$Scan, levels = c("NB", "RS", "FI"))

pcs_data_libd <- read.csv("LIBD_PCs.txt")
pcs_data_libd$Scan <- factor(pcs_data_libd$Scan, levels = c("NB", "RS", "FI"))

heatmap_uniba1 <- ggplot(pcs_data_uniba1, aes(x = PC, y = MotionParameter, fill = Loadings)) +
  geom_tile() +
  scale_fill_gradient2(low = "green", high = "purple", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  facet_wrap(~Scan) +
  theme_minimal() +
  labs(title = "UNIBA1", x = "PC", y = "MotionParameter", fill = "Loadings")+
  theme(plot.title = element_text(size = 20,face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 270, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y.right = element_text(size = 18),
        axis.title.y.right = element_text(size = 20),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right")

heatmap_uniba2 <- ggplot(pcs_data_uniba2, aes(x = PC, y = MotionParameter, fill = Loadings)) +
  geom_tile() +
  scale_fill_gradient2(low = "green", high = "purple", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  facet_wrap(~Scan) +
  theme_minimal() +
  labs(title = "UNIBA2", x = "PC", y = "MotionParameter", fill = "Loadings")+
  theme(plot.title = element_text(size = 20,face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 270, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y.right = element_text(size = 18),
        axis.title.y.right = element_text(size = 20),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right")

heatmap_libd <- ggplot(pcs_data_libd, aes(x = PC, y = MotionParameter, fill = Loadings)) +
  geom_tile() +
  scale_fill_gradient2(low = "green", high = "purple", mid = "white", midpoint = 0, limit = c(-1, 1)) +
  facet_wrap(~Scan) +
  theme_minimal() +
  labs(title = "LIBD", x = "PC", y = "MotionParameter", fill = "Loadings")+
  theme(plot.title = element_text(size = 20,face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 270, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.y.right = element_text(size = 18),
        axis.title.y.right = element_text(size = 20),
        strip.text = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right")

# Example adjustment for the heatmap of UNIBA1
heatmap_uniba1 <- heatmap_uniba1 +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5))

# Example adjustment for the heatmap of UNIBA2
heatmap_uniba2 <- heatmap_uniba2 +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5))

# Example adjustment for the heatmap of LIBD
heatmap_libd <- heatmap_libd +
  guides(fill = guide_colorbar(barwidth = 1, barheight = 5))

combined_heatmaps <- (heatmap_uniba1 / heatmap_uniba2 / heatmap_libd)
#ggsave("Figure1_Loadings.png", combined_heatmaps, width = 10, height = 15, units = "in", dpi = 300)

####################################################################################
####################################################################################
####################################################################################
### Correlation between factor loadings
####################################################################################
####################################################################################
####################################################################################

# Function to calculate correlations and p-values, adjusted for FDR
calculate_cor_pvals <- function(data) {
  # Convert the data to numeric, handling non-numeric data gracefully
  data <- data.frame(lapply(data, function(x) as.numeric(as.character(x))))
  if (any(sapply(data, is.na))) {
    stop("Data contains non-numeric values that could not be converted to numeric.")
  }
  
  # Create an empty data frame to store the results
  results <- expand.grid(Var1 = names(data), Var2 = names(data))
  results$Correlation <- NA
  results$P_value <- NA
  results$Adjusted_P_value <- NA  # Add a column for adjusted p-values
  results$Sig <- ""  # Initialize with an empty string
  
  # Calculate the correlations and p-values efficiently
  for (i in seq_len(ncol(data))) {
    for (j in i:ncol(data)) {
      test_result <- cor.test(data[[i]], data[[j]], method = "pearson")
      # Fill both the upper and the lower triangles
      idx <- results$Var1 == names(data)[i] & results$Var2 == names(data)[j]
      idx_symmetric <- results$Var1 == names(data)[j] & results$Var2 == names(data)[i]
      results$Correlation[idx] <- test_result$estimate
      results$Correlation[idx_symmetric] <- test_result$estimate
      results$P_value[idx] <- test_result$p.value
      results$P_value[idx_symmetric] <- test_result$p.value
    }
  }
  
  # Adjust the p-values for multiple comparisons using FDR
  results$Adjusted_P_value <- p.adjust(results$P_value, method = "BH")
  
  # Mark significant correlations based on adjusted p-values
  results$Sig[results$Adjusted_P_value < 0.05] <- "*"
  
  return(results)
}

# Function to plot correlation matrix with significant markings
plot_cor_matrix <- function(data, title) {
  # Ensure variables are in the right format for comparison
  data$Var1 <- as.character(data$Var1)
  data$Var2 <- as.character(data$Var2)
  
  # Filter the data to keep only the lower left triangle
  lower_triangle <- data %>%
    filter(Var1 > Var2)
  
  ggplot(lower_triangle, aes(x = Var2, y = Var1, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Sig), vjust = 0.5, hjust = 0.5, size = 10, color = "black", data = subset(lower_triangle, Sig != "")) +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white", midpoint = 0, 
      limit = c(-1, 1), space = "Lab", name = "Pearson r",
      guide = guide_colourbar(
        title.position = "top", title.hjust = 0.5, 
        label.position = "right", barwidth = 1, barheight = 5
      )
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_blank(),
      axis.text.x.top = element_text(angle = 270, vjust = 0.4, size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.x = element_blank(),
      axis.ticks.x.top = element_blank(),
      axis.title.x = element_blank(),
      axis.title.x.top = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 20, r = 10, b = 5, l = 10, unit = "pt")
    ) +
    scale_x_discrete(position = "top") +
    coord_fixed(ratio = 1)
}

# Process datasets and generate plots using FDR-corrected p-values
# You would need to set up a list of datasets and their titles and file paths as needed.
datasets <- list(
  uniba1 = list(file = "UNIBA1_Loadings.Corr.txt", 
                order = c("NB.PC1", "NB.PC2", "NB.PC3", "NB.PC4", "NB.PC5", 
                          "FE.PC1", "FE.PC2", "FE.PC3", "FE.PC4", "FE.PC5", 
                          "FI.PC1", "FI.PC2", "FI.PC3", "FI.PC4", "FI.PC5"),
                title = "UNIBA1"),
  uniba2 = list(file = "UNIBA2_Loadings.Corr.txt", 
                order = c("NB.PC1", "NB.PC2", "NB.PC3", "NB.PC4", "NB.PC5",
                          "RS.PC1", "RS.PC2", "RS.PC3", "RS.PC4", "RS.PC5",
                          "FI.PC1", "FI.PC2", "FI.PC3", "FI.PC4", "FI.PC5"),
                title = "UNIBA2"),
  libd = list(file = "LIBD_Loadings.Corr.txt", 
              order = c("NB.PC1", "NB.PC2", "NB.PC3", "NB.PC4", "NB.PC5",
                        "RS.PC1", "RS.PC2", "RS.PC3", "RS.PC4", "RS.PC5",
                        "FI.PC1", "FI.PC2", "FI.PC3", "FI.PC4", "FI.PC5"),
              title = "LIBD")
)

# Generate and save plots
plot_list <- lapply(datasets, function(ds) {
  data <- read.csv(ds$file, header = TRUE)
  selected_data <- data[, ds$order]
  cor_pvals_data <- calculate_cor_pvals(selected_data)
  plot_cor_matrix(cor_pvals_data, ds$title)
})


# Assign names to the list of plots for easy access
names(plot_list) <- c("plot_uniba1", "plot_uniba2", "plot_libd")


# Combine all plots into a single 3x3 grid
all_plots_combined <- p1/ heatmap_uniba1 / plot_list[["plot_uniba1"]] / 
  p2 / heatmap_uniba2 / plot_list[["plot_uniba2"]] / 
  p3 / heatmap_libd / plot_list[["plot_libd"]] +
  plot_layout(ncol = 3, nrow = 3)

# Save the combined plot to a file
ggsave("Figure1_Complete_Panels.png", all_plots_combined, width = 20, height = 20, units = "in")