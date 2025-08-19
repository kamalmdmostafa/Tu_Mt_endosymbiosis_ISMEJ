# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(pcaMethods)  # For SVD imputation

# Function to perform PCA with SVD imputation
perform_svd_pca <- function(data) {
  data <- as.data.frame(data)
  gene_names <- rownames(data)
  
  # Transpose the data to have samples as rows
  data_t <- t(data)
  sample_names <- rownames(data_t)
  
  # Perform SVD imputation PCA
  pca_result <- pca(data_t, method = "svdImpute", nPcs = 2)
  pca_df <- as.data.frame(pca_result@scores)
  colnames(pca_df) <- paste0("PC", 1:ncol(pca_df))
  pca_df$sample <- sample_names
  return(list(pca_result = pca_result, pca_df = pca_df))
}

# Function to plot PCA results
plot_pca <- function(pca_df, pca_result, title) {
  # Extract conditions and time points from sample names
  pca_df$condition <- ifelse(grepl("^g", pca_df$sample), "Symbiotic", "Aposymbiotic")
  pca_df$time_point <- sub(".*_(AM\\d+|PM\\d+)_.*", "\\1", pca_df$sample)
  
  # Set the order of time points
  pca_df$time_point <- factor(pca_df$time_point, levels = c("AM8", "AM11", "PM2", "PM5", "PM8", "PM11"))
  
  # Create a unique identifier for each time point and condition combination
  pca_df$shape_group <- interaction(pca_df$condition, pca_df$time_point, drop = TRUE)
  
  # Create a vector of shape values for symbiotic (filled) and aposymbiotic (hollow) conditions
  symbiotic_shapes <- c(16, 17, 15, 18, 14, 19)  # Filled shapes
  aposymbiotic_shapes <- c(1, 2, 0, 5, 6, 4)     # Hollow shapes
  
  # Assign shapes based on the condition and time point
  levels_symbiotic <- unique(pca_df$shape_group[pca_df$condition == "Symbiotic"])
  levels_aposymbiotic <- unique(pca_df$shape_group[pca_df$condition == "Aposymbiotic"])
  
  shape_values <- c(setNames(symbiotic_shapes[1:length(levels_symbiotic)], levels_symbiotic),
                    setNames(aposymbiotic_shapes[1:length(levels_aposymbiotic)], levels_aposymbiotic))
  
  # Extract variance explained and format the axis labels
  R2 <- pca_result@R2
  variance_explained <- round(100 * R2, 2)
  x_label <- paste0("PC1 (", variance_explained[1], "%)")
  y_label <- paste0("PC2 (", variance_explained[2], "%)")
  
  # Create a color palette for time points
  unique_time_points <- unique(pca_df$time_point)
  color_palette <- scales::hue_pal()(length(unique_time_points))
  
  # Create the plot
  ggplot(pca_df, aes(x = PC1, y = PC2, label = sample)) +
    geom_point(aes(shape = shape_group, color = time_point), 
               size = 6) +
    geom_text_repel(hjust = 0.5, vjust = -1, size = 3, max.overlaps = Inf) +
    scale_shape_manual(values = shape_values, name = "Time Point & Condition") +
    scale_color_manual(values = setNames(color_palette, unique_time_points), name = "Time Point") +
    theme_minimal() +
    labs(title = title, x = x_label, y = y_label) +
    theme(legend.position = "right") +
    guides(shape = guide_legend(override.aes = list(size = 6)),
           color = guide_legend(override.aes = list(shape = 16, size = 6)))
}

# Load your data
TPM_of_Sig_DEG <- read.csv("TPM.csv", row.names = 1)

# Perform SVD PCA
results_svd <- perform_svd_pca(TPM_of_Sig_DEG)

# Plot the results
plot_pca(results_svd$pca_df, results_svd$pca_result, "SVD PCA")

print("Analysis complete. Use RStudio's plot export options to save your plot.")