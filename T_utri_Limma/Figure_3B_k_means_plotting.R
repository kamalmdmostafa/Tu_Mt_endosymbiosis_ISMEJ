# Load necessary libraries
requiredPackages <- c(
  "ComplexHeatmap", "dplyr", "circlize", "RColorBrewer", 
  "reshape2", "gridExtra", "grid", "ggplot2", 
  "data.table", "ClusterR", "writexl", "doParallel", 
  "parallel", "scales"
)
newPackages <- requiredPackages[!(requiredPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)
lapply(requiredPackages, library, character.only = TRUE)

# Placeholder for actual data loading
# Load your dataset here
data <- KMC_input_from_IDEP # Adjust this line to load your dataset correctly
dataframeName <- "Kallisto_four_clusters" # Set once to maintain consistency

# Specify the number of clusters you want
specified_clusters <- 4  # Change this number as needed

# Choose light condition pattern
# Set this to TRUE for the 12-sample pattern, FALSE for the 6-sample pattern
use_12_sample_pattern <- TRUE  # Change this to FALSE if you want to use the 6-sample pattern

# Add light condition information based on the chosen pattern
if (use_12_sample_pattern) {
  light_cond <- c("Dark", "Light", "Light", "Light", "Light", "Dark", "Dark", "Light", "Light", "Light", "Light", "Dark")
} else {
  light_cond <- c("Dark", "Light", "Light", "Light", "Light", "Dark")
}

# Prepare data for clustering
expression_data <- as.matrix(data[, -which(names(data) %in% c("Gene.Name", "Cluster"))])
rownames(expression_data) <- data$Gene.Name

# Function to perform fast K-means clustering
improved_kmeans <- function(data, k, nstart = 30, max_iters = 300, sample_size = NULL) {
  set.seed(123)  # for reproducibility
  
  if (is.null(sample_size)) {
    sample_size <- min(5000, nrow(data))
  }
  
  initial_centers <- data[sample(nrow(data), sample_size), ]
  
  km <- KMeans_rcpp(data, clusters = k, num_init = nstart, 
                    max_iters = max_iters, initializer = "kmeans++",
                    verbose = FALSE)
  
  return(km)
}

# Determine number of clusters and perform clustering
num_clusters <- specified_clusters

print(paste("Performing K-means clustering with", num_clusters, "clusters"))

# Perform fast K-means clustering
kmeans_result <- improved_kmeans(expression_data, num_clusters, nstart = 30, max_iters = 300)

# Add the new cluster assignments to the data
data$NewCluster <- kmeans_result$clusters

# Convert data to data.table
dt <- as.data.table(data)

# Ensure NewCluster is a factor with levels in the correct order
dt$NewCluster <- factor(dt$NewCluster, levels = 1:num_clusters)

# Reorder data based on new cluster assignments
dt_long <- melt(dt, id.vars = c("Gene.Name", "Cluster", "NewCluster"), 
                variable.name = "Sample", value.name = "Expression")

# Reorder genes within clusters
dt_long[, PeakValue := max(abs(Expression)), by = .(Gene.Name, NewCluster)]
setorder(dt_long, NewCluster, -PeakValue)
reordered_gene_order <- unique(dt_long[, .(Gene.Name, NewCluster)])$Gene.Name

# Reconstruct the heatmap_data with reordered genes
heatmap_data <- expression_data[reordered_gene_order, ]

# Create annotation for clusters
cluster_anno <- dt$NewCluster[match(rownames(heatmap_data), dt$Gene.Name)]
names(cluster_anno) <- rownames(heatmap_data)

# Create a vector for column split
column_split <- ifelse(grepl("^g", colnames(heatmap_data)), "Symbiotic", "Aposymbiotic")

# Ensure "Symbiotic" is placed before "Aposymbiotic"
column_split_factor <- factor(column_split, levels = c("Symbiotic", "Aposymbiotic"))

# Choose palette
palette_name <- "RdBu"
brewer_palette <- rev(brewer.pal(11, palette_name))

# Generate color-blind friendly palette for clusters
generate_cb_friendly_palette <- function(palette_name, n) {
  base_palette <- brewer.pal(9, palette_name)
  colorRampPalette(c(base_palette[c(1,3,5,7,9)]))(n)
}

cluster_colors <- generate_cb_friendly_palette(palette_name, num_clusters)
names(cluster_colors) <- 1:num_clusters

# Light condition colors
light_cond_colors <- c("Dark" = brewer_palette[2], "Light" = brewer_palette[10])

# Use RColorBrewer for column annotation
col_anno <- HeatmapAnnotation(
  CellType = column_split_factor,
  LightCond = factor(light_cond, levels = c("Dark", "Light")),
  col = list(
    CellType = c("Symbiotic" = brewer_palette[3], "Aposymbiotic" = brewer_palette[9]),
    LightCond = light_cond_colors
  ),
  annotation_height = unit(c(4, 1), "cm")
)

# Define RColorBrewer colors for heatmap
brewer_colors <- colorRamp2(
  seq(-2, 2, length.out = 11),
  brewer_palette
)

# Calculate gene counts for each cluster
gene_count <- dt[, .(GeneCount = .N), by = NewCluster]
setorder(gene_count, NewCluster)  # Ensure gene_count is in the correct order

# Create a legend for clusters with gene counts
legend_labels <- paste0("Cluster ", gene_count$NewCluster, " (n=", gene_count$GeneCount, ")")
legend_annotations <- Legend(
  labels = legend_labels,
  legend_gp = gpar(fill = cluster_colors, col = "black"),
  title = "Clusters",
  title_position = "topcenter",
  grid_height = unit(0.5, "cm"),
  grid_width = unit(0.5, "cm")
)

# Create number labels without "Cluster" for left annotation
number_labels_left <- as.character(unique(cluster_anno))

# Create the heatmap with RColorBrewer colors
heatmap <- Heatmap(
  heatmap_data,
  name = "Expression",
  col = brewer_colors,
  show_row_names = FALSE,
  show_column_names = TRUE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = factor(cluster_anno, levels = unique(cluster_anno)),
  column_split = column_split_factor,
  top_annotation = col_anno,
  row_title = "Genes",
  column_title = "Samples",
  heatmap_legend_param = list(
    title = "Z-scaled Expression",
    title_position = "leftcenter-rot",
    legend_direction = "vertical",
    legend_height = unit(4, "cm"),
    legend_width = unit(1, "cm")
  ),
  row_gap = unit(2, "mm"),
  column_gap = unit(3, "mm"),
  width = unit(8, "cm"), # Adjust width to make the heatmap narrower
  height = unit(12, "cm"), # Adjust height to make the heatmap shorter
  use_raster = TRUE,
  raster_quality = 4,
  left_annotation = rowAnnotation(
    Cluster = anno_block(
      gp = gpar(fill = cluster_colors, col = "black"), 
      labels_gp = gpar(col = "black", fontsize = 12),
      labels = number_labels_left,
      labels_rot = 90 # Rotate labels to be vertical
    )
  )
)

# Display the heatmap in RStudio
draw(heatmap, merge_legends = TRUE, annotation_legend_list = list(legend_annotations))

# Save the heatmap as a PDF
pdf(paste0(dataframeName, "_new_kmeans_", num_clusters, "_clusters_heatmap_", palette_name, ".pdf"), width = 10, height = 14)
draw(heatmap, merge_legends = TRUE, annotation_legend_list = list(legend_annotations))
dev.off()

# Create a summary table of gene counts per cluster
cluster_summary <- dt %>%
  group_by(NewCluster) %>%
  summarise(GeneCount = n()) %>%
  arrange(NewCluster)
# Print the summary table
print(cluster_summary)

# Create a sequential color palette
color_palette <- colorRampPalette(c("#FFFFFF", brewer_palette[9]))(100)

# Create a bar plot of gene counts per cluster with sequential color intensity
cluster_plot <- ggplot(cluster_summary, aes(x = factor(NewCluster), y = GeneCount)) +
  geom_bar(stat = "identity", 
           aes(fill = GeneCount), 
           width = 0.7,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = GeneCount), 
            vjust = -0.5, 
            size = 4,
            fontface = "bold") +
  scale_fill_gradientn(
    colors = color_palette,
    limits = c(0, max(cluster_summary$GeneCount)),
    breaks = pretty_breaks(n = 5)(c(0, max(cluster_summary$GeneCount))),
    guide = guide_colorbar(ticks = TRUE, nbin = 100, barwidth = 1, barheight = 10)
  ) +
  labs(title = "Number of Genes per Cluster",
       x = "Cluster",
       y = "Number of Genes",
       fill = "Gene Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
print(cluster_plot)
# Save the plot
#pdf(paste0(dataframeName, "_gene_counts_per_cluster_sequential_", palette_name, ".pdf"), width = 12, height = 8)

#dev.off()


# Optionally, save the summary table as a CSV file
#write.csv(cluster_summary, file = paste0(dataframeName, "_gene_counts_per_cluster.csv"), row.names = FALSE)

# Prepare and save output data (including NewCluster)
# Calculate row means for reordering
sample_columns <- grep("Sample", colnames(dt), value = TRUE) # Adjust the pattern to match your sample columns
dt[, RowMean := rowMeans(.SD, na.rm = TRUE), .SDcols = sample_columns]

# Order data based on NewCluster and RowMean
ordered_output_data <- dt[order(NewCluster, -abs(RowMean))]

# Remove RowMean column before saving
ordered_output_data[, RowMean := NULL]

# Save the output data as CSV and Excel files
fwrite(ordered_output_data, file = paste0(dataframeName, "_with_", num_clusters, "_new_clusters.csv"))
write_xlsx(as.data.frame(ordered_output_data), path = paste0(dataframeName, "_with_", num_clusters, "_new_clusters.xlsx"))

print("Output data has been saved as CSV and Excel files.")


# After adding NewCluster to the data, let's clean it up
data <- data %>%
  select(Gene.Name, NewCluster, everything()) %>%
  select(-Cluster) # Remove the original Cluster column and any NA column

# Convert data to data.table
dt <- as.data.table(data)

# Print column names to verify
print("Column names in the data:")
print(colnames(dt))

# Prepare data for plotting expression patterns
data_long <- melt(dt, id.vars = c("Gene.Name", "NewCluster"), 
                  variable.name = "TimePoint", value.name = "Expression")

# Remove any rows with NA in TimePoint (if any)
data_long <- data_long[!is.na(data_long$TimePoint),]

# Print unique values in TimePoint to verify
print("Unique values in TimePoint:")
print(unique(data_long$TimePoint))

data_long$CellType <- ifelse(grepl("^g", data_long$TimePoint), "Symbiotic", "Aposymbiotic")
data_long$CellType <- factor(data_long$CellType, levels = c("Symbiotic", "Aposymbiotic"))

# Define the correct time point order
time_order <- c("gAM8", "gAM11", "gPM2", "gPM5", "gPM8", "gPM11", 
                "wAM8", "wAM11", "wPM2", "wPM5", "wPM8", "wPM11")
data_long$TimePoint <- factor(data_long$TimePoint, levels = time_order)

# Calculate the mean expression per cluster at each time point for visualization
average_expression <- data_long %>%
  group_by(NewCluster, TimePoint, CellType) %>%
  summarise(Expression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Print the first few rows of average_expression to verify
print("First few rows of average_expression:")
print(head(average_expression))

# Plot Average Expression Pattern Trajectories by Cluster
p_avg <- ggplot(average_expression, aes(x = TimePoint, y = Expression, group = NewCluster, color = as.factor(NewCluster))) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~ CellType, scales = "free_x", ncol = 2) +
  theme_minimal() +
  labs(title = "Average Expression Pattern Trajectories by NewCluster",
       x = "Time Point", y = "Z-Scaled Expression", color = "NewCluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_color_brewer(palette = "Set1")

# Display the average expression plot
print(p_avg)

# Save the average expression plot
ggsave(paste0(dataframeName, "_new_kmeans_", num_clusters ,"_average_expression_patterns.pdf"), plot = p_avg, width = 12, height = 8)
print("Average expression plot saved.")

# Create a list to store plots for each cluster
plot_list <- list()

# Loop through each cluster to create expression pattern plots
for (i in sort(unique(data_long$NewCluster))) {
  cluster_data <- subset(data_long, NewCluster == i)
  
  p <- ggplot(cluster_data, aes(x = TimePoint, y = Expression, group = Gene.Name)) +
    geom_line(aes(color = Gene.Name), alpha = 0.2, linewidth = 0.2) +
    stat_summary(fun = mean, geom = "line", color = "blue", linewidth = 1, aes(group = 1)) +
    facet_wrap(~ CellType, scales = "free_x", ncol = 2) +
    theme_minimal() +
    labs(title = paste("NewCluster", i, "Expression Patterns"), x = "Time Point", y = "Z-Scaled Expression") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  plot_list[[i]] <- p
}

# Create the collage
collage <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3))

# Display the collage
#print(do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3)))

# Create and save the collage
collageFilename <- paste0(dataframeName, "_new_kmeans_", num_clusters, "_expression_patterns_collage.pdf")
total_height <- 6 * ceiling(length(plot_list) / 3)
pdf(collageFilename, width = 18, height = total_height)
do.call(gridExtra::grid.arrange, c(plot_list, ncol = 3))
dev.off()
print(paste("Expression patterns collage saved as", collageFilename))


print("All analyses and visualizations have been completed.")