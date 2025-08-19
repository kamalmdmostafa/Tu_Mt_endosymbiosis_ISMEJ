# AAI Matrix Heatmap Script
library(ggplot2)
library(reshape2)
library(dplyr)

# Create the AAI matrix data
aai_matrix <- matrix(c(
  100, 89.3, 89.0, 82.0, 70.0,
  89.3, 100, 89.5, 82.7, 70.6,
  89.0, 89.5, 100, 82.8, 70.7,
  82.0, 82.7, 82.8, 100, 69.7,
  70.0, 70.6, 70.7, 69.7, 100
), nrow = 5, byrow = TRUE)

# Set row and column names in the correct order from the figure
species_names <- c("T. utriculariae", "T. malaccensis", "T. thermophila", "T. elliotti", "T. borealis")
rownames(aai_matrix) <- species_names
colnames(aai_matrix) <- species_names

# Convert matrix to long format for ggplot
aai_long <- melt(aai_matrix, varnames = c("Species1", "Species2"), value.name = "AAI")

# Create species factor with proper order - REVERSED for y-axis to match figure
aai_long$Species1 <- factor(aai_long$Species1, levels = rev(species_names))
aai_long$Species2 <- factor(aai_long$Species2, levels = species_names)

# Round AAI values for display
aai_long$AAI_display <- round(aai_long$AAI, 0)

# Create color scheme similar to the image (blue-yellow gradient)
color_palette <- c("#FFEB3B", "#FFC107", "#FF9800", "#9E9E9E", "#607D8B", "#3F51B5", "#1A237E")

# Create the heatmap
aai_plot <- ggplot(aai_long, aes(x = Species2, y = Species1, fill = AAI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = AAI_display), color = "white", size = 4, fontface = "bold") +
  scale_fill_gradientn(
    colors = color_palette,
    values = scales::rescale(c(70, 75, 80, 85, 90, 95, 100)),
    limits = c(70, 100),
    name = "Percentage",
    breaks = seq(70, 100, 5),
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5),
      barwidth = 1,
      barheight = 8
    )
  ) +
  labs(
    title = "D",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    # Text formatting
    text = element_text(color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0, vjust = 1),
    
    # Axis text formatting
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, 
                               color = c("#4CAF50", "#E91E63", "#2196F3", "black", "black")),
    axis.text.y = element_text(size = 11, hjust = 1,
                               color = c("black", "black", "#2196F3", "#E91E63", "#4CAF50")),
    
    # Remove panel elements
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # Legend formatting
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    
    # Axis ticks
    axis.ticks = element_blank(),
    
    # Plot margins
    plot.margin = margin(10, 40, 10, 10)
  ) +
  coord_fixed(ratio = 1)

# Format species names with italics (approximated with color coding)
# Green for T. utriculariae, Pink for T. malaccensis, Blue for T. thermophila

print(aai_plot)

# Save the plot
ggsave("aai_matrix_heatmap.png", aai_plot, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("aai_matrix_heatmap.pdf", aai_plot, width = 8, height = 6, bg = "white")

# Alternative version with better text formatting using expression
aai_plot_formatted <- ggplot(aai_long, aes(x = Species2, y = Species1, fill = AAI)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = AAI_display), color = "white", size = 4, fontface = "bold") +
  scale_fill_gradientn(
    colors = color_palette,
    values = scales::rescale(c(70, 75, 80, 85, 90, 95, 100)),
    limits = c(70, 100),
    name = "Percentage",
    breaks = seq(70, 100, 5),
    guide = guide_colorbar(
      title.position = "right",
      title.theme = element_text(angle = 90, hjust = 0.5),
      barwidth = 1,
      barheight = 8
    )
  ) +
  scale_x_discrete(
    labels = c(
      expression(italic("T. utriculariae")),
      expression(italic("T. malaccensis")),
      expression(italic("T. thermophila")),
      expression(italic("T. elliotti")),
      expression(italic("T. borealis"))
    )
  ) +
  scale_y_discrete(
    labels = c(
      expression(italic("T. borealis")),
      expression(italic("T. elliotti")),
      expression(italic("T. thermophila")),
      expression(italic("T. malaccensis")),
      expression(italic("T. utriculariae"))
    )
  ) +
  labs(
    title = "D",
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    # Text formatting
    text = element_text(color = "black"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0, vjust = 1),
    
    # Axis text formatting with colors
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
    axis.text.y = element_text(size = 11, hjust = 1),
    
    # Remove panel elements
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    # Legend formatting
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    
    # Axis ticks
    axis.ticks = element_blank(),
    
    # Plot margins
    plot.margin = margin(10, 40, 10, 10)
  ) +
  coord_fixed(ratio = 1)

print(aai_plot_formatted)

# Save the formatted version
ggsave("aai_matrix_heatmap_formatted.png", aai_plot_formatted, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("aai_matrix_heatmap_formatted.pdf", aai_plot_formatted, width = 8, height = 6, bg = "white")

cat("AAI Matrix Heatmap Complete!\n")
cat("Files saved:\n")
cat("- aai_matrix_heatmap.png (with colored species names)\n")
cat("- aai_matrix_heatmap.pdf\n")
cat("- aai_matrix_heatmap_formatted.png (with italic species names)\n")
cat("- aai_matrix_heatmap_formatted.pdf\n")