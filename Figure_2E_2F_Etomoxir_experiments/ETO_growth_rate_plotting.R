# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggsignif)
library(ggrepel)
library(extrafont)
library(scales)

# If extrafont is not installed, uncomment and run the following line:
# install.packages("extrafont")

# Import fonts (run this once, then you can comment it out)
# font_import()
# loadfonts(device = "win")  # For Windows
# loadfonts(device = "mac")  # For Mac

# Set the default font to Times New Roman
theme_update(text = element_text(family = "Times New Roman"))

# ============================
# Configuration: Define Prefix and Filenames
# ============================

# Define the prefix for all output files
# It's recommended to avoid spaces in filenames. Replace spaces with underscores or hyphens.
output_prefix <- "ETO_Growth"  # Change this to your desired prefix without spaces

# Optional: Define output directory to organize your plots
output_dir <- "plots"  # Change or set to NULL if you don't want a separate directory
if (!is.null(output_dir)) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
}

# Function to generate filenames with the prefix
generate_filename <- function(plot_type, extension = "pdf") {
  if (!is.null(output_dir)) {
    paste0(output_dir, "/", output_prefix, "_", plot_type, ".", extension)
  } else {
    paste0(output_prefix, "_", plot_type, ".", extension)
  }
}

# Customization options (modify these as needed)
FONT_SIZE <- 12          # Base font size for the plot
LINE_WIDTH <- 1          # Width of lines in the plots
POINT_SIZE <- 2          # Size of points in scatter plots
ERROR_BAR_WIDTH <- 0.1   # Width of error bars (horizontal extent)
ERROR_BAR_THICKNESS <- 1 # Thickness of error bars (line thickness)
LEGEND_POSITION <- "bottom"  # Position of the legend ("top", "bottom", "left", "right", or c(x,y) coordinates)
PLOT_TITLE_SIZE <- 16    # Font size for the main plot title
AXIS_TITLE_SIZE <- 14    # Font size for axis titles
PLOT_WIDTH <- 6          # Width of saved plots in inches
PLOT_HEIGHT <- 4         # Height of saved plots in inches

# Function to create a base theme with Times New Roman
theme_times_new_roman <- function(base_size = FONT_SIZE, base_family = "Times New Roman") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      text = element_text(family = base_family),
      plot.title = element_text(size = PLOT_TITLE_SIZE, face = "bold"),
      axis.title = element_text(size = AXIS_TITLE_SIZE),
      legend.position = LEGEND_POSITION
    )
}

# ============================
# Data Preparation
# ============================

# Read the data
data <- read.csv("ETO_growth.csv")  # Ensure this path is correct

# Get the original order of Cell Types
original_cell_types <- names(data)[-1] %>% 
  sub("_R[123]$", "", .) %>% 
  unique()

# Reshape the data, keeping individual replicates
data_long <- data %>%
  pivot_longer(cols = -Day, 
               names_to = c("Cell_Type", "Replicate"), 
               names_pattern = "(.+)_(R[123])",
               values_to = "Cell_Count") %>%
  mutate(Day = factor(Day, levels = unique(Day)),
         Cell_Type = factor(Cell_Type, levels = original_cell_types))

# Calculate mean and standard error
data_summary <- data_long %>%
  group_by(Day, Cell_Type) %>%
  summarise(Mean_Count = mean(Cell_Count),
            SD = sd(Cell_Count),
            SE = SD / sqrt(n()),
            .groups = 'drop')

# Define a color palette with distinct colors for each cell type
n_colors <- length(original_cell_types)
color_palette <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(n_colors), original_cell_types)

# Check if the y-axis should be log-scaled
y_is_log <- all(data_summary$Mean_Count > 0) && 
  (max(data_summary$Mean_Count) / min(data_summary$Mean_Count) > 100)

# Determine y-axis label
y_label <- if(y_is_log) {
  "Average Cell Count (log scale)"
} else {
  "Average Cell Count (linear scale)"
}

# Get the last day for labeling
last_day_data <- data_summary %>%
  group_by(Cell_Type) %>%
  slice(which.max(as.integer(Day)))

# ============================
# Plot Generation
# ============================

# 1. Line plot of growth curves with end labels
line_plot <- ggplot(data_summary, aes(x = Day, y = Mean_Count, color = Cell_Type, group = Cell_Type)) +
  geom_line(size = LINE_WIDTH) +
  geom_point(size = POINT_SIZE) +
  geom_errorbar(aes(ymin = Mean_Count - SE, ymax = Mean_Count + SE), width = 0.2) +
  geom_text_repel(
    data = last_day_data,
    aes(label = Cell_Type),
    nudge_x = 0.5,
    hjust = 0,
    segment.color = NA,
    show.legend = FALSE,
    family = "Times New Roman"
  ) +
  scale_color_manual(values = color_palette) +
  labs(title = "Cell Growth Over Time",
       x = "Days",
       y = "Average Cell Count") +
  theme_times_new_roman() +
  theme(legend.position = "none") +
  scale_x_discrete(expand = expansion(add = c(0, 1.5)))

print(line_plot)

# ============================
# Saving Plots as PDF Using the Prefix
# ============================

ggsave(generate_filename("line_plot"), line_plot, width = PLOT_WIDTH, height = PLOT_HEIGHT, device = cairo_pdf)

# To modify saved plot dimensions:
# - Change PLOT_WIDTH and PLOT_HEIGHT at the top of the script
# - Or, specify width and height directly in ggsave() function calls
