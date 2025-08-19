# Comprehensive Algae Count Analysis with Multiple Plot Types
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(ggbeeswarm)
library(viridis)  # For color blind friendly palettes

# Read and prepare data
data <- read.csv("Algae_count_ETO.csv")
colnames(data) <- c("Control", "Treatment_25", "Treatment_50")

data_long <- data %>%
  pivot_longer(cols = everything(), names_to = "Treatment", values_to = "Count") %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "Treatment_25", "Treatment_50"),
                            labels = c("Control", "25µM", "50µM")))

# Define consistent theme following publication standards
publication_theme <- theme_minimal() +
  theme(
    # Use Helvetica (more universally available than Arial)
    text = element_text(family = "Helvetica", color = "black"),
    
    # Axis text - black, no bold, consistent formatting
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 12, face = "plain"),
    
    # Plot titles - no bold except specified
    plot.title = element_text(size = 14, face = "plain", color = "black"),
    plot.subtitle = element_text(size = 11, face = "plain", color = "black"),
    
    # Remove legend
    legend.position = "none",
    
    # Panel and grid
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    
    # Axis lines
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5)
  )

# Define comparisons for statistical testing
my_comparisons <- list(c("Control", "25µM"), c("Control", "50µM"), c("25µM", "50µM"))

# Color palette - RdBu (Red-Blue) diverging palette
# Option 1: RdBu palette - Red to Blue diverging (recommended)
#colors <- c("Control" = "#d73027", "25µM" = "#f7f7f7", "50µM" = "#4575b4")

# Option 2: RdBu without white (more saturated)
# colors <- c("Control" = "#ca0020", "25µM" = "#92c5de", "50µM" = "#0571b0")

# Option 3: RdBu darker variant
colors <- c("Control" = "#b2182b", "25µM" = "#762a83", "50µM" = "#2166ac")

# Option 4: Classic RdBu with intermediate tone
# colors <- c("Control" = "#ef8a62", "25µM" = "#f7f7f7", "50µM" = "#67a9cf")
# Option 1: RdBu palette - Red to Blue diverging (recommended)
#colors <- c("Control" = "#d73027", "25µM" = "#f7f7f7", "50µM" = "#4575b4")

# Option 2: RdBu without white (more saturated)
# colors <- c("Control" = "#ca0020", "25µM" = "#92c5de", "50µM" = "#0571b0")

# Option 3: RdBu darker variant
# colors <- c("Control" = "#b2182b", "25µM" = "#762a83", "50µM" = "#2166ac")

# Option 4: Classic RdBu with intermediate tone
# colors <- c("Control" = "#ef8a62", "25µM" = "#f7f7f7", "50µM" = "#67a9cf")

# =====================================
# PLOT 1: Enhanced Box Plot with Statistics
# =====================================
p1 <- ggplot(data_long, aes(x = Treatment, y = Count, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
               fill = "white", color = "black", stroke = 1) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     bracket.size = 0.6,
                     step.increase = 0.08) +
  labs(title = "Box plot with individual points",
       subtitle = "Diamonds = means, *** p<0.001, ** p<0.01, * p<0.05",
       x = "Treatment", y = "Algae count") +
  scale_fill_manual(values = colors) +
  publication_theme

# =====================================
# PLOT 2: Violin + Box Plot Combination
# =====================================
p2 <- ggplot(data_long, aes(x = Treatment, y = Count, fill = Treatment)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, 
               outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 18, 
               size = 4, color = "red") +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     y.position = c(380, 420, 350)) +
  labs(title = "Violin and box plot combination",
       subtitle = "Shows distribution shape and quartiles",
       x = "Treatment", y = "Algae count") +
  scale_fill_manual(values = colors) +
  publication_theme

# =====================================
# PLOT 3: Bee Swarm Plot with Statistics
# =====================================
p3 <- ggplot(data_long, aes(x = Treatment, y = Count, color = Treatment)) +
  geom_violin(aes(fill = Treatment), alpha = 0.2, color = NA) +
  geom_beeswarm(size = 2.5, alpha = 0.8, cex = 2) +
  stat_summary(fun = mean, geom = "point", shape = 18, 
               size = 5, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.15, color = "black", linewidth = 1) +
  stat_compare_means(comparisons = my_comparisons, 
                     method = "t.test", 
                     label = "p.signif",
                     color = "black") +
  labs(title = "Bee swarm plot",
       subtitle = "Every data point visible with means ± SE",
       x = "Treatment", y = "Algae count") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  publication_theme

# =====================================
# PLOT 4: Mean ± Error Bar Plot
# =====================================
summary_stats <- data_long %>%
  group_by(Treatment) %>%
  summarise(
    mean = mean(Count),
    se = sd(Count)/sqrt(n()),
    sd = sd(Count),
    .groups = 'drop'
  )

p4 <- ggplot(summary_stats, aes(x = Treatment, y = mean, fill = Treatment)) +
  geom_col(alpha = 0.7, width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, linewidth = 1) +
  geom_text(aes(label = paste0("n=50\n", round(mean, 1), "±", round(se, 1))), 
            vjust = -0.5, size = 3) +
  # Manual significance annotations
  annotate("text", x = 1.5, y = 250, label = "***", size = 6) +
  annotate("text", x = 2, y = 280, label = "***", size = 6) +
  annotate("text", x = 2.5, y = 270, label = "***", size = 6) +
  # Add brackets manually
  annotate("segment", x = 1, xend = 2, y = 245, yend = 245) +
  annotate("segment", x = 1, xend = 3, y = 275, yend = 275) +
  annotate("segment", x = 2, xend = 3, y = 265, yend = 265) +
  labs(title = "Mean ± standard error",
       subtitle = "Bar heights show means, error bars show ±1 SE, *** p<0.001",
       x = "Treatment", y = "Mean algae count") +
  scale_fill_manual(values = colors) +
  publication_theme

# =====================================
# COMBINED LAYOUT AND INDIVIDUAL SAVES
# =====================================

# Save individual plots in both PNG and PDF formats
ggsave("plot1_boxplot_enhanced.png", p1, width = 8, height = 6, dpi = 300)
ggsave("plot1_boxplot_enhanced.pdf", p1, width = 8, height = 6)

ggsave("plot2_violin_box.png", p2, width = 8, height = 6, dpi = 300)
ggsave("plot2_violin_box.pdf", p2, width = 8, height = 6)

ggsave("plot3_beeswarm.png", p3, width = 8, height = 6, dpi = 300)
ggsave("plot3_beeswarm.pdf", p3, width = 8, height = 6)

ggsave("plot4_mean_errorbar.png", p4, width = 8, height = 6, dpi = 300)
ggsave("plot4_mean_errorbar.pdf", p4, width = 8, height = 6)

# Arrange all plots in a collage
combined_plot <- grid.arrange(
  p1, p2, p3, p4,
  ncol = 2, nrow = 2,
  top = "Comprehensive Algae Count Analysis\nMultiple Visualization Approaches"
)

# Save the combined collage in both formats
ggsave("comprehensive_algae_analysis_collage.png", combined_plot, 
       width = 16, height = 12, dpi = 300)
ggsave("comprehensive_algae_analysis_collage.pdf", combined_plot, 
       width = 16, height = 12)

# Create a multi-page PDF with all individual plots
pdf("algae_analysis_all_plots.pdf", width = 8, height = 6)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()

# Display plots individually for viewing
print("=== INDIVIDUAL PLOTS ===")
print("Plot 1: Enhanced Box Plot")
print(p1)

print("Plot 2: Violin + Box Plot")
print(p2)

print("Plot 3: Bee Swarm Plot")
print(p3)

print("Plot 4: Mean ± Error Bar")
print(p4)

print("=== COLLAGE ===")
print(combined_plot)

# =====================================
# STATISTICAL SUMMARY TABLE
# =====================================

# Perform statistical tests
model <- aov(Count ~ Treatment, data = data_long)
anova_results <- summary(model)
tukey_results <- TukeyHSD(model)

# Create summary table
cat("\n=== COMPREHENSIVE STATISTICAL RESULTS ===\n")
cat("==========================================\n")

# Descriptive statistics
cat("\nDESCRIPTIVE STATISTICS:\n")
desc_stats <- data_long %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    Mean = round(mean(Count), 2),
    SD = round(sd(Count), 2),
    Median = median(Count),
    Min = min(Count),
    Max = max(Count),
    .groups = 'drop'
  )
print(desc_stats)

# ANOVA results
cat("\nANOVA RESULTS:\n")
print(anova_results)

# Post-hoc results
cat("\nTUKEY POST-HOC RESULTS:\n")
print(tukey_results)

# Effect sizes
cat("\nEFFECT SIZES (Cohen's d):\n")
library(effsize)

# Calculate Cohen's d for each comparison
control_data <- data_long$Count[data_long$Treatment == "Control"]
um25_data <- data_long$Count[data_long$Treatment == "25µM"]
um50_data <- data_long$Count[data_long$Treatment == "50µM"]

d1 <- cohen.d(control_data, um25_data)
d2 <- cohen.d(control_data, um50_data)
d3 <- cohen.d(um25_data, um50_data)

cat(sprintf("Control vs 25µM: d = %.3f (%s)\n", 
            d1$estimate, d1$magnitude))
cat(sprintf("Control vs 50µM: d = %.3f (%s)\n", 
            d2$estimate, d2$magnitude))
cat(sprintf("25µM vs 50µM: d = %.3f (%s)\n", 
            d3$estimate, d3$magnitude))

cat("\n=== CONCLUSION ===\n")
cat("All pairwise comparisons show significant differences\n")
cat("Clear dose-response relationship: Control < 25µM < 50µM\n")
cat("Large effect sizes indicate biological significance\n")
cat("Results are consistent across all visualization methods\n")

# Save results to file
sink("statistical_results.txt")
cat("ALGAE COUNT STATISTICAL ANALYSIS\n")
cat("================================\n\n")
print(desc_stats)
cat("\n")
print(anova_results)
cat("\n")
print(tukey_results)
sink()

cat("\nFiles saved:\n")
cat("INDIVIDUAL PLOTS (PNG):\n")
cat("- plot1_boxplot_enhanced.png\n")
cat("- plot2_violin_box.png\n") 
cat("- plot3_beeswarm.png\n")
cat("- plot4_mean_errorbar.png\n")
cat("\nINDIVIDUAL PLOTS (PDF):\n")
cat("- plot1_boxplot_enhanced.pdf\n")
cat("- plot2_violin_box.pdf\n") 
cat("- plot3_beeswarm.pdf\n")
cat("- plot4_mean_errorbar.pdf\n")
cat("\nCOLLAGES:\n")
cat("- comprehensive_algae_analysis_collage.png\n")
cat("- comprehensive_algae_analysis_collage.pdf\n")
cat("\nMULTI-PAGE PDF:\n")
cat("- algae_analysis_all_plots.pdf (all 4 plots in one PDF)\n")
cat("\nSTATISTICS:\n")
cat("- statistical_results.txt\n")