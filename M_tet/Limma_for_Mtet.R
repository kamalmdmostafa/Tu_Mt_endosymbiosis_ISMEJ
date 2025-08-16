# Load necessary libraries
library(limma)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(viridis)
library(gridExtra)
library(grid)
library(ggrepel)
library(RColorBrewer)
library(circlize)

# Create directory structure
dir.create("results", showWarnings = FALSE)
dir.create("results/plots", showWarnings = FALSE)
dir.create("results/plots/individual", showWarnings = FALSE)

# Set plot parameters
PLOT_PARAMS <- list(
  point_size = 2,              # Size of points in scatter plots
  point_alpha = 0.6,           # Transparency of points
  title_size = 14,             # Size of plot titles
  axis_title_size = 12,        # Size of axis titles
  axis_text_size = 10,         # Size of axis text
  legend_text_size = 10,       # Size of legend text
  annotation_text_size = 3,    # Size of annotation text
  pdf_width = 10,             # Width of PDF output (increased for annotation space)
  pdf_height = 8,             # Height of PDF output (increased for annotation space)
  heatmap_height = 12,        # Height for heatmap
  heatmap_width = 10,         # Width for heatmap
  collage_width = 20,         # Width for collage plot
  collage_height = 24,        # Height for collage plot
  colors = list(
    upregulated = "#D55E00",
    downregulated = "#0072B2",
    nonsignificant = "#999999",
    conditions = c("Mt_Apo" = "#E69F00", "Mt_Sym" = "#56B4E9")
  ),
  # Add annotation panel settings
  annotation_panel = list(
    width = 0.3,              # Width of annotation panel relative to plot
    background_alpha = 0.9,   # Background transparency
    border_color = "grey70"   # Border color for annotation box
  )
)

# Analysis parameters
ANALYSIS_PARAMS <- list(
  tpm_threshold = 0.005,      # Minimum TPM value
  sample_percentage = 0.75,    # Minimum percentage of samples passing threshold
  pvalue_threshold = 0.05,    # Adjusted p-value threshold (changed to match your output)
  logfc_threshold = 1         # Log2 fold change threshold
)

# Function to preprocess data
preprocess_data <- function(raw_data, tpm_threshold, sample_percentage) {
  # Convert tibble to data frame
  raw_data <- as.data.frame(raw_data)
  
  # Ensure the first column is gene_id
  if (!all(colnames(raw_data)[1] == "gene_id")) {
    stop("First column of raw_data must be 'gene_id'")
  }
  
  # Set gene_id as row names and remove the column
  rownames(raw_data) <- raw_data$gene_id
  raw_data <- raw_data[,-1]
  
  # Filter genes
  keep <- rowSums(raw_data >= tpm_threshold) >= (ncol(raw_data) * sample_percentage)
  filtered_data <- raw_data[keep, ]
  
  # Log2 transform
  log_data <- log2(filtered_data + 1)
  
  return(log_data)
}

# Read TPM data
tpm_data_file <- "kallisto_tpm_matrix.tsv"  # Replace with your TPM data file path

# Check if file exists
if (!file.exists(tpm_data_file)) {
  stop("TPM data file not found: ", tpm_data_file)
}

# Read TPM data
tpm_data <- try(read.delim(tpm_data_file, sep="\t", header=TRUE))
if (inherits(tpm_data, "try-error")) {
  stop("Error reading TPM data file: ", tpm_data_file)
}

# Create and prepare metadata
metadata <- data.frame(
  Sample_name = c("Endo_M_2", "Endo_M_4", "Endo_M_5", 
                  "Exo_M_1", "Exo_M_2", "Exo_M_5"),
  Condition = factor(c(rep("Mt_Sym", 3), rep("Mt_Apo", 3)), 
                     levels = c("Mt_Apo", "Mt_Sym"))
)

# Preprocess data using analysis parameters
tpm_data <- preprocess_data(tpm_data, 
                            tpm_threshold = ANALYSIS_PARAMS$tpm_threshold, 
                            sample_percentage = ANALYSIS_PARAMS$sample_percentage)

# Ensure the columns of tpm_data match the order of samples in metadata
tpm_data <- tpm_data[, metadata$Sample_name]

# Create design matrix
design <- model.matrix(~Condition, data = metadata)
colnames(design) <- c("Intercept", "Mt_Sym_vs_Apo")

# Fit model using limma
fit <- lmFit(tpm_data, design)
fit <- eBayes(fit)

# Extract differentially expressed genes
results <- topTable(fit, coef="Mt_Sym_vs_Apo", adjust.method="BH", n=Inf, sort.by="none")
results$gene_id <- rownames(results)
results$Category <- ifelse(results$adj.P.Val <= ANALYSIS_PARAMS$pvalue_threshold & 
                             results$logFC >= ANALYSIS_PARAMS$logfc_threshold, "Upregulated",
                           ifelse(results$adj.P.Val <= ANALYSIS_PARAMS$pvalue_threshold & 
                                    results$logFC <= -ANALYSIS_PARAMS$logfc_threshold, "Downregulated", 
                                  "None"))

# Write results to file
write.csv(results, "results/Mt_Sym_vs_Apo_DE_genes.csv", row.names = FALSE)

# Calculate DE gene counts and create stats label (used across multiple plots)
up_genes <- sum(results$Category == "Upregulated")
down_genes <- sum(results$Category == "Downregulated")
stats_label <- sprintf("Up-regulated: %d\nDown-regulated: %d\nFDR = %.2f, |log2FC| = %.1f", 
                       up_genes, down_genes, 
                       ANALYSIS_PARAMS$pvalue_threshold, 
                       ANALYSIS_PARAMS$logfc_threshold)

# Create volcano plot
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
  geom_point(alpha = PLOT_PARAMS$point_alpha, size = PLOT_PARAMS$point_size) +
  scale_color_manual(values = c("Upregulated" = PLOT_PARAMS$colors$upregulated, 
                                "Downregulated" = PLOT_PARAMS$colors$downregulated, 
                                "None" = PLOT_PARAMS$colors$nonsignificant)) +
  ggtitle("Mt_Sym vs Mt_Apo Differential Expression") +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "grey70", size = 0.2),
    legend.title = element_text(size = PLOT_PARAMS$axis_title_size),
    legend.text = element_text(size = PLOT_PARAMS$legend_text_size),
    plot.title = element_text(size = PLOT_PARAMS$title_size, face = "bold"),
    axis.title = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.text = element_text(size = PLOT_PARAMS$axis_text_size),
    plot.margin = margin(1, 2, 1, 1, "cm")  # Increased right margin for annotation
  ) +
  geom_hline(yintercept = -log10(ANALYSIS_PARAMS$pvalue_threshold), 
             linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = c(-ANALYSIS_PARAMS$logfc_threshold, ANALYSIS_PARAMS$logfc_threshold), 
             linetype = "dashed", color = "grey50") +
  annotate("text", 
           x = min(results$logFC) * 0.8, 
           y = max(-log10(results$adj.P.Val)) * 0.95,
           label = stats_label,
           hjust = 0,
           size = PLOT_PARAMS$annotation_text_size,
           color = "black",
           bg = "white") +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted P-value")

# Display volcano plot
print(volcano_plot)

# Save volcano plot
pdf(file.path("results/plots/individual", "volcano_plot.pdf"), 
    width = PLOT_PARAMS$pdf_width, height = PLOT_PARAMS$pdf_height)
print(volcano_plot)
dev.off()

# Create MA plot
ma_plot <- ggplot(results, aes(x = AveExpr, y = logFC, color = Category)) +
  geom_point(alpha = PLOT_PARAMS$point_alpha, size = PLOT_PARAMS$point_size) +
  scale_color_manual(values = c("Upregulated" = PLOT_PARAMS$colors$upregulated, 
                                "Downregulated" = PLOT_PARAMS$colors$downregulated, 
                                "None" = PLOT_PARAMS$colors$nonsignificant)) +
  ggtitle("MA Plot: Mt_Sym vs Mt_Apo") +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.75),
    legend.background = element_rect(fill = "white", color = "grey70", size = 0.2),
    legend.title = element_text(size = PLOT_PARAMS$axis_title_size),
    legend.text = element_text(size = PLOT_PARAMS$legend_text_size),
    plot.title = element_text(size = PLOT_PARAMS$title_size, face = "bold"),
    axis.title = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.text = element_text(size = PLOT_PARAMS$axis_text_size),
    plot.margin = margin(1, 2, 1, 1, "cm")  # Increased right margin for annotation
  ) +
  geom_hline(yintercept = c(-ANALYSIS_PARAMS$logfc_threshold, 0, ANALYSIS_PARAMS$logfc_threshold), 
             linetype = "dashed", color = "grey50") +
  annotate("text", 
           x = min(results$AveExpr) * 1.2, 
           y = max(results$logFC) * 0.95,
           label = stats_label,
           hjust = 0,
           size = PLOT_PARAMS$annotation_text_size,
           color = "black",
           bg = "white") +
  xlab("Average Expression") +
  ylab("log2 Fold Change")

# Display MA plot
print(ma_plot)

# Save MA plot
pdf(file.path("results/plots/individual", "ma_plot.pdf"), 
    width = PLOT_PARAMS$pdf_width, height = PLOT_PARAMS$pdf_height)
print(ma_plot)
dev.off()

# Create summary plot
summary_plot <- ggplot(filter(results %>% group_by(Category) %>% 
                                summarise(Count = n()), Category != "None"), 
                       aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Upregulated" = PLOT_PARAMS$colors$upregulated, 
                               "Downregulated" = PLOT_PARAMS$colors$downregulated)) +
  labs(
    title = "Significant DE Genes",
    subtitle = sprintf("FDR = %.2f, |log2FC| = %.1f", 
                       ANALYSIS_PARAMS$pvalue_threshold, 
                       ANALYSIS_PARAMS$logfc_threshold)
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = PLOT_PARAMS$title_size, face = "bold"),
    plot.subtitle = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.title = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.text = element_text(size = PLOT_PARAMS$axis_text_size)
  ) +
  geom_text(aes(label = Count), 
            vjust = -0.5, 
            position = position_dodge(width = 0.9), 
            size = PLOT_PARAMS$annotation_text_size, 
            fontface = "bold") +
  xlab("Category") + 
  ylab("Count")

# Display summary plot
print(summary_plot)

# Save summary plot
pdf(file.path("results/plots/individual", "summary_plot.pdf"), 
    width = PLOT_PARAMS$pdf_width, height = PLOT_PARAMS$pdf_height)
print(summary_plot)
dev.off()

# Perform PCA and create plot
pca_result <- prcomp(t(tpm_data), scale. = TRUE)
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

pca_plot <- as.data.frame(pca_result$x) %>%
  mutate(Sample = rownames(.)) %>%
  merge(metadata, by.x = "Sample", by.y = "Sample_name") %>%
  ggplot(aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = PLOT_PARAMS$point_size) +
  geom_text_repel(size = PLOT_PARAMS$annotation_text_size, 
                  show.legend = FALSE, 
                  max.overlaps = Inf,
                  box.padding = 0.75, 
                  point.padding = 0.75, 
                  segment.size = 0.3) +
  scale_color_manual(values = PLOT_PARAMS$colors$conditions) +
  theme_minimal() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "grey70", size = 0.2),
    plot.title = element_text(size = PLOT_PARAMS$title_size, face = "bold"),
    plot.subtitle = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.title = element_text(size = PLOT_PARAMS$axis_title_size),
    axis.text = element_text(size = PLOT_PARAMS$axis_text_size)
  ) +
  labs(
    title = "PCA of Gene Expression Data",
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  )

# Display PCA plot
print(pca_plot)

# Save PCA plot
pdf(file.path("results/plots/individual", "pca_plot.pdf"), 
    width = PLOT_PARAMS$pdf_width, height = PLOT_PARAMS$pdf_height)
print(pca_plot)
dev.off()

# Prepare heatmap data (using all significant DEGs)
significant_genes <- results %>%
  filter(adj.P.Val <= ANALYSIS_PARAMS$pvalue_threshold & 
           abs(logFC) >= ANALYSIS_PARAMS$logfc_threshold) %>%
  arrange(desc(abs(logFC))) %>%
  pull(gene_id)

heatmap_data <- tpm_data[significant_genes, ]

# Create heatmap annotation
column_ha <- HeatmapAnnotation(
  Condition = metadata$Condition,
  col = list(Condition = PLOT_PARAMS$colors$conditions),
  annotation_legend_param = list(
    Condition = list(
      title = "Condition",
      at = c("Mt_Apo", "Mt_Sym"),
      labels = c("Mt_Apo", "Mt_Sym")
    )
  )
)

# Create heatmap
ht <- Heatmap(t(scale(t(heatmap_data))),
              name = "Z-score",
              col = viridis(100),
              show_row_names = FALSE,  # Hide gene names
              cluster_rows = TRUE,
              cluster_columns = FALSE,
              top_annotation = column_ha,
              column_title = paste0("All Significant DE Genes (n=", length(significant_genes), ")"),
              column_title_side = "top",
              show_column_names = TRUE,
              clustering_method_rows = "ward.D2",
              column_title_gp = gpar(fontsize = PLOT_PARAMS$title_size),
              heatmap_legend_param = list(
                title = "Z-score",
                legend_height = unit(4, "cm"),
                title_position = "lefttop-rot"
              ))

# Display heatmap
draw(ht)

# Save heatmap
pdf(file.path("results/plots/individual", "heatmap.pdf"), 
    width = PLOT_PARAMS$heatmap_width, height = PLOT_PARAMS$heatmap_height)
draw(ht)
dev.off()

# Create and save plot collage
pdf(file.path("results/plots", "DE_analysis_collage.pdf"), 
    width = PLOT_PARAMS$collage_width, 
    height = PLOT_PARAMS$collage_height)

# Arrange plots in a grid
grid.arrange(
  volcano_plot, ma_plot,
  summary_plot, pca_plot,
  grid.grabExpr(draw(ht)),
  layout_matrix = rbind(
    c(1,2),
    c(3,4),
    c(5,5)
  ),
  heights = c(1, 1, 1.2)
)
dev.off()

# Save session info and analysis summary
writeLines(capture.output(sessionInfo()), file.path("results", "session_info.txt"))

analysis_summary <- c(
  "Analysis Summary:",
  "----------------",
  sprintf("Total genes analyzed: %d", nrow(results)),
  sprintf("Up-regulated genes: %d", up_genes),
  sprintf("Down-regulated genes: %d", down_genes),
  sprintf("FDR threshold: %.3f", ANALYSIS_PARAMS$pvalue_threshold),
  sprintf("Log2FC threshold: %.2f", ANALYSIS_PARAMS$logfc_threshold),
  sprintf("TPM threshold: %.3f", ANALYSIS_PARAMS$tpm_threshold),
  sprintf("Sample percentage threshold: %.2f", ANALYSIS_PARAMS$sample_percentage)
)

writeLines(analysis_summary, file.path("results", "analysis_summary.txt"))
cat("\n", paste(analysis_summary, collapse = "\n"), "\n")