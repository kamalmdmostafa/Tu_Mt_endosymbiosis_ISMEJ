# Load necessary libraries
library(limma)

# Load your data
# Assuming your log-transformed TPM data is in 'processed_data.csv'
tpm_data <- processed_data
tpm_data <- as.data.frame(tpm_data)
metadata <- updated_metadata
# Set the first column (Gene.Name) as row names and remove it from the data frame
rownames(tpm_data) <- tpm_data$Gene.Name
tpm_data <- tpm_data[,-1] # Remove the Gene.Name column

# Ensure the columns of tpm_data match the order of samples in metadata
tpm_data <- tpm_data[, metadata$Sample_name]

# Prepare factors in metadata for the design matrix
metadata$Time <- factor(metadata$Time)
metadata$Cell_type <- factor(metadata$Cell_type)
metadata$Interaction <- interaction(metadata$Cell_type, metadata$Time)

# Design matrix
design <- model.matrix(~0 + metadata$Interaction)
colnames(design) <- make.names(colnames(design))


# Fit model using limma
fit <- lmFit(tpm_data, design)

# Apply empirical Bayes smoothing
fit <- eBayes(fit, trend=TRUE)


# Correcting the contrast matrix creation
# Expanding the contrast matrix to include all time points
contrast.matrix <- makeContrasts(
  # Comparisons between cell types at each time point
  g_vs_w_0h = metadata.Interactiong.0h - metadata.Interactionw.0h,
  g_vs_w_3h = metadata.Interactiong.3h - metadata.Interactionw.3h,
  g_vs_w_6h = metadata.Interactiong.6h - metadata.Interactionw.6h,
  g_vs_w_9h = metadata.Interactiong.9h - metadata.Interactionw.9h,
  g_vs_w_12h = metadata.Interactiong.12h - metadata.Interactionw.12h,
  g_vs_w_15h = metadata.Interactiong.15h - metadata.Interactionw.15h,
  # Assuming transitions were also of interest; keep these if still relevant
  g_0h_to_3h = metadata.Interactiong.0h - metadata.Interactiong.3h,
  g_12h_to_15h = metadata.Interactiong.12h - metadata.Interactiong.15h,
  w_0h_to_3h = metadata.Interactionw.0h - metadata.Interactionw.3h,
  w_12h_to_15h = metadata.Interactionw.12h - metadata.Interactionw.15h,
  levels = design
)

# Apply the contrasts to the fit
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract differentially expressed genes for each contrast
# Loop through each coefficient in the contrast matrix
for(i in 1:ncol(fit2$coefficients)) {
  coef_name <- colnames(fit2$coefficients)[i]
  de_genes <- topTable(fit2, coef=i, adjust.method="BH", n=Inf)
  file_name <- paste("DE_genes_", coef_name, ".csv", sep="")
  write.csv(de_genes, file_name)
}

# Initialize an empty data frame to store combined DE genes from all contrasts
combined_de_genes <- data.frame()

# Initialize an empty data frame to store combined DE genes from all contrasts
combined_de_genes <- data.frame()

# Loop through each coefficient in the contrast matrix to extract DE genes
for(i in 1:ncol(fit2$coefficients)) {
  coef_name <- colnames(fit2$coefficients)[i]
  de_genes <- topTable(fit2, coef=i, adjust.method="BH", n=Inf, sort.by="none")
  
  # Check if de_genes is not empty
  if (nrow(de_genes) > 0) {
    # Convert row names (Gene Names) into a column
    de_genes$Gene.Name <- rownames(de_genes)
    
    # Add a column to indicate the contrast
    de_genes$Contrast <- coef_name
    
    # Categorize genes based on log2FC and FDR
    de_genes$Category <- ifelse(de_genes$adj.P.Val <= 0.01 & de_genes$logFC >= 0.585, "Upregulated",
                                ifelse(de_genes$adj.P.Val <= 0.01 & de_genes$logFC <= -0.585, "Downregulated", 
                                       "None"))
    
    # Combine with the main data frame
    combined_de_genes <- rbind(combined_de_genes, de_genes)
  }
}

# After loop check to ensure combined_de_genes is not empty
if (nrow(combined_de_genes) == 0) {
  cat("No differentially expressed genes found under the specified criteria.\n")
} else {
  # Write the combined and categorized DE gene list to a CSV file
  write.csv(combined_de_genes, "combined_categorized_DE_genes.csv", row.names = FALSE)
}

#####################extracting_significant_gene_expression_data#######################
library(dplyr)
# Specify the contrasts of interest
time_point_contrasts <- c("g_vs_w_0h", "g_vs_w_3h", "g_vs_w_6h", "g_vs_w_9h", "g_vs_w_12h", "g_vs_w_15h")

# Filter for significant upregulated or downregulated genes within these contrasts
significant_genes <- combined_de_genes %>%
  filter(Contrast %in% time_point_contrasts & (Category == "Upregulated" | Category == "Downregulated"))

# Get the names of significant genes
significant_gene_names <- significant_genes$Gene.Name

# Extract expression data for significant genes
significant_gene_data <- tpm_data[rownames(tpm_data) %in% significant_gene_names, ]

write.csv(significant_gene_data, "TPM_of_Sig_DEG.csv", row.names = TRUE)



#DE statistics
# Calculate the count of each category within each contrast
de_stats <- aggregate(. ~ Contrast + Category, data = combined_de_genes, FUN = length)

# The aggregate function might give the count column a name based on the column it counted. Let's rename it to 'Count'.
# This step might need adjustment based on your actual column names.
names(de_stats)[names(de_stats) == "Gene.Name"] <- "Count"

library(ggplot2)

# Define the contrast orders
time_point_contrasts <- c("g_vs_w_0h", "g_vs_w_3h", "g_vs_w_6h", "g_vs_w_9h", "g_vs_w_12h", "g_vs_w_15h")
transition_contrasts <- c("g_0h_to_3h", "g_12h_to_15h", "w_0h_to_3h", "w_12h_to_15h")

# Ensure all categories are present for coloring
de_stats$Category <- factor(de_stats$Category, levels = c("Upregulated", "Downregulated", "None"))

# Filter and prepare the dataset for each plot
time_point_data_sig <- subset(de_stats, Contrast %in% time_point_contrasts & Category != "None")
time_point_data_none <- subset(de_stats, Contrast %in% time_point_contrasts & Category == "None")
transition_data_sig <- subset(de_stats, Contrast %in% transition_contrasts & Category != "None")
transition_data_none <- subset(de_stats, Contrast %in% transition_contrasts & Category == "None")

#Plot 1: Significant Time-point Comparisons
plot1 <- ggplot(time_point_data_sig, aes(x = Contrast, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "darkgreen")) +
  scale_x_discrete(limits = time_point_contrasts) +
  labs(title = "Significant DE Genes Over Time Points", x = "Contrast", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot 2: Non-significant Time-point Comparisons
plot2 <- ggplot(time_point_data_none, aes(x = Contrast, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("None" = "grey")) +
  scale_x_discrete(limits = time_point_contrasts) +
  labs(title = "Non-significant DE Genes Over Time Points", x = "Contrast", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot 3: Significant Transition Comparisons
plot3 <- ggplot(transition_data_sig, aes(x = Contrast, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "darkgreen")) +
  scale_x_discrete(limits = transition_contrasts) +
  labs(title = "Significant DE Genes Over Transitions", x = "Contrast", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot 4: Non-significant Transition Comparisons
plot4 <- ggplot(transition_data_none, aes(x = Contrast, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("None" = "grey")) +
  scale_x_discrete(limits = transition_contrasts) +
  labs(title = "Non-significant DE Genes Over Transitions", x = "Contrast", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Adjusting plot commands to include counts above bars
plot1 <- plot1 +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size = 3)

plot2 <- plot2 +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size = 3)

plot3 <- plot3 +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size = 3)

plot4 <- plot4 +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.9), size = 3)


#Displaying the Plots
## If you want to display them individually
print(plot1)
print(plot2)
print(plot3)
print(plot4)

# Or use gridExtra to arrange them in a grid if you prefer seeing them together
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)

grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)


###################################volcano_plots###################################################

# Specify the contrast order
time_point_contrasts <- c("g_vs_w_0h", "g_vs_w_3h", "g_vs_w_6h", "g_vs_w_9h", "g_vs_w_12h", "g_vs_w_15h", "g_0h_to_3h", "g_12h_to_15h", "w_0h_to_3h", "w_12h_to_15h")

# Initialize lists to store plots
volcano_plots_ordered <- list()
ma_plots_ordered <- list()

# Generate plots for each contrast in the specified order
for (contrast in time_point_contrasts) {
  data <- subset(combined_de_genes, Contrast == contrast)
  
  # Categorization for coloring
  data$Category <- with(data, ifelse(adj.P.Val < 0.05 & logFC > 0, "Upregulated",
                                     ifelse(adj.P.Val < 0.05 & logFC < 0, "Downregulated", "None")))
  
  # Volcano Plot
  volcano_plot <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), color = Category)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "darkgreen", "None" = "grey")) +
    ggtitle(paste("Volcano Plot:", contrast)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  volcano_plots_ordered[[contrast]] <- volcano_plot
  
  # MA Plot
  ma_plot <- ggplot(data, aes(x = AveExpr, y = logFC, color = Category)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "darkgreen", "None" = "grey")) +
    ggtitle(paste("MA Plot:", contrast)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  ma_plots_ordered[[contrast]] <- ma_plot
}


library(patchwork)

# Combine and arrange volcano plots
volcano_collage <- wrap_plots(volcano_plots_ordered, ncol = 2)  # Adjust 'ncol' as needed for your layout

# Combine and arrange MA plots
ma_collage <- wrap_plots(ma_plots_ordered, ncol = 2)  # Adjust 'ncol' as needed for your layout

# Display the collages
volcano_collage
ma_collage
