# Load required libraries
library(clusterProfiler)
library(org.Mtet3.eg.db)  # Adjust with your specific organism database
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(extras)


# Function to read gene list without requiring specific headers
read_gene_list <- function(file_path) {
  tryCatch({
    data <- read.csv(file_path, header = FALSE, stringsAsFactors = FALSE)
    genes <- as.character(data[,1])
    return(genes)
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(NULL)
  })
}

# Load gene data for all clusters
clusters <- list(
  #cluster_1 = read_gene_list("Kallisto_Sig_new_4_clusters_NewCluster_1_genelist.txt"),
  #cluster_2 = read_gene_list("Kallisto_Sig_new_4_clusters_NewCluster_2_genelist.txt"),
  Mt_up = read_gene_list("Kallisto_sig_up_genes.txt"),
  Mt_down = read_gene_list("Kallisto_sig_down_genes.txt")
  )

# Remove any NULL entries (failed to load)
clusters <- clusters[!sapply(clusters, is.null)]

# Load the universe gene set
universe_data <- read.delim("mt_universe.txt", sep = "\t", header = FALSE)
universe_genes <- as.character(universe_data$V1)

# Function to perform compareCluster analysis
perform_compareCluster <- function(gene_lists, ontology) {
  tryCatch({
    compare_result <- compareCluster(
      geneCluster = gene_lists,
      fun = "enrichGO",
      OrgDb = org.Mtet3.eg.db,
      keyType = "GID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.01,
      universe = universe_genes
    )
    
    # Simplify the results
    compare_result_simplified <- clusterProfiler::simplify(compare_result, cutoff = 0.7, by = "p.adjust", select_fun = min)
    
    return(list(full = compare_result, simplified = compare_result_simplified))
  }, error = function(e) {
    warning(paste("Error in compareCluster for", ontology, "ontology:", e$message))
    return(NULL)
  })
}

# Function to generate and visualize plots
generate_plots <- function(compare_result, ont, type) {
  if (is.null(compare_result) || nrow(compare_result@compareClusterResult) == 0) {
    warning(paste("No enrichment results for", ont, "ontology", type, ". Skipping plot generation."))
    return()
  }
  
  # 1. Dot plot
  tryCatch({
    dot_plot <- dotplot(compare_result, 
                        showCategory = 30) +  # Adjust this value to show more/fewer terms
      ggtitle(paste("GO Terms Dot Plot -", ont, "-", type)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 12),  # Adjust title size
        axis.text.y = element_text(size = 6, margin = margin(r = 5)),  # Adjust y-axis text size
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),  # Adjust x-axis text size and angle
        axis.title = element_text(size = 10),  # Adjust axis title size
        legend.text = element_text(size = 6),  # Adjust legend text size
        legend.title = element_text(size = 8),  # Adjust legend title size
        plot.margin = margin(5, 15, 5, 5, "pt"),  # Adjust plot margins
        panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.3),
        axis.ticks.length = unit(0, "pt"),
        panel.spacing = unit(0.1, "lines"),
        panel.spacing.x = unit(0, "lines")  # Adjust space between clusters (0 for minimum)
      ) +
      scale_y_discrete(labels = function(x) stringr::str_trunc(x, 45)) +  # Truncate long y-axis labels
      coord_cartesian(clip = "off") +
      scale_size(range = c(1, 4)) +  # Adjust dot size range
      scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Adjust x-axis expansion
      scale_y_discrete(expand = expansion(mult = c(0.01, 0.05)))  # Adjust y-axis expansion
    
    print(dot_plot)
    ggsave(paste0("GO_", ont, "_", type, "_dotplot.pdf"), plot = dot_plot, 
           width = 8, height = 10, limitsize = FALSE)  # Adjust width and height as needed
    cat("Dot plot generated for", ont, "ontology", type, "\n")
  }, error = function(e) {
    warning(paste("Error generating dot plot for", ont, "ontology", type, ":", e$message))
  })
  
  # 2. Enrichment Map
  tryCatch({
    # Check if there are enough significant terms for the enrichment map
    sig_terms <- compare_result@compareClusterResult[compare_result@compareClusterResult$p.adjust <= 0.05, ]
    if(nrow(sig_terms) < 2) {
      warning("Not enough significant terms to generate enrichment map.")
      return()
    }
    
    # Calculate pairwise similarity between terms
    compare_result_with_sim <- pairwise_termsim(compare_result)
    
    em_plot <- emapplot(compare_result_with_sim, 
                        showCategory = 30,  # Adjust this value to show more/fewer terms
                        max.overlaps = 20) +  # Increase max.overlaps to show more labels
      ggtitle(paste("GO Terms Enrichment Map -", ont, "-", type)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12))  # Adjust title size
    
    print(em_plot)
    ggsave(paste0("GO_", ont, "_", type, "_enrichment_map.pdf"), plot = em_plot, 
           width = 12, height = 10, limitsize = FALSE)  # Adjust width and height as needed
    cat("Enrichment map generated for", ont, "ontology", type, "\n")
  }, error = function(e) {
    warning(paste("Error generating enrichment map for", ont, "ontology", type, ":", e$message))
    print(e)  # Print the full error message for debugging
  })
  
  # 3. Category Net Plot (cnetplot)
  tryCatch({
    cnet_plot <- cnetplot(compare_result, 
                          categorySize = "pvalue", 
                          showCategory = 5,  # Adjust this value to show more/fewer categories
                          colorEdge = TRUE)
    
    # Extract plot data for label adjustment
    plot_data <- cnet_plot$data
    go_terms <- plot_data[nchar(as.character(plot_data$name)) > 10, ]
    
    final_cnet_plot <- cnet_plot +
      geom_label(data = go_terms, 
                 aes(x = x, y = y, label = name),
                 fill = "white",
                 color = "black",
                 label.size = 0.25,
                 label.padding = unit(0.2, "lines"),
                 size = 5) +  # Adjust this value for larger/smaller category names
      scale_size_continuous(range = c(5, 15)) +  # Adjust this range for larger/smaller nodes
      scale_color_continuous(low = "red", high = "blue", name = "p.adjust") +
      ggtitle(paste("GO Terms Category Net Plot -", ont, "-", type)) +
      theme(plot.title = element_text(hjust = 0.5, size = 12),
            # Reduce space between x-axis KMC clusters
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8, margin = margin(t = 0, r = 0, b = 0, l = 0)),
            panel.spacing.x = unit(0.1, "lines"))  # Adjust this value to change space between facets
    
    print(final_cnet_plot)
    ggsave(paste0("GO_", ont, "_", type, "_cnetplot.pdf"), plot = final_cnet_plot, 
           width = 14, height = 12, limitsize = FALSE)  # Adjust width and height as needed
    cat("Category net plot generated for", ont, "ontology", type, "\n")
  }, error = function(e) {
    warning(paste("Error generating category net plot for", ont, "ontology", type, ":", e$message))
  })
}

# Function to export results to Excel
export_to_excel <- function(compare_results, ontology) {
  if (is.null(compare_results$full) || nrow(compare_results$full@compareClusterResult) == 0) {
    warning(paste("No enrichment results for", ontology, "ontology. Skipping Excel export."))
    return()
  }
  
  tryCatch({
    wb <- createWorkbook()
    
    # All results
    addWorksheet(wb, "All_Enrichment_Results")
    writeData(wb, "All_Enrichment_Results", as.data.frame(compare_results$full))
    
    # Significant results
    significant <- compare_results$full@compareClusterResult[compare_results$full@compareClusterResult$p.adjust <= 0.05, ]
    addWorksheet(wb, "Significant_Enrichment_Results")
    writeData(wb, "Significant_Enrichment_Results", significant)
    
    # Simplified significant results
    if (!is.null(compare_results$simplified)) {
      simplified_significant <- compare_results$simplified@compareClusterResult[compare_results$simplified@compareClusterResult$p.adjust <= 0.05, ]
      addWorksheet(wb, "Simplified_Significant_Results")
      writeData(wb, "Simplified_Significant_Results", simplified_significant)
    }
    
    saveWorkbook(wb, paste0("GO_", ontology, "_enrichment_results.xlsx"), overwrite = TRUE)
  }, error = function(e) {
    warning(paste("Error exporting results to Excel for", ontology, "ontology:", e$message))
  })
}

# Perform analysis for each ontology
ontologies <- c("BP", "MF", "CC")

for (ont in ontologies) {
  cat("Performing analysis for", ont, "ontology...\n")
  
  compare_results <- perform_compareCluster(clusters, ont)
  
  if (!is.null(compare_results) && 
      !is.null(compare_results$full) && 
      nrow(compare_results$full@compareClusterResult) > 0) {
    
    # Generate plots for full results
    generate_plots(compare_results$full, ont, "Full")
    
    # Generate plots for simplified results
    if (!is.null(compare_results$simplified)) {
      generate_plots(compare_results$simplified, ont, "Simplified")
    }
    
    # Export results to Excel
    export_to_excel(compare_results, ont)
    
    cat("Analysis for", ont, "ontology completed.\n\n")
  } else {
    cat("No enrichment found for", ont, "ontology. Skipping to next ontology.\n\n")
  }
}

cat("All analyses completed. Check the Excel files and PDF plot files for results.\n")