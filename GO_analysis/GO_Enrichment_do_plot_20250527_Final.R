# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)
library(writexl)
library(extrafont)

# Load Arial font
loadfonts(device = "win")  # Use "mac" instead of "win" if you're on a Mac

# Define distinct colors for 30 categories (original vector preserved)
distinct_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#1B9E77", "#D95F02",
  "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D",
  "#00CED1", "#FF1493", "#32CD32", "#9370DB", "#8B4513",
  "#DC143C", "#00FF7F", "#4B0082", "#CD853F", "#800000",
  "#008B8B", "#9932CC", "#FFB6C1", "#556B2F", "#2F4F4F"
)

# Function to read data and calculate Enrichment Score
calculate_enrichment_score <- function(file_path) {
  data <- read_excel(file_path)
  
  required_columns <- c("Cluster", "Description", "GeneRatio", "BgRatio", "p.adjust", "Count")
  missing_columns <- setdiff(required_columns, colnames(data))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))
  }
  
  data <- data %>%
    mutate(
      GeneRatio_k = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 1)),
      GeneRatio_n = as.numeric(sapply(strsplit(GeneRatio, "/"), `[`, 2)),
      BgRatio_K = as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 1)),
      BgRatio_N = as.numeric(sapply(strsplit(BgRatio, "/"), `[`, 2)),
      GeneRatio_numeric = GeneRatio_k / GeneRatio_n,
      BgRatio_numeric = BgRatio_K / BgRatio_N,
      EnrichmentScore = GeneRatio_numeric / BgRatio_numeric,
      p.adjust = as.numeric(p.adjust),
      Count = as.numeric(Count)
    )
  
  return(data)
}

# Function to get top n terms per cluster with sorting option
get_top_n_terms <- function(data, n, sort_by = "p.adjust") {
  data %>%
    group_by(Cluster) %>%
    slice_min(order_by = if(sort_by == "p.adjust") p.adjust else -EnrichmentScore, 
              n = n, with_ties = FALSE) %>%
    ungroup()
}

# Function to filter by enrichment score
filter_by_enrichment <- function(data, cutoff) {
  filtered_data <- data %>%
    filter(EnrichmentScore >= cutoff)
  
  # Print summary statistics with na.rm = TRUE
  cat("\nEnrichment Score Summary:")
  cat("\nTotal terms before filtering:", nrow(data))
  cat("\nTerms passing cutoff (ES ≥", format(cutoff, digits = 2), "):", nrow(filtered_data))
  cat("\nRange of Enrichment Scores:")
  cat("\n  Min:", format(min(data$EnrichmentScore, na.rm = TRUE), digits = 2))
  cat("\n  Max:", format(max(data$EnrichmentScore, na.rm = TRUE), digits = 2))
  cat("\n  Median:", format(median(data$EnrichmentScore, na.rm = TRUE), digits = 2), "\n")
  
  return(filtered_data)
}

# Function to filter by gene count
filter_by_gene_count <- function(data, min_count) {
  filtered_data <- data %>%
    filter(Count >= min_count)
  
  # Print summary statistics
  cat("\nGene Count Summary:")
  cat("\nTotal terms before filtering:", nrow(data))
  cat("\nTerms passing minimum count (≥", min_count, "):", nrow(filtered_data))
  cat("\nRange of Gene Counts:")
  cat("\n  Min:", min(data$Count, na.rm = TRUE))
  cat("\n  Max:", max(data$Count, na.rm = TRUE))
  cat("\n  Median:", median(data$Count, na.rm = TRUE), "\n")
  
  return(filtered_data)
}

# Filter by selected clusters
filter_by_clusters <- function(data, selected_clusters) {
  filtered_data <- data %>%
    filter(Cluster %in% selected_clusters)
  
  # Print summary statistics
  cat("\nCluster Selection Summary:")
  cat("\nTotal terms before filtering:", nrow(data))
  cat("\nTerms in selected clusters:", nrow(filtered_data))
  cat("\nSelected clusters:", paste(selected_clusters, collapse = ", "), "\n")
  
  return(filtered_data)
}

# Function to create plot with category ordering options and distinct colors
create_plot <- function(data, width = 4, height = 8, category_order = "data", category_display = "legend") {
  # Truncate long descriptions
  data <- data %>%
    mutate(Description = str_trunc(Description, 50, "right"))
  
  # Preserve exact cluster order
  cluster_order <- unique(data$Cluster)
  data$Cluster <- factor(data$Cluster, levels = cluster_order)
  
  if ("User_Defined_Category" %in% colnames(data)) {
    # Fix for category ordering
    if (category_order == "data") {
      # Preserve the order as it appears in the data
      category_order_levels <- unique(data$User_Defined_Category)
    } else if (category_order == "plot") {
      # Order by max enrichment score
      data_ordered <- data %>%
        group_by(User_Defined_Category) %>%
        summarise(max_score = max(EnrichmentScore, na.rm = TRUE)) %>%
        arrange(desc(max_score))
      category_order_levels <- data_ordered$User_Defined_Category
    }
    
    # Create the ordered dataset - fixing ordering issues
    data <- data %>%
      mutate(User_Defined_Category = factor(User_Defined_Category, levels = category_order_levels))
    
    # Order descriptions by enrichment score within categories
    ordered_data <- data %>%
      arrange(User_Defined_Category, desc(EnrichmentScore))
    
    # Set the factor levels for descriptions
    data <- data %>%
      mutate(Description = factor(Description, levels = rev(unique(ordered_data$Description))))
    
    # Set up category aesthetics
    tile_aes <- aes(fill = User_Defined_Category)
    
    # Set up color palette with distinct colors
    fill_scale <- scale_fill_manual(
      name = "Category",
      values = distinct_colors[1:length(category_order_levels)],
      limits = category_order_levels
    )
    
    # Additional plot elements for right-side text display
    if (category_display == "text") {
      # Create a dataframe with representative rows for each category
      # Get coordinates for each category's position
      category_data <- data %>%
        group_by(User_Defined_Category) %>%
        summarise(
          y_position = which(levels(data$Description) == as.character(Description[1]))[1]
        ) %>%
        ungroup()
      
      # Account for possible NA values
      category_data <- category_data[!is.na(category_data$y_position), ]
      
      # Map category to colors
      category_colors <- setNames(
        distinct_colors[1:length(unique(data$User_Defined_Category))],
        unique(data$User_Defined_Category)
      )
    }
  } else {
    data <- data %>%
      arrange(desc(EnrichmentScore)) %>%
      mutate(Description = factor(Description, levels = unique(Description)))
    
    tile_aes <- aes(fill = "lightblue")
    fill_scale <- scale_fill_manual(values = c("lightblue"), name = "Category")
    category_display <- "legend" # Force legend display if no categories
  }
  
  # Base plot
  plot <- ggplot(data, aes(x = Cluster, y = Description)) +
    geom_tile(tile_aes, alpha = 0.3, width = Inf) +
    geom_point(aes(color = p.adjust, size = EnrichmentScore)) +
    fill_scale +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted p-value") +
    scale_size(range = c(1, 4), name = "Enrichment Score") +
    ggtitle("GO Terms Dot Plot - BP - Simplified")
  
  # Complete the plot with appropriate theme
  plot <- plot + 
    theme_minimal() +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(hjust = 0.5, size = 12),
      axis.text.y = element_text(size = 7, margin = margin(r = 2), color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 6),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 6),
      plot.margin = margin(5, ifelse(category_display == "text", 80, 15), 5, 5, "pt"),
      panel.grid.major.y = element_line(color = "lightgrey", linewidth = 0.2),
      axis.ticks.length = unit(0, "pt"),
      panel.spacing = unit(0.05, "lines"),
      panel.spacing.x = unit(0, "lines"),
      legend.position = "right",
      legend.key.size = unit(0.5, "lines"),
      legend.box.just = "left",
      legend.justification = "left",
      legend.box.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.box = "vertical",
      legend.spacing = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0)
    ) +
    guides(
      color = guide_colorbar(order = 1, title = "p.adjust"),
      size = guide_legend(order = 2, title = "Enrichment Score"),
      fill = guide_legend(order = 3, title = "Category", ncol = 1)
    ) +
    coord_cartesian(clip = "off") +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0)))
  
  # Add annotations for right-side text display
  if (category_display == "text" && "User_Defined_Category" %in% colnames(data)) {
    # Add an empty column for the category labels
    last_cluster <- tail(levels(data$Cluster), 1)
    data_x_pos <- max(as.numeric(factor(last_cluster))) + 0.5
    
    # Add annotations for each category
    for (i in 1:nrow(category_data)) {
      cat_name <- as.character(category_data$User_Defined_Category[i])
      y_pos <- category_data$y_position[i]
      cat_color <- category_colors[cat_name]
      
      # Add annotation
      plot <- plot + 
        annotate("text", 
                 x = data_x_pos, 
                 y = y_pos, 
                 label = cat_name, 
                 color = cat_color,
                 hjust = -0.1,
                 fontface = "bold")
    }
    
    # Add "Category" title
    plot <- plot + 
      annotate("text", 
               x = data_x_pos, 
               y = nrow(data) + 1,
               label = "Category", 
               fontface = "bold",
               hjust = -0.1)
  }
  
  return(plot)
}

# Function to export data with custom filename
export_data <- function(data, prompt = TRUE) {
  if (prompt) {
    filename <- readline("Enter filename for data export (default: enrichment_data.xlsx): ")
    if (filename == "") filename <- "enrichment_data.xlsx"
    if (!grepl("\\.xlsx$", filename)) filename <- paste0(filename, ".xlsx")
  } else {
    filename <- "enrichment_data.xlsx"
  }
  
  tryCatch({
    write_xlsx(data, path = filename)
    cat("Data exported to", filename, "\n")
  }, error = function(e) {
    cat("Error exporting data:", e$message, "\n")
  })
}

# Function to handle visualization loop with original menu structure
handle_visualization <- function(enriched_data) {
  current_plot <- NULL
  current_filtered_data <- NULL
  current_cutoff <- NULL
  current_min_count <- NULL
  current_category_order <- "data"
  current_sort_by <- "p.adjust"  # Default sorting by p.adjust
  current_selected_clusters <- NULL  # Track selected clusters
  current_category_display <- "legend"  # Default category display method
  
  while (TRUE) {
    # Show current enrichment score range with na.rm = TRUE
    cat("\nEnrichment Score Statistics:")
    cat("\nMin:", format(min(enriched_data$EnrichmentScore, na.rm = TRUE), digits = 2))
    cat("\nMax:", format(max(enriched_data$EnrichmentScore, na.rm = TRUE), digits = 2))
    cat("\nMedian:", format(median(enriched_data$EnrichmentScore, na.rm = TRUE), digits = 2))
    
    if (!is.null(current_cutoff)) {
      cat("\nCurrent cutoff: ≥", format(current_cutoff, digits = 2))
    }
    if (!is.null(current_min_count)) {
      cat("\nCurrent minimum gene count: ≥", current_min_count)
    }
    if (!is.null(current_selected_clusters)) {
      cat("\nSelected clusters: ", paste(current_selected_clusters, collapse = ", "))
    }
    cat("\nCurrent sorting method:", current_sort_by)
    cat("\nCurrent category display:", ifelse(current_category_display == "legend", "Legend", "Right side text"))
    
    cat("\n\nEnter action:")
    cat("\n1: Enter new enrichment score cutoff")
    cat("\n2: Set minimum gene count")
    cat("\n3: Change sorting method (p.adjust/EnrichmentScore)")
    cat("\n4: Change category order (data order/plot appearance order)")
    cat("\n5: Toggle category display (legend/right side text)")
    cat("\n6: Select specific clusters to visualize")
    cat("\n7: Plot in R Studio")
    cat("\n8: Save current plot")
    cat("\n9: Export data")
    cat("\n10: Quit visualization\n")
    
    action <- readline("Choose action (1-10): ")
    
    if (action == "1") {
      # Enrichment score cutoff
      cutoff_input <- readline("Enter enrichment score cutoff: ")
      cutoff <- as.numeric(cutoff_input)
      
      if (is.na(cutoff) || cutoff < 0) {
        cat("\nInvalid cutoff. Please enter a non-negative number.\n")
        next
      }
      
      filtered_data <- filter_by_enrichment(enriched_data, cutoff)
      if (!is.null(current_min_count)) {
        filtered_data <- filter_by_gene_count(filtered_data, current_min_count)
      }
      if (!is.null(current_selected_clusters)) {
        filtered_data <- filter_by_clusters(filtered_data, current_selected_clusters)
      }
      
      if (nrow(filtered_data) == 0) {
        cat("\nNo terms pass this cutoff. Please try a lower value.\n")
        next
      }
      
      current_filtered_data <- filtered_data
      top_terms <- get_top_n_terms(filtered_data, 30, current_sort_by)
      current_plot <- create_plot(top_terms, category_order = current_category_order, 
                                  category_display = current_category_display)
      current_cutoff <- cutoff
      
      cat("\nCutoff applied. Use option 7 to display the plot in RStudio.\n")
      
    } else if (action == "2") {
      # Minimum gene count
      count_input <- readline("Enter minimum number of genes per GO category: ")
      min_count <- as.numeric(count_input)
      
      if (is.na(min_count) || min_count < 0) {
        cat("\nInvalid count. Please enter a non-negative number.\n")
        next
      }
      
      filtered_data <- enriched_data
      if (!is.null(current_cutoff)) {
        filtered_data <- filter_by_enrichment(filtered_data, current_cutoff)
      }
      if (!is.null(current_selected_clusters)) {
        filtered_data <- filter_by_clusters(filtered_data, current_selected_clusters)
      }
      filtered_data <- filter_by_gene_count(filtered_data, min_count)
      
      if (nrow(filtered_data) == 0) {
        cat("\nNo terms pass this gene count threshold. Please try a lower value.\n")
        next
      }
      
      current_filtered_data <- filtered_data
      top_terms <- get_top_n_terms(filtered_data, 30, current_sort_by)
      current_plot <- create_plot(top_terms, category_order = current_category_order, 
                                  category_display = current_category_display)
      current_min_count <- min_count
      
      cat("\nGene count filter applied. Use option 7 to display the plot in RStudio.\n")
      
    } else if (action == "3") {
      # Sorting method
      cat("\nSelect sorting method:")
      cat("\n1: Sort by adjusted p-value (lowest first)")
      cat("\n2: Sort by enrichment score (highest first)\n")
      sort_choice <- readline("Choose sorting method (1-2): ")
      
      if (sort_choice == "1") {
        current_sort_by <- "p.adjust"
      } else if (sort_choice == "2") {
        current_sort_by <- "EnrichmentScore"
      } else {
        cat("\nInvalid choice. Keeping current sorting method.\n")
        next
      }
      
      if (!is.null(current_filtered_data)) {
        top_terms <- get_top_n_terms(current_filtered_data, 30, current_sort_by)
        current_plot <- create_plot(top_terms, category_order = current_category_order,
                                    category_display = current_category_display)
        cat("\nPlot updated with new sorting method. Use option 7 to display.\n")
      } else {
        cat("\nNo filtered data available. Please set a cutoff first using option 1 or 2.\n")
      }
      
    } else if (action == "4") {
      # Category order
      cat("\nSelect category order:")
      cat("\n1: Order by data appearance")
      cat("\n2: Order by plot appearance\n")
      order_choice <- readline("Choose order (1-2): ")
      
      if (order_choice == "1") {
        current_category_order <- "data"
      } else if (order_choice == "2") {
        current_category_order <- "plot"
      } else {
        cat("\nInvalid choice. Keeping current order.\n")
        next
      }
      
      if (!is.null(current_filtered_data)) {
        top_terms <- get_top_n_terms(current_filtered_data, 30, current_sort_by)
        current_plot <- create_plot(top_terms, category_order = current_category_order,
                                    category_display = current_category_display)
        cat("\nPlot updated with new category order. Use option 7 to display.\n")
      } else {
        cat("\nNo filtered data available. Please set a cutoff first using option 1 or 2.\n")
      }
      
    } else if (action == "5") {
      # Toggle category display
      cat("\nSelect category display method:")
      cat("\n1: Display categories in legend (default)")
      cat("\n2: Display categories as text on right side of plot\n")
      
      display_choice <- readline("Choose display method (1-2): ")
      
      if (display_choice == "1") {
        current_category_display <- "legend"
        cat("\nCategory display set to legend.\n")
      } else if (display_choice == "2") {
        current_category_display <- "text"
        cat("\nCategory display set to right side text.\n")
      } else {
        cat("\nInvalid choice. Keeping current display method.\n")
        next
      }
      
      if (!is.null(current_filtered_data)) {
        top_terms <- get_top_n_terms(current_filtered_data, 30, current_sort_by)
        current_plot <- create_plot(top_terms, category_order = current_category_order, 
                                    category_display = current_category_display)
        cat("\nPlot updated with new category display. Use option 7 to display.\n")
      } else {
        cat("\nNo filtered data available. Please set a cutoff first using option 1 or 2.\n")
      }
      
    } else if (action == "6") {
      # Cluster selection
      all_clusters <- unique(enriched_data$Cluster)
      cat("\nAvailable clusters:\n")
      for (i in seq_along(all_clusters)) {
        cat(i, ": ", all_clusters[i], "\n", sep = "")
      }
      
      cat("\nSelect action:")
      cat("\n1: Select specific clusters")
      cat("\n2: Select all clusters (reset selection)\n")
      
      cluster_action <- readline("Choose action (1-2): ")
      
      if (cluster_action == "1") {
        selections <- readline("Enter cluster numbers (comma-separated, e.g., 1,3,5): ")
        selection_indices <- as.numeric(strsplit(selections, ",")[[1]])
        
        # Validate selection indices
        valid_indices <- selection_indices[!is.na(selection_indices) & 
                                             selection_indices > 0 & 
                                             selection_indices <= length(all_clusters)]
        
        if (length(valid_indices) == 0) {
          cat("\nNo valid cluster selected. Please try again.\n")
          next
        }
        
        selected_clusters <- all_clusters[valid_indices]
        cat("\nSelected clusters:", paste(selected_clusters, collapse = ", "), "\n")
        
        filtered_data <- enriched_data
        if (!is.null(current_cutoff)) {
          filtered_data <- filter_by_enrichment(filtered_data, current_cutoff)
        }
        if (!is.null(current_min_count)) {
          filtered_data <- filter_by_gene_count(filtered_data, current_min_count)
        }
        
        filtered_data <- filter_by_clusters(filtered_data, selected_clusters)
        
        if (nrow(filtered_data) == 0) {
          cat("\nSelected clusters contain no terms that pass other filters. Please try different clusters or adjust other filters.\n")
          next
        }
        
        current_filtered_data <- filtered_data
        top_terms <- get_top_n_terms(filtered_data, 30, current_sort_by)
        current_plot <- create_plot(top_terms, category_order = current_category_order,
                                    category_display = current_category_display)
        current_selected_clusters <- selected_clusters
        
        cat("\nCluster selection applied. Use option 7 to display the plot in RStudio.\n")
        
      } else if (cluster_action == "2") {
        current_selected_clusters <- NULL
        
        filtered_data <- enriched_data
        if (!is.null(current_cutoff)) {
          filtered_data <- filter_by_enrichment(filtered_data, current_cutoff)
        }
        if (!is.null(current_min_count)) {
          filtered_data <- filter_by_gene_count(filtered_data, current_min_count)
        }
        
        current_filtered_data <- filtered_data
        top_terms <- get_top_n_terms(filtered_data, 30, current_sort_by)
        current_plot <- create_plot(top_terms, category_order = current_category_order,
                                    category_display = current_category_display)
        
        cat("\nAll clusters selected. Use option 7 to display the plot in RStudio.\n")
        
      } else {
        cat("\nInvalid choice. Returning to main menu.\n")
      }
      
    } else if (action == "7") {
      # Plot in R Studio
      if (is.null(current_plot)) {
        cat("\nNo plot available. Please create a plot first using options 1 or 2.\n")
        next
      }
      
      tryCatch({
        dev.off()
      }, error = function(e) {})
      
      graphics.off()
      dev.new()
      print(current_plot)
      cat("\nPlot displayed in RStudio.\n")
      
    } else if (action == "8") {
      # Save plot
      if (is.null(current_plot)) {
        cat("\nNo plot available. Please create a plot first using options 1 or 2.\n")
        next
      }
      
      filename <- readline("Enter filename for plot (default: GO_BP_Simplified_dotplot.pdf): ")
      if (filename == "") filename <- "GO_BP_Simplified_dotplot.pdf"
      if (!grepl("\\.pdf$", filename)) filename <- paste0(filename, ".pdf")
      
      tryCatch({
        ggsave(filename, plot = current_plot, width = 6, height = 6, limitsize = FALSE)
        cat("Plot saved as", filename, "\n")
      }, error = function(e) {
        cat("Error saving plot:", e$message, "\n")
      })
      
    } else if (action == "9") {
      # Export data
      export_data(enriched_data, prompt = TRUE)
      cat("\nData exported.\n")
      
    } else if (action == "10") {
      # Quit
      return(list(continue = FALSE, plot = current_plot, filtered_data = current_filtered_data))
      
    } else {
      cat("\nInvalid choice. Please enter a number between 1 and 10.\n")
    }
  }
}

# Main execution
main <- function() {
  file_path <- 'Revigo_reduced_for_enrichment.xlsx'
  
  cat("Step 1: Calculating Enrichment Score\n")
  enriched_data <- calculate_enrichment_score(file_path)
  
  cat("\nFirst few rows of the enriched data:\n")
  print(head(enriched_data))
  cat("\nColumns in the enriched data:\n")
  print(colnames(enriched_data))
  
  cat("\nStep 2: Interactive Visualization\n")
  vis_results <- handle_visualization(enriched_data)
  
  if (vis_results$continue) {
    cat("\nStep 3: Exporting Data\n")
    export_data(enriched_data, prompt = TRUE)
    cat("\nAnalysis complete.\n")
  } else {
    cat("\nAnalysis stopped after visualization.\n")
  }
  
  return(list(
    complete_data = enriched_data,
    filtered_data = vis_results$filtered_data,
    last_plot = vis_results$plot
  ))
}

# Run the main function
results <- main()