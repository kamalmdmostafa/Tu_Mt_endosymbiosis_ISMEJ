library(shiny)
library(tidyverse)
library(ggplot2)
library(viridis)
library(svglite)

# Function to calculate chlorophyll content (μg/mL)
calculate_chlorophyll <- function(A652, A665) {
  chl_a <- 16.29 * A665 - 8.54 * A652
  chl_b <- 30.66 * A652 - 13.58 * A665
  total_chl <- chl_a + chl_b
  
  return(tibble(chl_a = chl_a, chl_b = chl_b, total_chl = total_chl))
}

ui <- fluidPage(
  titlePanel("Chlorophyll Content Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      radioButtons("input_type", "Input Type:",
                   choices = c("Manual Entry" = "manual", "CSV Upload" = "csv")),
      
      conditionalPanel(
        condition = "input.input_type == 'manual'",
        numericInput("num_samples", "Number of Samples", value = 2, min = 1),
        uiOutput("manual_input_fields")
      ),
      
      conditionalPanel(
        condition = "input.input_type == 'csv'",
        fileInput("csv_file", "Choose CSV File",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
        tags$p("CSV file should have columns: Sample, A652, A665")
      ),
      
      actionButton("analyze", "Analyze Data"),
      
      br(), br(),
      
      downloadButton("download_svg_box", "Download Box Plot (SVG)"),
      downloadButton("download_pdf_box", "Download Box Plot (PDF)"),
      br(),
      downloadButton("download_svg_percent", "Download Percent Difference Plot (SVG)"),
      downloadButton("download_pdf_percent", "Download Percent Difference Plot (PDF)")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Instructions", 
                 tags$div(
                   tags$h3("How to Use This App"),
                   tags$p("1. Choose your input method: Manual Entry or CSV Upload."),
                   tags$h4("For Manual Entry:"),
                   tags$ul(
                     tags$li("Specify the number of samples you have."),
                     tags$li("For each sample, provide a name and paste A652 and A665 values."),
                     tags$li("You can paste values directly from Excel or Google Sheets."),
                     tags$li("Make sure the number of A652 values matches the number of A665 values for each sample.")
                   ),
                   tags$h4("For CSV Upload:"),
                   tags$ul(
                     tags$li("Prepare a CSV file with columns: Sample, A652, A665"),
                     tags$li("Upload your CSV file using the file input button.")
                   ),
                   tags$p("2. Click 'Analyze Data' to process your input."),
                   tags$p("3. View results in the Data, Box Plot, and Percent Difference Plot tabs."),
                   tags$p("4. Use the download buttons to save plots as SVG or PDF files.")
                 )),
        tabPanel("Data", tableOutput("data_table")),
        tabPanel("Box Plot", plotOutput("box_plot")),
        tabPanel("Percent Difference Plot", plotOutput("percent_plot")),
        tabPanel("Statistics", verbatimTextOutput("stats_output"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$manual_input_fields <- renderUI({
    num_samples <- input$num_samples
    
    lapply(1:num_samples, function(i) {
      tagList(
        h4(paste("Sample", i)),
        textInput(paste0("sample_name_", i), "Sample Name", value = paste0("Sample_", i)),
        textAreaInput(paste0("a652_", i), "A652 Values (paste directly from spreadsheet)", rows = 3),
        textAreaInput(paste0("a665_", i), "A665 Values (paste directly from spreadsheet)", rows = 3),
        hr()
      )
    })
  })
  
  data <- reactive({
    req(input$analyze)
    
    if (input$input_type == "manual") {
      num_samples <- input$num_samples
      
      data <- tibble(
        Sample = character(),
        A652 = numeric(),
        A665 = numeric()
      )
      
      for (i in 1:num_samples) {
        sample_name <- input[[paste0("sample_name_", i)]]
        a652_values <- as.numeric(strsplit(input[[paste0("a652_", i)]], "[\n\t,;]+")[[1]])
        a665_values <- as.numeric(strsplit(input[[paste0("a665_", i)]], "[\n\t,;]+")[[1]])
        
        if (length(a652_values) != length(a665_values)) {
          stop(paste("The number of A652 and A665 values for", sample_name, "do not match."))
        }
        
        for (j in seq_along(a652_values)) {
          data <- data %>% add_row(Sample = sample_name, A652 = a652_values[j], A665 = a665_values[j])
        }
      }
    } else {
      req(input$csv_file)
      data <- read_csv(input$csv_file$datapath)
    }
    
    data %>%
      mutate(calculate_chlorophyll(A652, A665)) %>%
      mutate(Group = str_extract(Sample, "^[^_]+"))
  })
  
  output$data_table <- renderTable({
    data()
  })
  
  plot_data <- reactive({
    data() %>%
      pivot_longer(cols = c(chl_a, chl_b, total_chl), 
                   names_to = "Chlorophyll_Type", 
                   values_to = "Value")
  })
  
  t_test_results <- reactive({
    plot_data() %>%
      group_by(Chlorophyll_Type) %>%
      summarize(p_value = t.test(Value ~ Group)$p.value,
                y_pos = max(Value) * 1.1)
  })
  
  create_box_plot <- function() {
    ggplot(plot_data(), aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.5) +
      facet_wrap(~ Chlorophyll_Type, scales = "free_y", 
                 labeller = labeller(Chlorophyll_Type = c(chl_a = "Chlorophyll a", 
                                                          chl_b = "Chlorophyll b", 
                                                          total_chl = "Total Chlorophyll"))) +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      labs(title = "Chlorophyll Content by Group",
           y = "Concentration (μg/mL)",
           x = "Group") +
      theme_minimal() +
      theme(
        legend.position = "bottom",
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      geom_text(data = t_test_results(), 
                aes(x = 1.5, y = y_pos, label = paste("p =", sprintf("%.5f", p_value))),
                inherit.aes = FALSE)
  }
  
  output$box_plot <- renderPlot({
    create_box_plot()
  })
  
  create_percent_plot <- function() {
    percent_diff <- data() %>%
      group_by(Group) %>%
      summarize(across(c(chl_a, chl_b, total_chl), mean)) %>%
      pivot_longer(cols = c(chl_a, chl_b, total_chl), names_to = "Chlorophyll_Type", values_to = "Value") %>%
      pivot_wider(names_from = Group, values_from = Value) %>%
      mutate(Percent_Difference = (.[[3]] - .[[2]]) / .[[2]] * 100)
    
    percent_diff_with_p <- percent_diff %>%
      left_join(t_test_results(), by = "Chlorophyll_Type")
    
    ggplot(percent_diff_with_p, aes(x = Chlorophyll_Type, y = Percent_Difference, fill = Chlorophyll_Type)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = sprintf("%.1f%%", Percent_Difference), y = Percent_Difference + 
                      sign(Percent_Difference) * 2), vjust = ifelse(percent_diff_with_p$Percent_Difference >= 0, -0.5, 1.5)) +
      geom_text(aes(label = paste("p =", sprintf("%.5f", p_value)), y = 0), vjust = -0.5) +
      scale_fill_viridis(discrete = TRUE, option = "D") +
      labs(title = "Percent Difference in Chlorophyll Content",
           subtitle = "With p-values from t-test",
           y = "Percent Difference (%)",
           x = "Chlorophyll Type") +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.line = element_line(color = "black", linewidth = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA)
      ) +
      scale_x_discrete(labels = c("Chlorophyll a", "Chlorophyll b", "Total Chlorophyll"))
  }
  
  output$percent_plot <- renderPlot({
    create_percent_plot()
  })
  
  output$stats_output <- renderPrint({
    t_test_results()
  })
  
  output$download_svg_box <- downloadHandler(
    filename = function() { "chlorophyll_boxplot.svg" },
    content = function(file) {
      svglite(file)
      print(create_box_plot())
      dev.off()
    }
  )
  
  output$download_pdf_box <- downloadHandler(
    filename = function() { "chlorophyll_boxplot.pdf" },
    content = function(file) {
      pdf(file)
      print(create_box_plot())
      dev.off()
    }
  )
  
  output$download_svg_percent <- downloadHandler(
    filename = function() { "chlorophyll_percent_difference.svg" },
    content = function(file) {
      svglite(file)
      print(create_percent_plot())
      dev.off()
    }
  )
  
  output$download_pdf_percent <- downloadHandler(
    filename = function() { "chlorophyll_percent_difference.pdf" },
    content = function(file) {
      pdf(file)
      print(create_percent_plot())
      dev.off()
    }
  )
}

shinyApp(ui = ui, server = server)