# ============================================================
# app.R
# Interactive NGC + UC v2
# Plots triaxial data
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(DT)

source("R/functions_NGC_UC.R")
source("R/functions_io.R")


# ============================================================
# UI
# ============================================================

ui <- fluidPage(
  
  titlePanel("Combined NGC + UC Failure Criterion"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      # ------------------------------------------------------
      # Input strengths
      # ------------------------------------------------------
      
      h4("Input strengths"),
      
      numericInput(
        inputId = "sigma_c",
        label = HTML("&sigma;<sub>c</sub> / UCS (MPa)"),
        value = 160,
        min = 1,
        step = 5
      ),
      
      numericInput(
        inputId = "sigma_t",
        label = HTML("&sigma;<sub>t</sub> / UTS (MPa, negative)"),
        value = -20,
        max = -0.1,
        step = -1
      ),
      
      hr(),
      
      
      # ------------------------------------------------------
      # Curve density
      # ------------------------------------------------------
      
      h4("Curve density"),
      
      sliderInput(
        inputId = "n_ngc",
        label = "Number of NGC points",
        min = 200,
        max = 5000,
        value = 1500,
        step = 100
      ),
      
      sliderInput(
        inputId = "n_uc",
        label = "Number of UC points",
        min = 200,
        max = 8000,
        value = 2500,
        step = 100
      ),
      
      hr(),
      
      
      # ------------------------------------------------------
      # UC truncation
      # ------------------------------------------------------
      
      h4("UC truncation"),
      
      checkboxInput(
        inputId = "use_ductile",
        label = "Use ductile-transition limit",
        value = TRUE
      ),
      
      numericInput(
        inputId = "ratio_limit",
        label = HTML(
          "&sigma;<sub>1</sub> / &sigma;<sub>3</sub> minimum ratio"
        ),
        value = 3.4,
        min = 1.1,
        step = 0.1
      ),
      
      hr(),
      
      
      # ------------------------------------------------------
      # Plot colours
      # ------------------------------------------------------
      
      h4("Plot colours"),
      
      selectInput(
        inputId = "color_ngc",
        label = "NGC colour",
        choices = c(
          "blue",
          "red",
          "darkgreen",
          "black",
          "purple"
        ),
        selected = "blue"
      ),
      
      selectInput(
        inputId = "color_uc",
        label = "UC colour",
        choices = c(
          "orange",
          "red",
          "darkgreen",
          "black",
          "purple"
        ),
        selected = "orange"
      ),
      
      hr(),
      
      
      # ------------------------------------------------------
      # Plot limits
      # ------------------------------------------------------
      
      h4("Plot limits"),
      
      checkboxInput(
        inputId = "manual_limits",
        label = "Use manual plot limits",
        value = FALSE
      ),
      
      numericInput(
        "x_sigma_tau_min",
        "σn-τ x min",
        value = -25
      ),
      
      numericInput(
        "x_sigma_tau_max",
        "σn-τ x max",
        value = 160
      ),
      
      numericInput(
        "y_sigma_tau_min",
        "σn-τ y min",
        value = 0
      ),
      
      numericInput(
        "y_sigma_tau_max",
        "σn-τ y max",
        value = 120
      ),
      
      numericInput(
        "x_s3_s1_min",
        "σ3-σ1 x min",
        value = -40
      ),
      
      numericInput(
        "x_s3_s1_max",
        "σ3-σ1 x max",
        value = 100
      ),
      
      numericInput(
        "y_s3_s1_min",
        "σ3-σ1 y min",
        value = 0
      ),
      
      numericInput(
        "y_s3_s1_max",
        "σ3-σ1 y max",
        value = 480
      ),
      
      hr(),
      
      
      # ------------------------------------------------------
      # Experimental data
      # ------------------------------------------------------
      
      h4("Experimental data"),
      
      fileInput(
        inputId = "triaxial_file",
        label = "NGC-Ucar dataset",
        accept = c(".csv")
      ),
      
      verbatimTextOutput("dataset_info"),
      
      hr(),
      
      
      # ------------------------------------------------------
      # Downloads
      # ------------------------------------------------------
      
      downloadButton(
        outputId = "download_combined_csv",
        label = "Download combined curve CSV"
      ),
      
      br(),
      br(),
      
      downloadButton(
        outputId = "download_parameters_csv",
        label = "Download parameters CSV"
      )
    ),
    
    
    # ========================================================
    # Main panel
    # ========================================================
    
    mainPanel(
      
      tabsetPanel(
        
        tabPanel(
          "Plots",
          br(),
          plotOutput(
            "plot_sigma_tau",
            height = "600px"
          ),
          br(),
          plotOutput(
            "plot_sigma3_sigma1",
            height = "600px"
          )
        ),
        
        tabPanel(
          "Parameters",
          br(),
          DTOutput("parameters_table")
        ),
        
        tabPanel(
          "Combined curve",
          br(),
          DTOutput("combined_table")
        ),
        
        tabPanel(
          "About",
          br(),
          
          h4("Purpose"),
          
          p(
            "This prototype calculates and displays the combined New Griffith Criterion ",
            "(NGC) and Ucar Criterion (UC) failure envelope from user-defined tensile ",
            "and compressive rock strengths."
          ),
          
          p(
            "The calculation is based on the original Script 01, refactored into reusable ",
            "functions suitable for an interactive Shiny/Shinylive application."
          )
        )
      )
    )
  )
)


# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {
  
  
  # ----------------------------------------------------------
  # NGC + UC calculation
  # ----------------------------------------------------------
  
  result <- reactive({
    
    validate(
      need(
        input$sigma_c > 0,
        "sigma_c must be positive."
      ),
      
      need(
        input$sigma_t < 0,
        "sigma_t must be negative."
      ),
      
      need(
        abs(input$sigma_t) < input$sigma_c,
        "abs(sigma_t) must be lower than sigma_c."
      )
    )
    
    calculate_ngc_uc(
      sigma_t_adj = input$sigma_t,
      sigma_c_adj = input$sigma_c,
      n_ngc = input$n_ngc,
      n_uc = input$n_uc,
      use_ductile_transition_limit = input$use_ductile,
      sigma1_sigma3_min_ratio = input$ratio_limit
    )
  })
  
  
  # ----------------------------------------------------------
  # Experimental dataset
  # ----------------------------------------------------------
  
  experimental_dataset <- reactive({
    
    req(input$triaxial_file)
    
    read_ngc_toolkit_dataset(
      input$triaxial_file$datapath
    )
  })
  
  
  # ----------------------------------------------------------
  # Experimental dataset information
  # ----------------------------------------------------------
  
  output$dataset_info <- renderText({
    
    if (is.null(input$triaxial_file)) {
      return("No dataset loaded.")
    }
    
    ds <- experimental_dataset()
    
    paste(
      "Dataset :", ds$dataset_id,
      "\nUnits   :", ds$units,
      "\nTests   :", ds$n_tests
    )
  })
  
  
  # ----------------------------------------------------------
  # Plot limits
  # ----------------------------------------------------------
  
  sigma_tau_xlim <- reactive({
    
    if (isTRUE(input$manual_limits)) {
      
      c(
        input$x_sigma_tau_min,
        input$x_sigma_tau_max
      )
      
    } else {
      
      NULL
    }
  })
  
  
  sigma_tau_ylim <- reactive({
    
    if (isTRUE(input$manual_limits)) {
      
      c(
        input$y_sigma_tau_min,
        input$y_sigma_tau_max
      )
      
    } else {
      
      NULL
    }
  })
  
  
  sigma3_sigma1_xlim <- reactive({
    
    if (isTRUE(input$manual_limits)) {
      
      c(
        input$x_s3_s1_min,
        input$x_s3_s1_max
      )
      
    } else {
      
      NULL
    }
  })
  
  
  sigma3_sigma1_ylim <- reactive({
    
    if (isTRUE(input$manual_limits)) {
      
      c(
        input$y_s3_s1_min,
        input$y_s3_s1_max
      )
      
    } else {
      
      NULL
    }
  })
  
  
  # ----------------------------------------------------------
  # Plots
  # ----------------------------------------------------------
  
  output$plot_sigma_tau <- renderPlot({
    
    plot_ngc_uc_sigma_tau(
      result = result(),
      xlim = sigma_tau_xlim(),
      ylim = sigma_tau_ylim(),
      color_ngc = input$color_ngc,
      color_uc = input$color_uc
    )
  })
  
  
  output$plot_sigma3_sigma1 <- renderPlot({
    
    plot_ngc_uc_sigma3_sigma1(
      result = result(),
      xlim = sigma3_sigma1_xlim(),
      ylim = sigma3_sigma1_ylim(),
      color_ngc = input$color_ngc,
      color_uc = input$color_uc
    )
  })
  
  
  # ----------------------------------------------------------
  # Tables
  # ----------------------------------------------------------
  
  output$parameters_table <- renderDT({
    
    datatable(
      result()$parameters,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })
  
  
  output$combined_table <- renderDT({
    
    datatable(
      result()$curve_combined,
      rownames = FALSE,
      options = list(
        pageLength = 20,
        scrollX = TRUE
      )
    )
  })
  
  
  # ----------------------------------------------------------
  # Downloads
  # ----------------------------------------------------------
  
  output$download_combined_csv <- downloadHandler(
    
    filename = function() {
      
      paste0(
        "curve_combined_NGC_UC_sigma_c_",
        input$sigma_c,
        "_sigma_t_",
        abs(input$sigma_t),
        ".csv"
      )
    },
    
    content = function(file) {
      
      write.csv(
        result()$curve_combined,
        file,
        row.names = FALSE
      )
    }
  )
  
  
  output$download_parameters_csv <- downloadHandler(
    
    filename = function() {
      
      paste0(
        "parameters_NGC_UC_sigma_c_",
        input$sigma_c,
        "_sigma_t_",
        abs(input$sigma_t),
        ".csv"
      )
    },
    
    content = function(file) {
      
      write.csv(
        result()$parameters,
        file,
        row.names = FALSE
      )
    }
  )
}


# ============================================================
# Run application
# ============================================================

shinyApp(ui, server)