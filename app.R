# Pharmacokinetic Population Analysis Application
# Main Shiny Application

library(shiny)
library(DT)
library(ggplot2)
library(nlme)
library(scales)

# Source helper modules
source("R/data_validation.R")
source("R/nca_analysis.R")
source("R/pk_models.R")
source("R/plotting.R")

# ============================================================
# Species options with typical weight ranges
# ============================================================
species_options <- c(
  "Human", "Dog", "Cat", "Horse", "Cattle",
  "Rat", "Mouse", "Rabbit", "Pig", "Sheep", "Goat",
  "Non-human Primate", "Guinea Pig", "Ferret", "Other"
)

route_options <- c(
  "IV bolus" = "IV",
  "IV infusion" = "IV_INF",
  "Oral (PO)" = "PO",
  "Subcutaneous (SC)" = "SC",
  "Intramuscular (IM)" = "IM",
  "Intraperitoneal (IP)" = "IP",
  "Other extravascular" = "EV"
)

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  # Custom CSS
  tags$head(tags$style(HTML("
    body { background-color: #f7f9fc; }
    .main-header {
      background: linear-gradient(135deg, #1a3a5c 0%, #2d6a9f 100%);
      color: white; padding: 20px 30px; margin-bottom: 20px;
      border-radius: 0 0 8px 8px;
    }
    .main-header h2 { margin: 0 0 5px 0; font-weight: 700; }
    .main-header p { margin: 0; opacity: 0.85; font-size: 14px; }
    .well { background-color: #ffffff; border: 1px solid #dce3ec; border-radius: 6px; }
    .nav-tabs > li.active > a { border-top: 3px solid #2d6a9f; font-weight: 600; }
    .btn-primary { background-color: #2d6a9f; border-color: #1a3a5c; }
    .btn-primary:hover { background-color: #1a3a5c; }
    .section-title { color: #1a3a5c; border-bottom: 2px solid #2d6a9f;
                     padding-bottom: 8px; margin-bottom: 15px; font-weight: 600; }
    .status-box {
      padding: 12px 15px; border-radius: 6px; margin: 10px 0;
      border-left: 4px solid;
    }
    .status-info { background-color: #e8f0fe; border-color: #2d6a9f; color: #1a3a5c; }
    .status-success { background-color: #e6f4ea; border-color: #34a853; color: #1e4620; }
    .status-error { background-color: #fce8e6; border-color: #ea4335; color: #5f1412; }
    .status-warning { background-color: #fef7e0; border-color: #fbbc04; color: #5f4b08; }
  "))),

  # Header
  div(class = "main-header",
    h2("PopPK Analysis"),
    p("Population Pharmacokinetic Analysis - Non-Compartmental & Compartmental Modeling")
  ),

  sidebarLayout(
    # ---- Sidebar ----
    sidebarPanel(
      width = 3,

      h4("Data Upload", class = "section-title"),
      fileInput("data_file", "Upload PK Data",
                accept = c(".csv", ".xls", ".xlsx"),
                placeholder = "CSV or Excel file"),
      helpText("Required columns: ID, Time, Concentration"),

      hr(),
      h4("Study Information", class = "section-title"),

      textInput("drug_name", "Drug/Compound Name", placeholder = "e.g., Amoxicillin"),

      selectInput("species", "Species", choices = species_options, selected = "Dog"),

      selectInput("route", "Route of Administration", choices = route_options, selected = "IV"),

      numericInput("dose", "Dose", value = 10, min = 0.001, step = 0.1),

      fluidRow(
        column(6, selectInput("dose_unit", "Dose Unit",
                              choices = c("mg", "mg/kg", "ug", "ug/kg", "g", "g/kg", "nmol", "umol"),
                              selected = "mg/kg")),
        column(6, selectInput("conc_unit", "Conc Unit",
                              choices = c("ng/mL", "ug/mL", "mg/L", "ug/L", "nmol/L", "umol/L"),
                              selected = "ug/mL"))
      ),

      selectInput("time_unit", "Time Unit",
                  choices = c("h" = "h", "min" = "min", "days" = "days"),
                  selected = "h"),

      hr(),
      h4("Analysis Settings", class = "section-title"),

      selectInput("model_type", "PK Model",
                  choices = c(
                    "Non-Compartmental Analysis (NCA)" = "NCA",
                    "One-Compartment Model" = "1comp",
                    "Two-Compartment Model" = "2comp",
                    "Three-Compartment Model" = "3comp"
                  ),
                  selected = "NCA"),

      conditionalPanel(
        condition = "input.model_type != 'NCA'",
        helpText("Compartmental models use Non-Linear Mixed Effects (nlme) modeling with
                  between-subject random effects on key parameters.")
      ),

      hr(),
      actionButton("run_analysis", "Run Analysis",
                    class = "btn-primary btn-block",
                    icon = icon("play")),

      br(), br(),

      # Download buttons
      conditionalPanel(
        condition = "output.analysis_done",
        hr(),
        h4("Export Results", class = "section-title"),
        downloadButton("download_params", "Download Parameters (CSV)", class = "btn-block"),
        br(),
        downloadButton("download_plot", "Download Plot (PNG)", class = "btn-block")
      )
    ),

    # ---- Main Panel ----
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "main_tabs",

        # Tab 1: Data Preview
        tabPanel("Data",
          icon = icon("table"),
          br(),
          uiOutput("data_status"),
          conditionalPanel(
            condition = "output.data_loaded",
            fluidRow(
              column(4, wellPanel(
                h5("Dataset Summary"),
                uiOutput("data_summary_info")
              )),
              column(8, wellPanel(
                h5("Data Preview"),
                DTOutput("data_preview")
              ))
            )
          )
        ),

        # Tab 2: Concentration-Time Plots
        tabPanel("Conc-Time Plots",
          icon = icon("chart-line"),
          br(),
          conditionalPanel(
            condition = "output.data_loaded",
            wellPanel(
              fluidRow(
                column(4, checkboxInput("log_y_obs", "Log-transformed Y-axis", value = FALSE)),
                column(4, radioButtons("plot_type_obs", "Plot Type",
                                        choices = c("Mean \u00B1 SEM" = "mean",
                                                    "Individual" = "individual",
                                                    "Both" = "both"),
                                        inline = TRUE, selected = "mean"))
              )
            ),
            uiOutput("obs_plots_ui")
          )
        ),

        # Tab 3: Analysis Results
        tabPanel("PK Parameters",
          icon = icon("flask"),
          br(),
          uiOutput("analysis_status"),
          conditionalPanel(
            condition = "output.analysis_done",
            uiOutput("results_ui")
          )
        ),

        # Tab 4: Model Fit (compartmental only)
        tabPanel("Model Fit",
          icon = icon("project-diagram"),
          br(),
          conditionalPanel(
            condition = "output.is_compartmental",
            wellPanel(
              fluidRow(
                column(4, checkboxInput("log_y_fit", "Log-transformed Y-axis", value = FALSE)),
                column(4, checkboxInput("show_diagnostics", "Show Diagnostic Plots", value = TRUE))
              )
            ),
            plotOutput("model_fit_plot", height = "500px"),
            conditionalPanel(
              condition = "input.show_diagnostics",
              hr(),
              h4("Goodness-of-Fit Diagnostics", class = "section-title"),
              fluidRow(
                column(6, plotOutput("diag_obs_pred", height = "350px")),
                column(6, plotOutput("diag_resid", height = "350px"))
              ),
              fluidRow(
                column(6, plotOutput("diag_qq", height = "350px")),
                column(6, plotOutput("diag_hist", height = "350px"))
              )
            )
          ),
          conditionalPanel(
            condition = "!output.is_compartmental",
            div(class = "status-box status-info",
                "Model fit plots are available for compartmental models (1-, 2-, 3-compartment).
                 Select a compartmental model and run the analysis to view fits.")
          )
        )
      )
    )
  )
)

# ============================================================
# Server
# ============================================================
server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(
    pk_data = NULL,
    data_summary = NULL,
    nca_results = NULL,
    nca_summary_results = NULL,
    comp_results = NULL,
    analysis_complete = FALSE,
    analysis_type = NULL,
    error_msg = NULL
  )

  # ---- Data Loading ----
  observeEvent(input$data_file, {
    req(input$data_file)
    result <- load_pk_data(input$data_file$datapath, input$data_file$name)

    if (!is.null(result$error)) {
      rv$pk_data <- NULL
      rv$data_summary <- NULL
      rv$error_msg <- result$error
    } else {
      rv$pk_data <- result$data
      rv$data_summary <- pk_data_summary(result$data)
      rv$error_msg <- NULL
      rv$analysis_complete <- FALSE

      showNotification(
        paste("Data loaded:", rv$data_summary$n_subjects, "subjects,",
              rv$data_summary$n_observations, "observations"),
        type = "message", duration = 5
      )
    }
  })

  # ---- Output flags for conditionalPanel ----
  output$data_loaded <- reactive({ !is.null(rv$pk_data) })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

  output$analysis_done <- reactive({ rv$analysis_complete })
  outputOptions(output, "analysis_done", suspendWhenHidden = FALSE)

  output$is_compartmental <- reactive({
    rv$analysis_complete && !is.null(rv$analysis_type) && rv$analysis_type != "NCA"
  })
  outputOptions(output, "is_compartmental", suspendWhenHidden = FALSE)

  # ---- Data Status ----
  output$data_status <- renderUI({
    if (!is.null(rv$error_msg)) {
      div(class = "status-box status-error", icon("exclamation-triangle"), rv$error_msg)
    } else if (is.null(rv$pk_data)) {
      div(class = "status-box status-info", icon("info-circle"),
          "Upload a CSV or Excel file with columns: ID, Time, Concentration")
    } else {
      div(class = "status-box status-success", icon("check-circle"),
          paste("Data loaded successfully:", rv$data_summary$n_subjects, "subjects,",
                rv$data_summary$n_observations, "observations"))
    }
  })

  # ---- Data Summary ----
  output$data_summary_info <- renderUI({
    req(rv$data_summary)
    s <- rv$data_summary
    tags$div(
      tags$p(tags$strong("Subjects: "), s$n_subjects),
      tags$p(tags$strong("Observations: "), s$n_observations),
      tags$p(tags$strong("Time range: "), round(s$time_range[1], 2), " - ", round(s$time_range[2], 2)),
      tags$p(tags$strong("Conc range: "), round(s$conc_range[1], 3), " - ", round(s$conc_range[2], 3)),
      tags$p(tags$strong("Subject IDs: "), paste(s$subjects, collapse = ", "))
    )
  })

  # ---- Data Preview Table ----
  output$data_preview <- renderDT({
    req(rv$pk_data)
    datatable(
      rv$pk_data[, c("ID", "Time", "Conc")],
      options = list(pageLength = 15, scrollX = TRUE),
      colnames = c("Subject ID", "Time", "Concentration"),
      rownames = FALSE
    ) %>%
      formatRound(columns = c("Time", "Conc"), digits = 4)
  })

  # ---- Observational Plots ----
  output$obs_plots_ui <- renderUI({
    req(rv$pk_data)
    plot_type <- input$plot_type_obs

    if (plot_type == "mean") {
      plotOutput("mean_plot", height = "500px")
    } else if (plot_type == "individual") {
      plotOutput("indiv_plot", height = "500px")
    } else {
      tagList(
        plotOutput("mean_plot", height = "450px"),
        br(),
        plotOutput("indiv_plot", height = "450px")
      )
    }
  })

  output$mean_plot <- renderPlot({
    req(rv$pk_data)
    plot_mean_conc_time(
      rv$pk_data,
      log_y = input$log_y_obs,
      species = input$species,
      drug_name = input$drug_name,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit
    )
  })

  output$indiv_plot <- renderPlot({
    req(rv$pk_data)
    plot_individual_conc_time(
      rv$pk_data,
      log_y = input$log_y_obs,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit
    )
  })

  # ---- Run Analysis ----
  observeEvent(input$run_analysis, {
    req(rv$pk_data)

    # Validate inputs
    if (is.na(input$dose) || input$dose <= 0) {
      showNotification("Please enter a valid dose > 0", type = "error")
      return()
    }

    # Determine if route is IV (for NCA and model selection)
    route_val <- input$route
    is_iv <- route_val %in% c("IV", "IV_INF")
    route_for_model <- if (is_iv) "IV" else "EV"

    rv$analysis_type <- input$model_type
    rv$analysis_complete <- FALSE
    rv$error_msg <- NULL

    # Show progress
    withProgress(message = "Running PK analysis...", value = 0, {
      if (input$model_type == "NCA") {
        # ---- Non-Compartmental Analysis ----
        incProgress(0.3, detail = "Computing NCA parameters...")

        rv$nca_results <- run_nca(rv$pk_data, input$dose, route_for_model)
        rv$nca_summary_results <- nca_summary(rv$nca_results)
        rv$comp_results <- NULL

        incProgress(0.7, detail = "Done!")
        rv$analysis_complete <- TRUE

        showNotification("NCA completed successfully", type = "message", duration = 5)

      } else {
        # ---- Compartmental Modeling ----
        incProgress(0.2, detail = "Estimating initial parameters...")

        result <- fit_compartmental_model(
          rv$pk_data, input$dose, input$model_type, route_for_model
        )

        if (!is.null(result$error)) {
          rv$error_msg <- result$error
          rv$comp_results <- NULL
          showNotification(result$error, type = "error", duration = 10)
        } else {
          rv$comp_results <- result
          rv$nca_results <- NULL
          rv$nca_summary_results <- NULL
          rv$analysis_complete <- TRUE

          incProgress(0.8, detail = "Done!")
          showNotification(
            paste(input$model_type, "model fitted successfully"),
            type = "message", duration = 5
          )
        }
      }
    })

    # Switch to results tab
    if (rv$analysis_complete) {
      updateTabsetPanel(session, "main_tabs", selected = "PK Parameters")
    }
  })

  # ---- Analysis Status ----
  output$analysis_status <- renderUI({
    if (!is.null(rv$error_msg) && !rv$analysis_complete) {
      div(class = "status-box status-error", icon("exclamation-triangle"), rv$error_msg)
    } else if (!rv$analysis_complete) {
      div(class = "status-box status-info", icon("info-circle"),
          "Upload data, configure settings, and click 'Run Analysis' to compute PK parameters.")
    } else {
      model_desc <- switch(rv$analysis_type,
        "NCA" = "Non-Compartmental Analysis",
        "1comp" = "One-Compartment Model (NLME)",
        "2comp" = "Two-Compartment Model (NLME)",
        "3comp" = "Three-Compartment Model (NLME)"
      )
      div(class = "status-box status-success", icon("check-circle"),
          paste("Analysis complete:", model_desc))
    }
  })

  # ---- Results UI ----
  output$results_ui <- renderUI({
    req(rv$analysis_complete)

    if (rv$analysis_type == "NCA") {
      # NCA Results
      tagList(
        h4("Individual NCA Parameters", class = "section-title"),
        DTOutput("nca_individual_table"),
        br(),
        h4("Summary Statistics", class = "section-title"),
        DTOutput("nca_summary_table")
      )
    } else {
      # Compartmental Results
      tagList(
        h4("Population (Fixed Effect) Parameters", class = "section-title"),
        DTOutput("comp_pop_table"),
        br(),
        h4("Individual Parameter Estimates", class = "section-title"),
        DTOutput("comp_indiv_table"),
        br(),
        wellPanel(
          h5("Model Diagnostics"),
          uiOutput("model_info")
        )
      )
    }
  })

  # ---- NCA Tables ----
  output$nca_individual_table <- renderDT({
    req(rv$nca_results)
    df <- rv$nca_results
    # Format numeric columns
    num_cols <- setdiff(names(df), "ID")
    datatable(df, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatSignif(columns = num_cols, digits = 4)
  })

  output$nca_summary_table <- renderDT({
    req(rv$nca_summary_results)
    datatable(rv$nca_summary_results, rownames = FALSE,
              options = list(scrollX = TRUE, dom = 't')) %>%
      formatSignif(columns = c("Mean", "SD", "SEM", "Median", "Min", "Max"), digits = 4)
  })

  # ---- Compartmental Tables ----
  output$comp_pop_table <- renderDT({
    req(rv$comp_results)
    df <- rv$comp_results$summary
    # Remove diagnostic columns for display
    display_df <- df[, c("Parameter", "Estimate")]
    display_df$Estimate <- signif(display_df$Estimate, 4)
    datatable(display_df, rownames = FALSE,
              options = list(dom = 't', scrollX = TRUE))
  })

  output$comp_indiv_table <- renderDT({
    req(rv$comp_results)
    df <- rv$comp_results$params
    num_cols <- setdiff(names(df), "ID")
    datatable(df, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatSignif(columns = num_cols, digits = 4)
  })

  output$model_info <- renderUI({
    req(rv$comp_results)
    pop <- rv$comp_results$summary
    aic_val <- pop$AIC[1]
    bic_val <- pop$BIC[1]
    ll_val <- pop$logLik[1]

    tags$div(
      tags$p(tags$strong("AIC: "), round(aic_val, 2)),
      tags$p(tags$strong("BIC: "), round(bic_val, 2)),
      tags$p(tags$strong("Log-Likelihood: "), round(ll_val, 2))
    )
  })

  # ---- Model Fit Plots ----
  output$model_fit_plot <- renderPlot({
    req(rv$comp_results, rv$comp_results$predictions)
    model_labels <- c(
      "1comp" = "One-Compartment",
      "2comp" = "Two-Compartment",
      "3comp" = "Three-Compartment"
    )
    route_label <- names(route_options)[route_options == input$route]

    plot_model_fit(
      rv$pk_data,
      rv$comp_results$predictions,
      log_y = input$log_y_fit,
      model_label = paste(model_labels[rv$analysis_type], "-", route_label),
      time_unit = input$time_unit,
      conc_unit = input$conc_unit
    )
  })

  # ---- Diagnostic Plots ----
  output$diag_obs_pred <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, rv$pk_data)
    plots$obs_vs_pred
  })

  output$diag_resid <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, rv$pk_data)
    plots$resid_vs_pred
  })

  output$diag_qq <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, rv$pk_data)
    plots$qq_plot
  })

  output$diag_hist <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, rv$pk_data)
    plots$resid_hist
  })

  # ---- Downloads ----
  output$download_params <- downloadHandler(
    filename = function() {
      paste0("pk_parameters_", rv$analysis_type, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (rv$analysis_type == "NCA") {
        # Write both individual and summary
        write.csv(rv$nca_results, file, row.names = FALSE)
      } else {
        # Write population + individual params
        pop <- rv$comp_results$summary
        indiv <- rv$comp_results$params
        write.csv(pop, file, row.names = FALSE)
      }
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() {
      paste0("pk_conc_time_", Sys.Date(), ".png")
    },
    content = function(file) {
      p <- plot_mean_conc_time(
        rv$pk_data,
        log_y = input$log_y_obs,
        species = input$species,
        drug_name = input$drug_name,
        time_unit = input$time_unit,
        conc_unit = input$conc_unit
      )
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
