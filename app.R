# Population PK Modeling Application
# Author: Hugues Beaufrere, DVM, PhD, DACZM

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
source("R/manuscript_writer.R")

# ============================================================
# Species options
# ============================================================
species_options <- c(
  "Orange-winged Amazon parrot", "Cockatiel", "Great horned owl",
  "Red-tailed hawk", "Rabbit", "Bearded dragon",
  "Other"
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
    .main-header .author { margin: 0; opacity: 0.75; font-size: 12px; font-style: italic; }
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
    .subject-selector .checkbox { margin-top: 2px; margin-bottom: 2px; }
    .manuscript-output {
      background-color: #fff; border: 1px solid #dce3ec; border-radius: 6px;
      padding: 25px 30px; margin-top: 15px; font-family: 'Georgia', 'Times New Roman', serif;
      font-size: 14px; line-height: 1.8; color: #222;
    }
    .manuscript-output h4 { color: #1a3a5c; font-weight: 700; margin-top: 20px;
                            margin-bottom: 10px; font-family: Arial, sans-serif; }
    .manuscript-section { margin-bottom: 20px; }
    .btn-manuscript {
      background: linear-gradient(135deg, #5b2d8e 0%, #7b3fa0 100%);
      color: white; border: none; padding: 12px 24px; font-size: 15px;
      font-weight: 600; border-radius: 6px; cursor: pointer;
    }
    .btn-manuscript:hover { background: linear-gradient(135deg, #4a2475 0%, #6a348d 100%); color: white; }
    .api-key-input { max-width: 500px; }
    .compound-section { border-top: 2px solid #2d6a9f; padding-top: 15px; margin-top: 25px; }
  "))),

  # Header
  div(class = "main-header",
    h2("Population PK Modeling"),
    tags$div(class = "author", "Hugues Beaufr\u00e8re, DVM, PhD, DACZM"),
    p("Non-Compartmental & Compartmental Modeling")
  ),

  sidebarLayout(
    # ---- Sidebar ----
    sidebarPanel(
      width = 3,

      h4("Data Upload", class = "section-title"),
      fileInput("data_file", "Upload PK Data",
                accept = c(".csv", ".xls", ".xlsx"),
                placeholder = "CSV or Excel file"),
      helpText("Format: first column = Time, each subsequent column = concentrations for one animal.",
               "Optionally include a 'Drug' column for multiple compounds."),
      actionButton("load_example", "Load Example Dataset",
                    class = "btn-default btn-block btn-sm",
                    icon = icon("flask")),
      helpText(style = "font-size: 11px; color: #777;",
        "Acetaminophen PO in Orange-winged Amazon parrots (n=8)"),
      actionButton("load_example_multi", "Load Example Dataset (Multiple Compounds)",
                    class = "btn-default btn-block btn-sm",
                    icon = icon("vials")),
      helpText(style = "font-size: 11px; color: #777;",
        "APAP, APAP+NAC, APAP+Silymarin PO in Amazon parrots (n=8 per group)"),

      hr(),
      h4("Study Information", class = "section-title"),

      textInput("drug_name", "Drug/Compound Name", placeholder = "e.g., Amoxicillin"),

      selectInput("species", "Species", choices = species_options, selected = "Orange-winged Amazon parrot"),

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

      fluidRow(
        column(6, numericInput("body_weight", "Body Weight (optional)",
                               value = NA, min = 0.001, step = 0.1)),
        column(6, selectInput("weight_unit", "Weight Unit",
                              choices = c("kg", "g"),
                              selected = "kg"))
      ),
      helpText("Enter mean body weight if dose is absolute (e.g., mg) and weight-normalized PK parameters are desired."),

      # Subject exclusion (appears when data is loaded)
      conditionalPanel(
        condition = "output.data_loaded",
        hr(),
        h4("Subject Selection", class = "section-title"),
        helpText("Uncheck subjects to exclude from plots and analysis."),
        div(class = "subject-selector",
          uiOutput("subject_checkboxes")
        )
      ),

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
        downloadButton("download_params", "Download Parameters (CSV)", class = "btn-block")
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
                column(3, radioButtons("plot_type_obs", "Plot Type",
                                        choices = c("Mean \u00B1 SEM" = "mean",
                                                    "Individual" = "individual",
                                                    "Both" = "both"),
                                        inline = FALSE, selected = "mean")),
                column(3,
                  checkboxInput("log_y_obs", "Log-transformed Y-axis", value = FALSE),
                  checkboxInput("show_ci", "Show confidence shading", value = TRUE),
                  checkboxInput("show_grid", "Show background grid", value = TRUE),
                  conditionalPanel(
                    condition = "output.has_multi_compound",
                    checkboxInput("superimpose_compounds", "Superimpose compounds", value = FALSE)
                  )
                ),
                column(3,
                  radioButtons("color_mode", "Color Scheme",
                               choices = c("Color" = "color",
                                           "Color (colorblind-friendly)" = "cb",
                                           "BW" = "bw"),
                               inline = FALSE, selected = "color"),
                  sliderInput("point_size", "Point Size", min = 1, max = 8, value = 2.5, step = 0.5),
                  sliderInput("line_width", "Line Width", min = 0.5, max = 4, value = 1, step = 0.25),
                  sliderInput("errorbar_width", "Error Bar Width", min = 0.3, max = 3, value = 0.5, step = 0.1)
                ),
                column(3,
                  downloadButton("download_mean_plot", "Download Mean Plot", class = "btn-block btn-sm"),
                  br(),
                  downloadButton("download_indiv_plot", "Download Individual Plot", class = "btn-block btn-sm")
                )
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
            uiOutput("model_fit_ui")
          ),
          conditionalPanel(
            condition = "!output.is_compartmental",
            div(class = "status-box status-info",
                "Model fit plots are available for compartmental models (1-, 2-, 3-compartment).
                 Select a compartmental model and run the analysis to view fits.")
          )
        ),

        # Tab 5: For Manuscript
        tabPanel("For Manuscript",
          icon = icon("file-alt"),
          br(),
          conditionalPanel(
            condition = "output.analysis_done",
            wellPanel(
              h4("AI-Powered Manuscript Section Writing", class = "section-title"),
              p("Generate publication-ready Materials & Methods, Results, and Interpretation sections ",
                "based on your analysis. Powered by Claude (Anthropic API)."),
              div(class = "api-key-input",
                passwordInput("claude_api_key", "Anthropic API Key",
                              placeholder = "sk-ant-..."),
                helpText("Your API key is used only for this request and is not stored. ",
                         "Get a key at ", tags$a("console.anthropic.com",
                         href = "https://console.anthropic.com/", target = "_blank"), ".")
              ),
              actionButton("generate_manuscript", "AI-Powered Manuscript Section Writing",
                            class = "btn-manuscript",
                            icon = icon("magic"))
            ),
            uiOutput("manuscript_output")
          ),
          conditionalPanel(
            condition = "!output.analysis_done",
            div(class = "status-box status-info",
                icon("info-circle"),
                "Run an analysis first (NCA or compartmental) to enable manuscript generation.")
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
    comp_summary_stats = NULL,
    analysis_complete = FALSE,
    analysis_type = NULL,
    error_msg = NULL,
    # Multi-compound analysis results (list keyed by drug name)
    multi_nca_results = NULL,
    multi_nca_summary = NULL,
    multi_comp_results = NULL,
    multi_comp_summary = NULL
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

      desc <- paste("Data loaded:", rv$data_summary$n_subjects, "subjects,",
                    rv$data_summary$n_observations, "observations")
      if (isTRUE(rv$data_summary$has_drug)) {
        desc <- paste0(desc, " (", rv$data_summary$n_drugs, " compounds)")
      }
      showNotification(desc, type = "message", duration = 5)
    }
  })

  # ---- Load Example Dataset (single compound) ----
  observeEvent(input$load_example, {
    example_path <- file.path("data", "example_acetaminophen_amazons.csv")
    if (!file.exists(example_path)) {
      showNotification("Example dataset file not found.", type = "error")
      return()
    }

    result <- load_pk_data(example_path, "example_acetaminophen_amazons.csv")

    if (!is.null(result$error)) {
      showNotification(result$error, type = "error")
      return()
    }

    rv$pk_data <- result$data
    rv$data_summary <- pk_data_summary(result$data)
    rv$error_msg <- NULL
    rv$analysis_complete <- FALSE

    # Pre-fill study information
    updateTextInput(session, "drug_name", value = "Acetaminophen")
    updateSelectInput(session, "species", selected = "Orange-winged Amazon parrot")
    updateSelectInput(session, "route", selected = "PO")
    updateNumericInput(session, "dose", value = 100)
    updateSelectInput(session, "dose_unit", selected = "mg/kg")
    updateSelectInput(session, "conc_unit", selected = "ng/mL")
    updateSelectInput(session, "time_unit", selected = "h")

    showNotification(
      paste("Example dataset loaded:", rv$data_summary$n_subjects, "birds,",
            rv$data_summary$n_observations, "observations"),
      type = "message", duration = 5
    )
  })

  # ---- Load Example Dataset (multiple compounds) ----
  observeEvent(input$load_example_multi, {
    example_path <- file.path("data", "example_multicompound_apap.csv")
    if (!file.exists(example_path)) {
      showNotification("Multi-compound example dataset file not found.", type = "error")
      return()
    }

    result <- load_pk_data(example_path, "example_multicompound_apap.csv")

    if (!is.null(result$error)) {
      showNotification(result$error, type = "error")
      return()
    }

    rv$pk_data <- result$data
    rv$data_summary <- pk_data_summary(result$data)
    rv$error_msg <- NULL
    rv$analysis_complete <- FALSE

    # Pre-fill study information
    updateTextInput(session, "drug_name", value = "Acetaminophen")
    updateSelectInput(session, "species", selected = "Orange-winged Amazon parrot")
    updateSelectInput(session, "route", selected = "PO")
    updateNumericInput(session, "dose", value = 100)
    updateSelectInput(session, "dose_unit", selected = "mg/kg")
    updateSelectInput(session, "conc_unit", selected = "ng/mL")
    updateSelectInput(session, "time_unit", selected = "h")

    showNotification(
      paste("Multi-compound example loaded:", rv$data_summary$n_subjects, "subjects,",
            rv$data_summary$n_drugs, "compounds,",
            rv$data_summary$n_observations, "observations"),
      type = "message", duration = 5
    )
  })

  # ---- Subject Selection Checkboxes ----
  output$subject_checkboxes <- renderUI({
    req(rv$data_summary)
    subjects <- rv$data_summary$subjects
    checkboxGroupInput("included_subjects", label = NULL,
                       choices = subjects, selected = subjects)
  })

  # ---- Filtered data (respects subject exclusion) ----
  filtered_data <- reactive({
    req(rv$pk_data)
    included <- input$included_subjects
    if (is.null(included) || length(included) == 0) {
      return(rv$pk_data[0, ])  # empty data frame
    }
    rv$pk_data[rv$pk_data$ID %in% included, ]
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

  output$has_multi_compound <- reactive({
    !is.null(rv$data_summary) && isTRUE(rv$data_summary$has_drug)
  })
  outputOptions(output, "has_multi_compound", suspendWhenHidden = FALSE)

  # ---- Data Status ----
  output$data_status <- renderUI({
    if (!is.null(rv$error_msg)) {
      div(class = "status-box status-error", icon("exclamation-triangle"), rv$error_msg)
    } else if (is.null(rv$pk_data)) {
      div(class = "status-box status-info", icon("info-circle"),
          "Upload a CSV or Excel file. First column = Time, each subsequent column = concentrations for one animal.",
          "Optionally include a 'Drug' column for multiple compounds.")
    } else {
      n_included <- length(input$included_subjects)
      n_total <- rv$data_summary$n_subjects
      excl_text <- if (n_included < n_total) paste0(" (", n_total - n_included, " excluded)") else ""
      drug_text <- ""
      if (isTRUE(rv$data_summary$has_drug)) {
        drug_text <- paste0(", ", rv$data_summary$n_drugs, " compounds: ",
                           paste(rv$data_summary$drugs, collapse = ", "))
      }
      div(class = "status-box status-success", icon("check-circle"),
          paste0("Data loaded: ", n_included, " of ", n_total, " subjects included", excl_text,
                 ", ", nrow(filtered_data()), " observations", drug_text))
    }
  })

  # ---- Data Summary ----
  output$data_summary_info <- renderUI({
    req(rv$data_summary)
    fdata <- filtered_data()
    s <- rv$data_summary
    n_included <- length(input$included_subjects)
    items <- tagList(
      tags$p(tags$strong("Total Subjects: "), s$n_subjects),
      tags$p(tags$strong("Included: "), n_included),
      tags$p(tags$strong("Observations: "), nrow(fdata)),
      tags$p(tags$strong("Time range: "), round(s$time_range[1], 2), " - ", round(s$time_range[2], 2))
    )
    if (nrow(fdata) > 0) {
      items <- tagList(items,
        tags$p(tags$strong("Conc range: "),
               round(min(fdata$Conc, na.rm = TRUE), 3), " - ",
               round(max(fdata$Conc, na.rm = TRUE), 3))
      )
    }
    if (isTRUE(s$has_drug)) {
      items <- tagList(items,
        tags$p(tags$strong("Compounds: "), paste(s$drugs, collapse = ", "))
      )
    }
    items <- tagList(items,
      tags$p(tags$strong("Subject IDs: "), paste(s$subjects, collapse = ", "))
    )
    items
  })

  # ---- Data Preview Table ----
  output$data_preview <- renderDT({
    req(rv$pk_data)
    has_drug <- "Drug" %in% names(rv$pk_data)
    if (has_drug) {
      display_cols <- c("ID", "Time", "Conc", "Drug")
      col_names <- c("Subject ID", "Time", "Concentration", "Compound")
    } else {
      display_cols <- c("ID", "Time", "Conc")
      col_names <- c("Subject ID", "Time", "Concentration")
    }
    datatable(
      rv$pk_data[, display_cols],
      options = list(pageLength = 15, scrollX = TRUE),
      colnames = col_names,
      rownames = FALSE
    ) %>%
      formatRound(columns = c("Time", "Conc"), digits = 4)
  })

  # ---- Observational Plots ----
  output$obs_plots_ui <- renderUI({
    req(filtered_data())
    plot_type <- input$plot_type_obs
    is_multi <- isTRUE(rv$data_summary$has_drug)
    superimpose <- isTRUE(input$superimpose_compounds) && is_multi

    if (superimpose) {
      # Superimposed multi-compound plots
      if (plot_type == "mean") {
        plotOutput("mean_plot_multi", height = "550px")
      } else if (plot_type == "individual") {
        plotOutput("indiv_plot_multi", height = "550px")
      } else {
        tagList(
          plotOutput("mean_plot_multi", height = "500px"),
          br(),
          plotOutput("indiv_plot_multi", height = "500px")
        )
      }
    } else if (is_multi && !superimpose) {
      # Per-compound plots stacked
      drugs <- levels(filtered_data()$Drug)
      if (is.null(drugs)) drugs <- unique(as.character(filtered_data()$Drug))
      plot_list <- tagList()
      for (i in seq_along(drugs)) {
        drug <- drugs[i]
        mean_id <- paste0("mean_plot_drug_", i)
        indiv_id <- paste0("indiv_plot_drug_", i)
        if (plot_type == "mean") {
          plot_list <- tagList(plot_list,
            h4(drug, class = "section-title"),
            plotOutput(mean_id, height = "450px"),
            br()
          )
        } else if (plot_type == "individual") {
          plot_list <- tagList(plot_list,
            h4(drug, class = "section-title"),
            plotOutput(indiv_id, height = "450px"),
            br()
          )
        } else {
          plot_list <- tagList(plot_list,
            h4(drug, class = "section-title"),
            plotOutput(mean_id, height = "400px"),
            br(),
            plotOutput(indiv_id, height = "400px"),
            br()
          )
        }
      }
      plot_list
    } else {
      # Single compound
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
    }
  })

  # Helper to build mean plot (reused for display and download) - single compound
  build_mean_plot <- reactive({
    req(filtered_data(), nrow(filtered_data()) > 0)
    plot_mean_conc_time(
      filtered_data(),
      log_y = input$log_y_obs,
      species = input$species,
      drug_name = input$drug_name,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit,
      color_mode = input$color_mode,
      show_ci = input$show_ci,
      point_size = input$point_size,
      line_width = input$line_width,
      errorbar_width = input$errorbar_width,
      show_grid = input$show_grid
    )
  })

  # Helper to build individual plot (reused for display and download) - single compound
  build_indiv_plot <- reactive({
    req(filtered_data(), nrow(filtered_data()) > 0)
    plot_individual_conc_time(
      filtered_data(),
      log_y = input$log_y_obs,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit,
      color_mode = input$color_mode,
      point_size = input$point_size,
      line_width = input$line_width,
      show_grid = input$show_grid
    )
  })

  # Multi-compound superimposed plots
  build_mean_plot_multi <- reactive({
    req(filtered_data(), nrow(filtered_data()) > 0)
    req("Drug" %in% names(filtered_data()))
    plot_mean_conc_time_by_drug(
      filtered_data(),
      log_y = input$log_y_obs,
      species = input$species,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit,
      color_mode = input$color_mode,
      show_ci = input$show_ci,
      point_size = input$point_size,
      line_width = input$line_width,
      errorbar_width = input$errorbar_width,
      show_grid = input$show_grid
    )
  })

  build_indiv_plot_multi <- reactive({
    req(filtered_data(), nrow(filtered_data()) > 0)
    req("Drug" %in% names(filtered_data()))
    plot_individual_conc_time_by_drug(
      filtered_data(),
      log_y = input$log_y_obs,
      time_unit = input$time_unit,
      conc_unit = input$conc_unit,
      color_mode = input$color_mode,
      point_size = input$point_size,
      line_width = input$line_width,
      show_grid = input$show_grid
    )
  })

  output$mean_plot <- renderPlot({ build_mean_plot() })
  output$indiv_plot <- renderPlot({ build_indiv_plot() })
  output$mean_plot_multi <- renderPlot({ build_mean_plot_multi() })
  output$indiv_plot_multi <- renderPlot({ build_indiv_plot_multi() })

  # Per-compound individual plots (dynamic rendering)
  observe({
    fdata <- filtered_data()
    req(fdata, nrow(fdata) > 0)
    if (!("Drug" %in% names(fdata))) return()

    drugs <- levels(fdata$Drug)
    if (is.null(drugs)) drugs <- unique(as.character(fdata$Drug))

    for (i in seq_along(drugs)) {
      local({
        drug_idx <- i
        drug_name_local <- drugs[drug_idx]

        output[[paste0("mean_plot_drug_", drug_idx)]] <- renderPlot({
          req(filtered_data())
          drug_data <- filtered_data()[filtered_data()$Drug == drug_name_local, ]
          req(nrow(drug_data) > 0)
          plot_mean_conc_time(
            drug_data,
            log_y = input$log_y_obs,
            species = input$species,
            drug_name = drug_name_local,
            time_unit = input$time_unit,
            conc_unit = input$conc_unit,
            color_mode = input$color_mode,
            show_ci = input$show_ci,
            point_size = input$point_size,
            line_width = input$line_width,
            errorbar_width = input$errorbar_width,
            show_grid = input$show_grid
          )
        })

        output[[paste0("indiv_plot_drug_", drug_idx)]] <- renderPlot({
          req(filtered_data())
          drug_data <- filtered_data()[filtered_data()$Drug == drug_name_local, ]
          req(nrow(drug_data) > 0)
          plot_individual_conc_time(
            drug_data,
            log_y = input$log_y_obs,
            time_unit = input$time_unit,
            conc_unit = input$conc_unit,
            color_mode = input$color_mode,
            point_size = input$point_size,
            line_width = input$line_width,
            show_grid = input$show_grid
          )
        })
      })
    }
  })

  # ---- Plot Downloads ----
  output$download_mean_plot <- downloadHandler(
    filename = function() {
      paste0("pk_mean_conc_time_", Sys.Date(), ".png")
    },
    content = function(file) {
      is_multi <- isTRUE(rv$data_summary$has_drug)
      superimpose <- isTRUE(input$superimpose_compounds) && is_multi
      if (superimpose) {
        p <- build_mean_plot_multi()
      } else {
        p <- build_mean_plot()
      }
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )

  output$download_indiv_plot <- downloadHandler(
    filename = function() {
      paste0("pk_individual_conc_time_", Sys.Date(), ".png")
    },
    content = function(file) {
      is_multi <- isTRUE(rv$data_summary$has_drug)
      superimpose <- isTRUE(input$superimpose_compounds) && is_multi
      if (superimpose) {
        p <- build_indiv_plot_multi()
      } else {
        p <- build_indiv_plot()
      }
      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )

  # ---- Run Analysis ----
  observeEvent(input$run_analysis, {
    fdata <- filtered_data()
    req(fdata, nrow(fdata) > 0)

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

    # Check if multi-compound
    is_multi <- "Drug" %in% names(fdata) && length(unique(fdata$Drug)) > 1

    # Show progress
    withProgress(message = "Running PK analysis...", value = 0, {
      if (is_multi) {
        # ---- Multi-compound analysis ----
        drugs <- unique(as.character(fdata$Drug))

        if (input$model_type == "NCA") {
          multi_nca <- list()
          multi_nca_sum <- list()

          for (di in seq_along(drugs)) {
            drug <- drugs[di]
            incProgress(0.8 / length(drugs), detail = paste("NCA for", drug, "..."))
            drug_data <- fdata[fdata$Drug == drug, ]

            nca_res <- run_nca(drug_data, input$dose, route_for_model)
            nca_sum <- nca_summary(nca_res)
            multi_nca[[drug]] <- nca_res
            multi_nca_sum[[drug]] <- nca_sum
          }

          rv$multi_nca_results <- multi_nca
          rv$multi_nca_summary <- multi_nca_sum
          rv$nca_results <- NULL
          rv$nca_summary_results <- NULL
          rv$comp_results <- NULL
          rv$multi_comp_results <- NULL
          rv$multi_comp_summary <- NULL

          incProgress(0.2, detail = "Done!")
          rv$analysis_complete <- TRUE
          showNotification(
            paste("NCA completed for", length(drugs), "compounds"),
            type = "message", duration = 5
          )

        } else {
          # Compartmental per compound
          multi_comp <- list()
          multi_comp_sum <- list()
          any_success <- FALSE

          for (di in seq_along(drugs)) {
            drug <- drugs[di]
            incProgress(0.8 / length(drugs), detail = paste("Fitting", input$model_type, "for", drug, "..."))
            drug_data <- fdata[fdata$Drug == drug, ]

            n_subjects <- length(unique(drug_data$ID))
            if (n_subjects < 2) {
              showNotification(
                paste("Skipping", drug, "- need at least 2 subjects for compartmental modeling"),
                type = "warning", duration = 5
              )
              next
            }

            result <- fit_compartmental_model(
              drug_data, input$dose, input$model_type, route_for_model
            )

            if (!is.null(result$error)) {
              showNotification(
                paste(drug, ":", result$error),
                type = "warning", duration = 10
              )
            } else {
              multi_comp[[drug]] <- result
              any_success <- TRUE
            }
          }

          if (any_success) {
            rv$multi_comp_results <- multi_comp
            rv$nca_results <- NULL
            rv$nca_summary_results <- NULL
            rv$comp_results <- NULL
            rv$multi_nca_results <- NULL
            rv$multi_nca_summary <- NULL
            rv$analysis_complete <- TRUE

            incProgress(0.2, detail = "Done!")
            showNotification(
              paste(input$model_type, "fitted for", length(multi_comp), "of", length(drugs), "compounds"),
              type = "message", duration = 5
            )
          } else {
            rv$error_msg <- "Compartmental modeling failed for all compounds."
            showNotification(rv$error_msg, type = "error", duration = 10)
          }
        }

      } else {
        # ---- Single compound analysis ----
        # Clear multi-compound results
        rv$multi_nca_results <- NULL
        rv$multi_nca_summary <- NULL
        rv$multi_comp_results <- NULL
        rv$multi_comp_summary <- NULL

        n_subjects <- length(unique(fdata$ID))
        if (n_subjects < 2 && input$model_type != "NCA") {
          showNotification("Compartmental models require at least 2 subjects. Include more subjects or use NCA.",
                           type = "error")
          return()
        }

        if (input$model_type == "NCA") {
          incProgress(0.3, detail = "Computing NCA parameters...")

          rv$nca_results <- run_nca(fdata, input$dose, route_for_model)
          rv$nca_summary_results <- nca_summary(rv$nca_results)
          rv$comp_results <- NULL

          incProgress(0.7, detail = "Done!")
          rv$analysis_complete <- TRUE

          showNotification("NCA completed successfully", type = "message", duration = 5)

        } else {
          incProgress(0.2, detail = "Estimating initial parameters...")

          result <- fit_compartmental_model(
            fdata, input$dose, input$model_type, route_for_model
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
      is_multi <- !is.null(rv$multi_nca_results) || !is.null(rv$multi_comp_results)
      if (is_multi) {
        n_comp <- max(length(rv$multi_nca_results), length(rv$multi_comp_results))
        model_desc <- paste0(model_desc, " (", n_comp, " compounds)")
      }
      div(class = "status-box status-success", icon("check-circle"),
          paste("Analysis complete:", model_desc))
    }
  })

  # ---- Results UI ----
  output$results_ui <- renderUI({
    req(rv$analysis_complete)

    is_multi_nca <- !is.null(rv$multi_nca_results)
    is_multi_comp <- !is.null(rv$multi_comp_results)

    if (is_multi_nca) {
      # Multi-compound NCA results - stacked
      drugs <- names(rv$multi_nca_results)
      result_panels <- tagList()
      for (i in seq_along(drugs)) {
        drug <- drugs[i]
        result_panels <- tagList(result_panels,
          div(class = "compound-section",
            h3(drug),
            h4("Individual NCA Parameters", class = "section-title"),
            DTOutput(paste0("nca_individual_table_", i)),
            br(),
            h4("Summary Statistics", class = "section-title"),
            DTOutput(paste0("nca_summary_table_", i))
          ),
          br()
        )
      }
      tagList(
        result_panels,
        h4("Parameter Guide", class = "section-title"),
        uiOutput("nca_param_guide")
      )

    } else if (is_multi_comp) {
      # Multi-compound compartmental results - stacked
      drugs <- names(rv$multi_comp_results)
      result_panels <- tagList()
      for (i in seq_along(drugs)) {
        drug <- drugs[i]
        result_panels <- tagList(result_panels,
          div(class = "compound-section",
            h3(drug),
            h4("Population (Fixed Effect) Parameters", class = "section-title"),
            DTOutput(paste0("comp_pop_table_", i)),
            br(),
            h4("Individual Parameter Estimates", class = "section-title"),
            DTOutput(paste0("comp_indiv_table_", i)),
            br(),
            h4("Summary Statistics of Individual Parameters", class = "section-title"),
            DTOutput(paste0("comp_summary_table_", i)),
            br(),
            wellPanel(
              h5("Model Diagnostics"),
              uiOutput(paste0("model_info_", i))
            )
          ),
          br()
        )
      }
      tagList(
        result_panels,
        h4("Parameter Guide", class = "section-title"),
        uiOutput("comp_param_guide")
      )

    } else if (rv$analysis_type == "NCA") {
      # Single compound NCA
      is_iv <- input$route %in% c("IV", "IV_INF")
      tagList(
        h4("Individual NCA Parameters", class = "section-title"),
        DTOutput("nca_individual_table"),
        br(),
        h4("Summary Statistics", class = "section-title"),
        DTOutput("nca_summary_table"),
        br(),
        h4("Parameter Guide", class = "section-title"),
        uiOutput("nca_param_guide")
      )
    } else {
      # Single compound Compartmental
      tagList(
        h4("Population (Fixed Effect) Parameters", class = "section-title"),
        DTOutput("comp_pop_table"),
        br(),
        h4("Individual Parameter Estimates", class = "section-title"),
        DTOutput("comp_indiv_table"),
        br(),
        h4("Summary Statistics of Individual Parameters", class = "section-title"),
        DTOutput("comp_summary_table"),
        br(),
        wellPanel(
          h5("Model Diagnostics"),
          uiOutput("model_info")
        ),
        br(),
        h4("Parameter Guide", class = "section-title"),
        uiOutput("comp_param_guide")
      )
    }
  })

  # ---- NCA Parameter Guide ----
  output$nca_param_guide <- renderUI({
    is_iv <- input$route %in% c("IV", "IV_INF")
    cl_label <- if (is_iv) "CL" else "CL/F"
    vd_label <- if (is_iv) "Vd" else "Vd/F"
    vss_label <- if (is_iv) "Vss" else "Vss/F"

    guide <- tags$table(
      class = "table table-condensed table-striped",
      style = "font-size: 12px;",
      tags$thead(tags$tr(
        tags$th("Parameter"), tags$th("Full Name"), tags$th("Interpretation")
      )),
      tags$tbody(
        tags$tr(tags$td("Cmax"), tags$td("Maximum concentration"),
                tags$td("Highest observed drug concentration in plasma")),
        tags$tr(tags$td("Tmax"), tags$td("Time to maximum concentration"),
                tags$td("Time at which Cmax occurs; reflects absorption rate for extravascular routes")),
        tags$tr(tags$td("AUC_last"), tags$td("Area under the curve to last observation"),
                tags$td("Total drug exposure from dosing to last measurable concentration")),
        tags$tr(tags$td("AUC_inf"), tags$td("Area under the curve extrapolated to infinity"),
                tags$td("Total drug exposure including extrapolated terminal phase")),
        tags$tr(tags$td("AUC_extrap_pct"), tags$td("Percent AUC extrapolated"),
                tags$td("Fraction of AUC_inf estimated by extrapolation; ideally < 20%")),
        tags$tr(tags$td(HTML("&lambda;z")), tags$td("Terminal elimination rate constant"),
                tags$td("Rate of drug elimination during the terminal log-linear phase")),
        tags$tr(tags$td(HTML("t&frac12;")), tags$td("Terminal half-life"),
                tags$td("Time for plasma concentration to decrease by 50% during terminal phase")),
        tags$tr(tags$td(HTML("&lambda;z R&sup2;")), tags$td("Terminal phase R-squared"),
                tags$td("Goodness of fit for terminal slope estimation; should be > 0.90")),
        tags$tr(tags$td("MRT"), tags$td("Mean residence time"),
                tags$td("Average time a drug molecule stays in the body")),
        tags$tr(tags$td(cl_label), tags$td(if (is_iv) "Clearance" else "Apparent clearance"),
                tags$td("Volume of plasma completely cleared of drug per unit time")),
        tags$tr(tags$td(vd_label), tags$td(if (is_iv) "Volume of distribution" else "Apparent volume of distribution"),
                tags$td("Theoretical volume needed to contain the total drug at the same concentration as plasma")),
        tags$tr(tags$td(vss_label), tags$td(if (is_iv) "Volume of distribution at steady state" else "Apparent Vss"),
                tags$td("Vd when drug distribution is at equilibrium between compartments"))
      )
    )
    guide
  })

  # ---- Compartmental Parameter Guide ----
  output$comp_param_guide <- renderUI({
    is_iv <- input$route %in% c("IV", "IV_INF")
    pref <- if (is_iv) "" else "/F"
    model <- rv$analysis_type

    if (model == "1comp") {
      rows <- list(
        tags$tr(tags$td(paste0("V", pref)), tags$td(if (is_iv) "Volume of distribution" else "Apparent volume"),
                tags$td("Theoretical volume needed to contain the total drug at the same concentration as plasma")),
        tags$tr(tags$td(paste0("CL", pref)), tags$td(if (is_iv) "Clearance" else "Apparent clearance"),
                tags$td("Volume of plasma completely cleared of drug per unit time")),
        tags$tr(tags$td("k"), tags$td("Elimination rate constant"),
                tags$td("Fraction of drug eliminated per unit time")),
        tags$tr(tags$td(HTML("t&frac12;")), tags$td("Half-life"),
                tags$td("Time for plasma concentration to decrease by 50%"))
      )
    } else {
      rows <- list(
        tags$tr(tags$td(paste0("CL", pref)), tags$td(if (is_iv) "Clearance" else "Apparent clearance"),
                tags$td("Volume of plasma completely cleared of drug per unit time")),
        tags$tr(tags$td(paste0("V1", pref)), tags$td(if (is_iv) "Central volume" else "Apparent central volume"),
                tags$td("Volume of the central (plasma) compartment"))
      )
    }

    if (model %in% c("2comp", "3comp")) {
      rows <- c(rows, list(
        tags$tr(tags$td(paste0("Q", pref)), tags$td(if (is_iv) "Inter-compartmental clearance" else "Apparent Q"),
                tags$td("Rate of drug transfer between central and peripheral compartments")),
        tags$tr(tags$td(paste0("V2", pref)), tags$td(if (is_iv) "Peripheral volume" else "Apparent peripheral volume"),
                tags$td("Volume of the first peripheral (tissue) compartment")),
        tags$tr(tags$td(paste0("Vss", pref)), tags$td(if (is_iv) "Volume at steady state" else "Apparent Vss"),
                tags$td("Sum of all compartment volumes; total distribution volume at equilibrium"))
      ))
    }

    if (model == "3comp") {
      rows <- c(rows, list(
        tags$tr(tags$td(paste0("Q3", pref)), tags$td("Deep inter-compartmental clearance"),
                tags$td("Rate of drug transfer between central and deep peripheral compartment")),
        tags$tr(tags$td(paste0("V3", pref)), tags$td("Deep peripheral volume"),
                tags$td("Volume of the second (deep tissue) peripheral compartment"))
      ))
    }

    if (model == "2comp") {
      rows <- c(rows, list(
        tags$tr(tags$td("Alpha"), tags$td("Distribution rate constant"),
                tags$td("Rate of the rapid initial distribution phase")),
        tags$tr(tags$td("Beta"), tags$td("Terminal elimination rate constant"),
                tags$td("Rate of the slower terminal elimination phase")),
        tags$tr(tags$td(HTML("t&frac12; alpha")), tags$td("Distribution half-life"),
                tags$td("Half-life of the rapid distribution phase")),
        tags$tr(tags$td(HTML("t&frac12; beta")), tags$td("Terminal half-life"),
                tags$td("Half-life of the terminal elimination phase"))
      ))
    }

    if (!is_iv) {
      rows <- c(rows, list(
        tags$tr(tags$td("ka"), tags$td("Absorption rate constant"),
                tags$td("Rate at which drug is absorbed from the administration site into systemic circulation"))
      ))
    }

    tags$table(
      class = "table table-condensed table-striped",
      style = "font-size: 12px;",
      tags$thead(tags$tr(
        tags$th("Parameter"), tags$th("Full Name"), tags$th("Interpretation")
      )),
      tags$tbody(rows)
    )
  })

  # ---- NCA Tables (single compound) ----
  output$nca_individual_table <- renderDT({
    req(rv$nca_results)
    df <- rv$nca_results
    num_cols <- setdiff(names(df), "ID")
    datatable(df, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatSignif(columns = num_cols, digits = 4)
  })

  output$nca_summary_table <- renderDT({
    req(rv$nca_summary_results)
    datatable(rv$nca_summary_results, rownames = FALSE,
              colnames = c("Parameter", "N", "Mean", "SD", "SEM",
                           "Geo Mean", "GSD", "Median", "Min", "Max"),
              options = list(scrollX = TRUE, dom = 't')) %>%
      formatSignif(columns = c("Mean", "SD", "SEM", "Geo_Mean", "GSD",
                                "Median", "Min", "Max"), digits = 4)
  })

  # ---- Multi-compound NCA Tables (dynamic) ----
  observe({
    req(rv$multi_nca_results)
    drugs <- names(rv$multi_nca_results)

    for (i in seq_along(drugs)) {
      local({
        drug_idx <- i
        drug <- drugs[drug_idx]

        output[[paste0("nca_individual_table_", drug_idx)]] <- renderDT({
          req(rv$multi_nca_results)
          df <- rv$multi_nca_results[[drug]]
          req(df)
          num_cols <- setdiff(names(df), "ID")
          datatable(df, rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 20)) %>%
            formatSignif(columns = num_cols, digits = 4)
        })

        output[[paste0("nca_summary_table_", drug_idx)]] <- renderDT({
          req(rv$multi_nca_summary)
          df <- rv$multi_nca_summary[[drug]]
          req(df)
          datatable(df, rownames = FALSE,
                    colnames = c("Parameter", "N", "Mean", "SD", "SEM",
                                 "Geo Mean", "GSD", "Median", "Min", "Max"),
                    options = list(scrollX = TRUE, dom = 't')) %>%
            formatSignif(columns = c("Mean", "SD", "SEM", "Geo_Mean", "GSD",
                                      "Median", "Min", "Max"), digits = 4)
        })
      })
    }
  })

  # ---- Compartmental Tables (single compound) ----
  output$comp_pop_table <- renderDT({
    req(rv$comp_results)
    df <- rv$comp_results$summary
    if (all(c("SE", "CI_lower", "CI_upper") %in% names(df))) {
      display_df <- df[, c("Parameter", "Estimate", "SE", "CI_lower", "CI_upper")]
      datatable(display_df, rownames = FALSE,
                colnames = c("Parameter", "Estimate", "SE", "95% CI Lower", "95% CI Upper"),
                options = list(dom = 't', scrollX = TRUE)) %>%
        formatSignif(columns = c("Estimate", "SE", "CI_lower", "CI_upper"), digits = 4)
    } else {
      display_df <- df[, c("Parameter", "Estimate")]
      display_df$Estimate <- signif(display_df$Estimate, 4)
      datatable(display_df, rownames = FALSE,
                options = list(dom = 't', scrollX = TRUE))
    }
  })

  output$comp_indiv_table <- renderDT({
    req(rv$comp_results)
    df <- rv$comp_results$params
    num_cols <- setdiff(names(df), "ID")
    datatable(df, rownames = FALSE,
              options = list(scrollX = TRUE, pageLength = 20)) %>%
      formatSignif(columns = num_cols, digits = 4)
  })

  output$comp_summary_table <- renderDT({
    req(rv$comp_results)
    df <- rv$comp_results$params
    req(df, nrow(df) > 1)
    num_cols <- setdiff(names(df), "ID")

    summary_rows <- list()
    for (col in num_cols) {
      vals <- df[[col]]
      if (is.numeric(vals) && length(vals) > 1) {
        n <- sum(!is.na(vals))
        m <- mean(vals, na.rm = TRUE)
        s <- sd(vals, na.rm = TRUE)
        sem <- s / sqrt(n)
        ci_lower <- m - qt(0.975, df = n - 1) * sem
        ci_upper <- m + qt(0.975, df = n - 1) * sem
        geo_vals <- vals[vals > 0 & !is.na(vals)]
        geo_mean <- if (length(geo_vals) > 0) exp(mean(log(geo_vals))) else NA
        gsd <- if (length(geo_vals) > 1) exp(sd(log(geo_vals))) else NA
        med <- median(vals, na.rm = TRUE)
        mn <- min(vals, na.rm = TRUE)
        mx <- max(vals, na.rm = TRUE)
        summary_rows[[length(summary_rows) + 1]] <- data.frame(
          Parameter = col, N = n, Mean = m, SD = s, SEM = sem,
          CI_lower = ci_lower, CI_upper = ci_upper,
          Geo_Mean = geo_mean, GSD = gsd,
          Median = med, Min = mn, Max = mx,
          stringsAsFactors = FALSE
        )
      }
    }
    summary_df <- do.call(rbind, summary_rows)
    rv$comp_summary_stats <- summary_df

    datatable(summary_df, rownames = FALSE,
              colnames = c("Parameter", "N", "Mean", "SD", "SEM",
                           "95% CI Lower", "95% CI Upper",
                           "Geo Mean", "GSD", "Median", "Min", "Max"),
              options = list(scrollX = TRUE, dom = 't')) %>%
      formatSignif(columns = c("Mean", "SD", "SEM", "CI_lower", "CI_upper",
                                "Geo_Mean", "GSD", "Median", "Min", "Max"), digits = 4)
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

  # ---- Multi-compound Compartmental Tables (dynamic) ----
  observe({
    req(rv$multi_comp_results)
    drugs <- names(rv$multi_comp_results)

    for (i in seq_along(drugs)) {
      local({
        drug_idx <- i
        drug <- drugs[drug_idx]

        output[[paste0("comp_pop_table_", drug_idx)]] <- renderDT({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result)
          df <- result$summary
          if (all(c("SE", "CI_lower", "CI_upper") %in% names(df))) {
            display_df <- df[, c("Parameter", "Estimate", "SE", "CI_lower", "CI_upper")]
            datatable(display_df, rownames = FALSE,
                      colnames = c("Parameter", "Estimate", "SE", "95% CI Lower", "95% CI Upper"),
                      options = list(dom = 't', scrollX = TRUE)) %>%
              formatSignif(columns = c("Estimate", "SE", "CI_lower", "CI_upper"), digits = 4)
          } else {
            display_df <- df[, c("Parameter", "Estimate")]
            display_df$Estimate <- signif(display_df$Estimate, 4)
            datatable(display_df, rownames = FALSE,
                      options = list(dom = 't', scrollX = TRUE))
          }
        })

        output[[paste0("comp_indiv_table_", drug_idx)]] <- renderDT({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result)
          df <- result$params
          num_cols <- setdiff(names(df), "ID")
          datatable(df, rownames = FALSE,
                    options = list(scrollX = TRUE, pageLength = 20)) %>%
            formatSignif(columns = num_cols, digits = 4)
        })

        output[[paste0("comp_summary_table_", drug_idx)]] <- renderDT({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result)
          df <- result$params
          req(df, nrow(df) > 1)
          num_cols <- setdiff(names(df), "ID")

          summary_rows <- list()
          for (col in num_cols) {
            vals <- df[[col]]
            if (is.numeric(vals) && length(vals) > 1) {
              n <- sum(!is.na(vals))
              m <- mean(vals, na.rm = TRUE)
              s <- sd(vals, na.rm = TRUE)
              sem <- s / sqrt(n)
              ci_lower <- m - qt(0.975, df = n - 1) * sem
              ci_upper <- m + qt(0.975, df = n - 1) * sem
              geo_vals <- vals[vals > 0 & !is.na(vals)]
              geo_mean <- if (length(geo_vals) > 0) exp(mean(log(geo_vals))) else NA
              gsd <- if (length(geo_vals) > 1) exp(sd(log(geo_vals))) else NA
              med <- median(vals, na.rm = TRUE)
              mn <- min(vals, na.rm = TRUE)
              mx <- max(vals, na.rm = TRUE)
              summary_rows[[length(summary_rows) + 1]] <- data.frame(
                Parameter = col, N = n, Mean = m, SD = s, SEM = sem,
                CI_lower = ci_lower, CI_upper = ci_upper,
                Geo_Mean = geo_mean, GSD = gsd,
                Median = med, Min = mn, Max = mx,
                stringsAsFactors = FALSE
              )
            }
          }
          summary_df <- do.call(rbind, summary_rows)

          datatable(summary_df, rownames = FALSE,
                    colnames = c("Parameter", "N", "Mean", "SD", "SEM",
                                 "95% CI Lower", "95% CI Upper",
                                 "Geo Mean", "GSD", "Median", "Min", "Max"),
                    options = list(scrollX = TRUE, dom = 't')) %>%
            formatSignif(columns = c("Mean", "SD", "SEM", "CI_lower", "CI_upper",
                                      "Geo_Mean", "GSD", "Median", "Min", "Max"), digits = 4)
        })

        output[[paste0("model_info_", drug_idx)]] <- renderUI({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result)
          pop <- result$summary
          aic_val <- pop$AIC[1]
          bic_val <- pop$BIC[1]
          ll_val <- pop$logLik[1]

          tags$div(
            tags$p(tags$strong("AIC: "), round(aic_val, 2)),
            tags$p(tags$strong("BIC: "), round(bic_val, 2)),
            tags$p(tags$strong("Log-Likelihood: "), round(ll_val, 2))
          )
        })
      })
    }
  })

  # ---- Model Fit Plots ----
  output$model_fit_ui <- renderUI({
    req(rv$analysis_complete)
    req(rv$analysis_type != "NCA")

    is_multi_comp <- !is.null(rv$multi_comp_results)

    if (is_multi_comp) {
      drugs <- names(rv$multi_comp_results)
      fit_panels <- tagList()
      for (i in seq_along(drugs)) {
        drug <- drugs[i]
        fit_panels <- tagList(fit_panels,
          div(class = "compound-section",
            h3(drug),
            plotOutput(paste0("model_fit_plot_", i), height = "500px"),
            conditionalPanel(
              condition = "input.show_diagnostics",
              hr(),
              h4("Goodness-of-Fit Diagnostics", class = "section-title"),
              fluidRow(
                column(6, plotOutput(paste0("diag_obs_pred_", i), height = "350px")),
                column(6, plotOutput(paste0("diag_resid_", i), height = "350px"))
              ),
              fluidRow(
                column(6, plotOutput(paste0("diag_qq_", i), height = "350px")),
                column(6, plotOutput(paste0("diag_hist_", i), height = "350px"))
              )
            )
          ),
          br()
        )
      }
      fit_panels
    } else {
      tagList(
        plotOutput("model_fit_plot", height = "500px"),
        conditionalPanel(
          condition = "input.show_diagnostics",
          hr(),
          h4("Goodness-of-Fit Diagnostics", class = "section-title"),
          fluidRow(
            column(6,
              plotOutput("diag_obs_pred", height = "350px"),
              helpText(style = "font-size: 11px; color: #555; margin-top: 5px;",
                tags$strong("Observed vs. Predicted:"),
                "Points should scatter closely around the line of identity (dashed diagonal). ",
                "Systematic deviations above or below the line indicate model bias. ",
                "A good fit shows points evenly distributed along the line with no trends."
              )
            ),
            column(6,
              plotOutput("diag_resid", height = "350px"),
              helpText(style = "font-size: 11px; color: #555; margin-top: 5px;",
                tags$strong("Residuals vs. Predicted:"),
                "Residuals should be randomly scattered around zero (horizontal dashed line) with no pattern. ",
                "A funnel shape suggests heteroscedasticity (variance changes with concentration). ",
                "Systematic curvature indicates model misspecification."
              )
            )
          ),
          fluidRow(
            column(6,
              plotOutput("diag_qq", height = "350px"),
              helpText(style = "font-size: 11px; color: #555; margin-top: 5px;",
                tags$strong("Q-Q Plot (Normal Quantiles):"),
                "Points should follow the diagonal reference line if residuals are normally distributed. ",
                "S-shaped deviations indicate heavy or light tails. ",
                "Departures at the extremes are common with small samples but large deviations may suggest outliers or model issues."
              )
            ),
            column(6,
              plotOutput("diag_hist", height = "350px"),
              helpText(style = "font-size: 11px; color: #555; margin-top: 5px;",
                tags$strong("Residual Histogram:"),
                "The distribution of residuals should be approximately bell-shaped and centered near zero. ",
                "Strong skewness or multimodality may indicate model misspecification or outliers. ",
                "With few subjects, some asymmetry is expected."
              )
            )
          )
        )
      )
    }
  })

  # Single compound model fit plot
  output$model_fit_plot <- renderPlot({
    req(rv$comp_results, rv$comp_results$predictions)
    model_labels <- c(
      "1comp" = "One-Compartment",
      "2comp" = "Two-Compartment",
      "3comp" = "Three-Compartment"
    )
    route_label <- names(route_options)[route_options == input$route]

    plot_model_fit(
      filtered_data(),
      rv$comp_results$predictions,
      log_y = input$log_y_fit,
      model_label = paste(model_labels[rv$analysis_type], "-", route_label),
      time_unit = input$time_unit,
      conc_unit = input$conc_unit
    )
  })

  # ---- Diagnostic Plots (single compound) ----
  output$diag_obs_pred <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, filtered_data())
    plots$obs_vs_pred
  })

  output$diag_resid <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, filtered_data())
    plots$resid_vs_pred
  })

  output$diag_qq <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, filtered_data())
    plots$qq_plot
  })

  output$diag_hist <- renderPlot({
    req(rv$comp_results, rv$comp_results$fit)
    plots <- plot_diagnostics(rv$comp_results$fit, filtered_data())
    plots$resid_hist
  })

  # ---- Multi-compound Model Fit & Diagnostic Plots (dynamic) ----
  observe({
    req(rv$multi_comp_results)
    drugs <- names(rv$multi_comp_results)
    fdata <- filtered_data()

    for (i in seq_along(drugs)) {
      local({
        drug_idx <- i
        drug <- drugs[drug_idx]

        output[[paste0("model_fit_plot_", drug_idx)]] <- renderPlot({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result, result$predictions)
          drug_data <- filtered_data()[filtered_data()$Drug == drug, ]
          req(nrow(drug_data) > 0)

          model_labels <- c(
            "1comp" = "One-Compartment",
            "2comp" = "Two-Compartment",
            "3comp" = "Three-Compartment"
          )
          route_label <- names(route_options)[route_options == input$route]

          plot_model_fit(
            drug_data,
            result$predictions,
            log_y = input$log_y_fit,
            model_label = paste(drug, "-", model_labels[rv$analysis_type], "-", route_label),
            time_unit = input$time_unit,
            conc_unit = input$conc_unit
          )
        })

        output[[paste0("diag_obs_pred_", drug_idx)]] <- renderPlot({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result, result$fit)
          drug_data <- filtered_data()[filtered_data()$Drug == drug, ]
          plots <- plot_diagnostics(result$fit, drug_data)
          plots$obs_vs_pred
        })

        output[[paste0("diag_resid_", drug_idx)]] <- renderPlot({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result, result$fit)
          drug_data <- filtered_data()[filtered_data()$Drug == drug, ]
          plots <- plot_diagnostics(result$fit, drug_data)
          plots$resid_vs_pred
        })

        output[[paste0("diag_qq_", drug_idx)]] <- renderPlot({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result, result$fit)
          drug_data <- filtered_data()[filtered_data()$Drug == drug, ]
          plots <- plot_diagnostics(result$fit, drug_data)
          plots$qq_plot
        })

        output[[paste0("diag_hist_", drug_idx)]] <- renderPlot({
          req(rv$multi_comp_results)
          result <- rv$multi_comp_results[[drug]]
          req(result, result$fit)
          drug_data <- filtered_data()[filtered_data()$Drug == drug, ]
          plots <- plot_diagnostics(result$fit, drug_data)
          plots$resid_hist
        })
      })
    }
  })

  # ---- Manuscript Generation ----
  rv_manuscript <- reactiveValues(content = NULL, error = NULL, generating = FALSE, table = NULL)

  observeEvent(input$generate_manuscript, {
    req(rv$analysis_complete)

    api_key <- input$claude_api_key
    if (is.null(api_key) || nchar(trimws(api_key)) == 0) {
      showNotification("Please enter your Anthropic API key.", type = "error")
      return()
    }

    # Gather study info
    n_subj <- length(input$included_subjects)
    study_info <- list(
      drug_name   = if (nchar(input$drug_name) > 0) input$drug_name else "Unknown drug",
      species     = input$species,
      route       = input$route,
      dose        = input$dose,
      dose_unit   = input$dose_unit,
      time_unit   = input$time_unit,
      conc_unit   = input$conc_unit,
      n_subjects  = n_subj,
      body_weight = input$body_weight,
      weight_unit = input$weight_unit
    )

    # For multi-compound, use first compound's results for manuscript generation
    nca_res <- rv$nca_results
    nca_sum <- rv$nca_summary_results
    comp_res <- rv$comp_results
    comp_sum <- rv$comp_summary_stats

    if (!is.null(rv$multi_nca_results)) {
      nca_res <- rv$multi_nca_results[[1]]
      nca_sum <- rv$multi_nca_summary[[1]]
    }
    if (!is.null(rv$multi_comp_results)) {
      comp_res <- rv$multi_comp_results[[1]]
    }

    rv_manuscript$content <- NULL
    rv_manuscript$error <- NULL
    rv_manuscript$table <- NULL
    rv_manuscript$generating <- TRUE

    # Run generation
    withProgress(message = "Generating manuscript sections...", value = 0.3, {
      result <- generate_manuscript(
        api_key    = api_key,
        study_info = study_info,
        analysis_type = rv$analysis_type,
        nca_results   = nca_res,
        nca_summary   = nca_sum,
        comp_results  = comp_res,
        comp_summary_stats = comp_sum
      )
      incProgress(0.7, detail = "Done!")
    })

    rv_manuscript$generating <- FALSE

    if (!is.null(result$error)) {
      rv_manuscript$error <- result$error
      showNotification(result$error, type = "error", duration = 10)
    } else {
      rv_manuscript$content <- result$content
      rv_manuscript$table <- result$table
      showNotification("Manuscript sections generated!", type = "message", duration = 5)
    }
  })

  output$manuscript_output <- renderUI({
    if (isTRUE(rv_manuscript$generating)) {
      div(class = "status-box status-info", icon("spinner", class = "fa-spin"),
          "Generating manuscript sections... This may take 15-30 seconds.")
    } else if (!is.null(rv_manuscript$error)) {
      div(class = "status-box status-error", icon("exclamation-triangle"),
          rv_manuscript$error)
    } else if (!is.null(rv_manuscript$content)) {
      # Parse the content and insert table between RESULTS and INTERPRETATION
      content <- rv_manuscript$content
      tbl <- rv_manuscript$table

      # Split content at INTERPRETATION section to insert table before it
      content_parts <- strsplit(content, "(?i)(\\n\\s*INTERPRETATION[:\\s]*\\n)", perl = TRUE)[[1]]

      if (length(content_parts) >= 2 && !is.null(tbl)) {
        # Extract the interpretation header from original text
        interp_match <- regmatches(content, regexpr("(?i)\\n\\s*INTERPRETATION[:\\s]*\\n", content, perl = TRUE))
        before_interp <- content_parts[1]
        after_interp <- paste(content_parts[-1], collapse = interp_match)

        section_html <- tagList(
          # Methods + Results text
          div(class = "manuscript-section",
            HTML(gsub("\n", "<br/>", htmltools::htmlEscape(before_interp)))
          ),
          # Table 1
          div(class = "manuscript-section", style = "margin: 24px 0;",
            p(tags$b(tbl$title), style = "font-size: 13px; margin-bottom: 6px;"),
            HTML(tbl$html),
            p(tags$em(tbl$legend), style = "font-size: 11px; color: #555; margin-top: 6px;")
          ),
          # Interpretation
          div(class = "manuscript-section",
            HTML(gsub("\n", "<br/>", htmltools::htmlEscape(paste0(interp_match, after_interp))))
          )
        )
      } else {
        # Fallback: append table at the end if we can't split
        section_html <- tagList(
          div(class = "manuscript-section",
            HTML(gsub("\n", "<br/>", htmltools::htmlEscape(content)))
          ),
          if (!is.null(tbl)) {
            div(class = "manuscript-section", style = "margin: 24px 0;",
              p(tags$b(tbl$title), style = "font-size: 13px; margin-bottom: 6px;"),
              HTML(tbl$html),
              p(tags$em(tbl$legend), style = "font-size: 11px; color: #555; margin-top: 6px;")
            )
          }
        )
      }

      tagList(
        div(class = "manuscript-output",
          h4(icon("file-alt"), " Generated Manuscript Sections"),
          hr(),
          section_html
        ),
        br(),
        downloadButton("download_manuscript", "Download as Text File", class = "btn-block")
      )
    }
  })

  output$download_manuscript <- downloadHandler(
    filename = function() {
      paste0("pk_manuscript_", rv$analysis_type, "_", Sys.Date(), ".txt")
    },
    content = function(file) {
      lines <- rv_manuscript$content
      # Insert Table 1 between RESULTS and INTERPRETATION if available
      tbl <- rv_manuscript$table
      if (!is.null(tbl)) {
        parts <- strsplit(lines, "(?i)(\\n\\s*INTERPRETATION[:\\s]*\\n)", perl = TRUE)[[1]]
        if (length(parts) >= 2) {
          interp_match <- regmatches(lines, regexpr("(?i)\\n\\s*INTERPRETATION[:\\s]*\\n", lines, perl = TRUE))
          table_block <- paste0(
            "\n\n", tbl$title, "\n\n",
            tbl$text, "\n\n",
            tbl$legend, "\n"
          )
          lines <- paste0(parts[1], table_block, interp_match, paste(parts[-1], collapse = interp_match))
        } else {
          # Append at end
          lines <- paste0(lines, "\n\n", tbl$title, "\n\n", tbl$text, "\n\n", tbl$legend, "\n")
        }
      }
      writeLines(lines, file)
    }
  )

  # ---- Downloads ----
  output$download_params <- downloadHandler(
    filename = function() {
      paste0("pk_parameters_", rv$analysis_type, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      if (!is.null(rv$multi_nca_results)) {
        # Multi-compound NCA: combine all with a Drug column
        all_results <- list()
        for (drug in names(rv$multi_nca_results)) {
          df <- rv$multi_nca_results[[drug]]
          df$Drug <- drug
          all_results[[drug]] <- df
        }
        combined <- do.call(rbind, all_results)
        rownames(combined) <- NULL
        write.csv(combined, file, row.names = FALSE)
      } else if (!is.null(rv$multi_comp_results)) {
        # Multi-compound compartmental: combine population summaries
        all_results <- list()
        for (drug in names(rv$multi_comp_results)) {
          df <- rv$multi_comp_results[[drug]]$summary
          df$Drug <- drug
          all_results[[drug]] <- df
        }
        combined <- do.call(rbind, all_results)
        rownames(combined) <- NULL
        write.csv(combined, file, row.names = FALSE)
      } else if (rv$analysis_type == "NCA") {
        write.csv(rv$nca_results, file, row.names = FALSE)
      } else {
        pop <- rv$comp_results$summary
        write.csv(pop, file, row.names = FALSE)
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
