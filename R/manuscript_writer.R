# Manuscript Section Writer using Claude API
# Generates Materials & Methods, Results, and Interpretation sections

#' Build the system prompt for manuscript writing
manuscript_system_prompt <- function() {
  paste0(
    "You are an expert pharmacokineticist and scientific writer. ",
    "You write publication-quality manuscript sections for pharmacokinetic studies. ",
    "Use formal scientific language appropriate for peer-reviewed veterinary or clinical pharmacology journals. ",
    "Be precise with statistical terminology. ",
    "Report parameter estimates with appropriate significant figures. ",
    "Use past tense for methods and results. ",
    "Do not use markdown headers - just write flowing paragraph text for each section. ",
    "When referencing parameters, use standard PK abbreviations (e.g., CL, Vd, t1/2, AUC, Cmax, Tmax). ",
    "For extravascular routes, use apparent parameters (CL/F, Vd/F, etc.)."
  )
}

#' Build the user prompt from analysis context
#'
#' @param study_info List with drug_name, species, route, dose, dose_unit,
#'   time_unit, conc_unit, n_subjects
#' @param analysis_type Character: "NCA", "1comp", "2comp", "3comp"
#' @param results_text Character summary of results (parameter table as text)
#' @return Character prompt
build_manuscript_prompt <- function(study_info, analysis_type, results_text) {
  route_desc <- switch(study_info$route,
    "IV"     = "intravenous bolus",
    "IV_INF" = "intravenous infusion",
    "PO"     = "oral",
    "SC"     = "subcutaneous",
    "IM"     = "intramuscular",
    "IP"     = "intraperitoneal",
    "EV"     = "extravascular"
  )

  model_desc <- switch(analysis_type,
    "NCA"   = "non-compartmental analysis (NCA)",
    "1comp" = "one-compartment population pharmacokinetic model",
    "2comp" = "two-compartment population pharmacokinetic model",
    "3comp" = "three-compartment population pharmacokinetic model"
  )

  is_compartmental <- analysis_type != "NCA"
  is_iv <- study_info$route %in% c("IV", "IV_INF")

  method_context <- if (is_compartmental) {
    paste0(
      "The data were analyzed using a ", model_desc,
      " fitted via nonlinear mixed-effects modeling (NLME) with the nlme package in R. ",
      "Parameters were estimated on the log scale to ensure positivity. ",
      if (!is_iv) "Since the route was extravascular, apparent parameters (divided by bioavailability F) are reported. " else "",
      "Between-subject variability was modeled using random effects. ",
      "Summary statistics for individual parameter estimates included mean, standard deviation (SD), ",
      "standard error of the mean (SEM), 95% confidence interval (CI) of the mean, geometric mean, ",
      "geometric standard deviation (GSD), median, and range."
    )
  } else {
    paste0(
      "The data were analyzed using ", model_desc,
      " with the linear trapezoidal rule for AUC calculation. ",
      "The terminal elimination rate constant (lambda_z) was estimated by log-linear regression ",
      "of the terminal phase concentrations. ",
      if (!is_iv) "Since the route was extravascular, apparent parameters (divided by bioavailability F) are reported. " else "",
      "Summary statistics included mean, standard deviation (SD), standard error of the mean (SEM), ",
      "95% confidence interval (CI) of the mean, geometric mean, geometric standard deviation (GSD), median, and range."
    )
  }

  software_context <- paste0(
    "All pharmacokinetic analyses were performed using R (R Foundation for Statistical Computing, Vienna, Austria; ",
    "https://www.R-project.org/) with the following packages: shiny (web application framework), ",
    "nlme (nonlinear mixed-effects modeling), ggplot2 (data visualization), DT (interactive tables), ",
    "readxl (data import), and scales (axis formatting). ",
    "A custom pharmacokinetic analysis application was developed with the assistance of Claude (Anthropic) ",
    "and is freely available at https://github.com/hbeaufrere/PK-Analysis for transparency and reproducibility."
  )

  # Body weight info
  weight_info <- ""
  if (!is.null(study_info$body_weight) && !is.na(study_info$body_weight) && study_info$body_weight > 0) {
    weight_info <- paste0("- Mean body weight: ", study_info$body_weight, " ", study_info$weight_unit, "\n")
  }

  prompt <- paste0(
    "Write three sections for a pharmacokinetic manuscript based on the following study:\n\n",
    "STUDY DETAILS:\n",
    "- Drug: ", study_info$drug_name, "\n",
    "- Species: ", study_info$species, "\n",
    "- Route: ", route_desc, "\n",
    "- Dose: ", study_info$dose, " ", study_info$dose_unit, "\n",
    weight_info,
    "- Number of subjects: ", study_info$n_subjects, "\n",
    "- Time unit: ", study_info$time_unit, "\n",
    "- Concentration unit: ", study_info$conc_unit, "\n",
    "- Analysis method: ", model_desc, "\n\n",
    "METHODOLOGICAL CONTEXT:\n", method_context, "\n\n",
    "SOFTWARE AND TOOLS:\n", software_context, "\n\n",
    "RESULTS:\n", results_text, "\n\n",
    "Please write the following three sections. Label each section clearly:\n\n",
    "STATISTICAL ANALYSIS (Materials and Methods):\n",
    "Write 1-2 paragraphs describing the statistical/pharmacokinetic analysis methods used. ",
    "Include the specific R software version and R packages used (nlme, ggplot2, shiny, etc.), ",
    "the modeling approach, how parameters were estimated, and how summary statistics were computed ",
    "(mean, SD, SEM, 95% CI of the mean, geometric mean, GSD, median, range). ",
    "Mention that the analysis was performed using a custom pharmacokinetic application ",
    "developed with the assistance of Claude (Anthropic, San Francisco, CA), available at ",
    "https://github.com/hbeaufrere/PK-Analysis.\n\n",
    "RESULTS:\n",
    "Write 2-3 paragraphs reporting the key pharmacokinetic parameters with their values. ",
    "Report mean \u00B1 SD (or SEM) and 95% CI of the mean, or geometric mean (GSD) as appropriate. ",
    "Describe the concentration-time profile and the key PK findings. ",
    "Reference 'Table 1' when directing the reader to the full parameter summary ",
    "(e.g., 'Pharmacokinetic parameters are summarized in Table 1.'). ",
    "A formatted Table 1 with all parameters will be automatically inserted after this section. ",
    "Include a paragraph discussing whether the chosen statistical model (", model_desc,
    ") was appropriate for these data. Evaluate goodness-of-fit metrics (AIC, BIC) if available, ",
    "comment on residual diagnostics, and discuss whether the model assumptions were reasonable ",
    "for this type of pharmacokinetic data.\n\n",
    "INTERPRETATION:\n",
    "Write 1-2 paragraphs interpreting the pharmacokinetic results. ",
    "Discuss what the parameters suggest about the drug's disposition in this species. ",
    "Compare to general pharmacokinetic expectations where appropriate. ",
    "Note any clinically relevant implications."
  )

  return(prompt)
}

#' Format analysis results into text for the prompt
#'
#' @param analysis_type Character
#' @param nca_results Data frame of individual NCA results (or NULL)
#' @param nca_summary Data frame of NCA summary stats (or NULL)
#' @param comp_results List with $summary and $params (or NULL)
#' @return Character string summarizing the results
format_results_for_prompt <- function(analysis_type, nca_results = NULL,
                                      nca_summary = NULL, comp_results = NULL) {
  if (analysis_type == "NCA") {
    if (is.null(nca_summary)) return("No summary results available.")
    lines <- c("NCA Summary Statistics (Mean \u00B1 SD [Min - Max]):")
    for (i in seq_len(nrow(nca_summary))) {
      row <- nca_summary[i, ]
      lines <- c(lines, paste0(
        "  ", row$Parameter, ": Mean = ", signif(row$Mean, 4),
        ", SD = ", signif(row$SD, 4),
        ", Median = ", signif(row$Median, 4),
        ", Range = [", signif(row$Min, 4), " - ", signif(row$Max, 4), "]",
        ", N = ", row$N
      ))
    }
    if (!is.null(nca_results)) {
      lines <- c(lines, paste0("\nNumber of subjects: ", nrow(nca_results)))
    }
    return(paste(lines, collapse = "\n"))
  } else {
    # Compartmental results
    if (is.null(comp_results)) return("No compartmental results available.")
    pop <- comp_results$summary
    indiv <- comp_results$params

    has_se <- "SE" %in% names(pop) && "CI_lower" %in% names(pop)
    lines <- c("Population (Fixed Effect) Parameter Estimates (Estimate [SE; 95% CI]):")
    for (i in seq_len(nrow(pop))) {
      if (pop$Parameter[i] %in% c("AIC", "BIC", "logLik")) next
      line <- paste0("  ", pop$Parameter[i], " = ", signif(pop$Estimate[i], 4))
      if (has_se && !is.na(pop$SE[i])) {
        line <- paste0(line, ", SE = ", signif(pop$SE[i], 4),
                       ", 95% CI = [", signif(pop$CI_lower[i], 4),
                       " - ", signif(pop$CI_upper[i], 4), "]")
      }
      lines <- c(lines, line)
    }

    # Add model fit statistics
    aic_row <- which(names(pop) == "AIC")
    if (length(aic_row) > 0 && "AIC" %in% names(pop)) {
      lines <- c(lines, paste0("\nModel Fit: AIC = ", round(pop$AIC[1], 2),
                                ", BIC = ", round(pop$BIC[1], 2),
                                ", Log-Likelihood = ", round(pop$logLik[1], 2)))
    }

    # Individual parameter summary with SEM and 95% CI
    if (!is.null(indiv) && nrow(indiv) > 0) {
      num_cols <- setdiff(names(indiv), "ID")
      lines <- c(lines, "\nIndividual Parameter Summary (Mean \u00B1 SD [SEM; 95% CI]):")
      for (col in num_cols) {
        vals <- indiv[[col]]
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
          lines <- c(lines, paste0(
            "  ", col, ": Mean = ", signif(m, 4),
            " \u00B1 ", signif(s, 4),
            ", SEM = ", signif(sem, 4),
            ", 95% CI = [", signif(ci_lower, 4), " - ", signif(ci_upper, 4), "]",
            if (!is.na(geo_mean)) paste0(", Geometric Mean = ", signif(geo_mean, 4)) else "",
            if (!is.na(gsd)) paste0(", GSD = ", signif(gsd, 4)) else ""
          ))
        }
      }
      lines <- c(lines, paste0("\nNumber of subjects: ", length(unique(indiv$ID))))
    }

    return(paste(lines, collapse = "\n"))
  }
}

#' Call Claude API to generate manuscript sections
#'
#' @param api_key Character API key for Claude
#' @param system_prompt Character system prompt
#' @param user_prompt Character user prompt
#' @param model Character model ID (default: claude-sonnet-4-5-20250929)
#' @return List with 'content' (character) or 'error' (character)
call_claude_api <- function(api_key, system_prompt, user_prompt,
                            model = "claude-sonnet-4-5-20250929") {
  if (!requireNamespace("httr", quietly = TRUE)) {
    return(list(content = NULL, error = "The 'httr' package is required. Install it with: install.packages('httr')"))
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    return(list(content = NULL, error = "The 'jsonlite' package is required. Install it with: install.packages('jsonlite')"))
  }

  body <- list(
    model = model,
    max_tokens = 4096,
    system = system_prompt,
    messages = list(
      list(role = "user", content = user_prompt)
    )
  )

  response <- tryCatch({
    httr::POST(
      url = "https://api.anthropic.com/v1/messages",
      httr::add_headers(
        `x-api-key` = api_key,
        `anthropic-version` = "2023-06-01",
        `content-type` = "application/json"
      ),
      body = jsonlite::toJSON(body, auto_unbox = TRUE),
      encode = "raw",
      httr::timeout(120)
    )
  }, error = function(e) {
    return(list(content = NULL, error = paste("API request failed:", e$message)))
  })

  if (inherits(response, "list") && !is.null(response$error)) {
    return(response)
  }

  status <- httr::status_code(response)
  resp_body <- httr::content(response, as = "text", encoding = "UTF-8")
  parsed <- tryCatch(jsonlite::fromJSON(resp_body), error = function(e) NULL)

  if (status != 200) {
    err_msg <- if (!is.null(parsed$error$message)) parsed$error$message else resp_body
    return(list(content = NULL, error = paste0("API error (", status, "): ", err_msg)))
  }

  if (is.null(parsed) || is.null(parsed$content)) {
    return(list(content = NULL, error = "Could not parse API response."))
  }

  # Extract text from content blocks
  text <- paste(parsed$content$text, collapse = "\n")
  return(list(content = text, error = NULL))
}

#' Build a publication-ready PK parameter table for the manuscript
#'
#' Generates both an HTML table (for display) and a plain-text table (for download).
#' For NCA: uses nca_summary. For compartmental: uses comp_summary_stats.
#'
#' @param study_info List with drug_name, species, route, dose, dose_unit, conc_unit, time_unit, n_subjects
#' @param analysis_type Character: "NCA", "1comp", "2comp", "3comp"
#' @param nca_summary Data frame from nca_summary() or NULL
#' @param comp_summary_stats Data frame of individual param summary stats or NULL
#' @return List with 'title' (character), 'legend' (character),
#'   'html' (character HTML string), 'text' (character plain-text table)
build_manuscript_table <- function(study_info, analysis_type,
                                    nca_summary = NULL, comp_summary_stats = NULL) {
  is_iv <- study_info$route %in% c("IV", "IV_INF")
  route_desc <- switch(study_info$route,
    "IV"     = "intravenous",
    "IV_INF" = "intravenous infusion",
    "PO"     = "oral",
    "SC"     = "subcutaneous",
    "IM"     = "intramuscular",
    "IP"     = "intraperitoneal",
    "EV"     = "extravascular",
    tolower(study_info$route)
  )
  model_desc <- switch(analysis_type,
    "NCA"   = "non-compartmental analysis",
    "1comp" = "one-compartment model",
    "2comp" = "two-compartment model",
    "3comp" = "three-compartment model"
  )

  # Table title
  title <- paste0(
    "Table 1. Pharmacokinetic parameters of ", study_info$drug_name,
    " following ", route_desc, " administration at ",
    study_info$dose, " ", study_info$dose_unit, " to ",
    study_info$n_subjects, " ", study_info$species,
    " determined by ", model_desc, "."
  )

  # Parameter label mapping (with units) for NCA
  conc_unit <- study_info$conc_unit
  time_unit <- study_info$time_unit
  cl_label <- if (is_iv) "CL" else "CL/F"
  vd_label <- if (is_iv) "Vd" else "Vd/F"
  vss_label <- if (is_iv) "Vss" else "Vss/F"

  nca_labels <- c(
    Cmax        = paste0("C\u2098\u2090\u2093 (", conc_unit, ")"),
    Tmax        = paste0("T\u2098\u2090\u2093 (", time_unit, ")"),
    AUC_last    = paste0("AUC\u2097\u2090\u209B\u209C (", conc_unit, "\u00B7", time_unit, ")"),
    AUC_inf     = paste0("AUC\u221E (", conc_unit, "\u00B7", time_unit, ")"),
    AUC_extrap_pct = "AUC extrapolated (%)",
    Lambda_z    = paste0("\u03BB\u1D63 (1/", time_unit, ")"),
    Half_life   = paste0("t\u00BD (", time_unit, ")"),
    Lambda_z_R2 = "Terminal R\u00B2",
    MRT         = paste0("MRT (", time_unit, ")"),
    CL          = paste0(cl_label, " (L/", time_unit, ")"),
    Vd          = paste0(vd_label, " (L)"),
    Vss         = paste0(vss_label, " (L)")
  )

  # Select data source
  if (analysis_type == "NCA" && !is.null(nca_summary)) {
    summ <- nca_summary
    label_map <- nca_labels
  } else if (!is.null(comp_summary_stats)) {
    summ <- comp_summary_stats
    label_map <- NULL  # compartmental params keep their own names
  } else {
    return(NULL)
  }

  # Build rows
  fmt <- function(x, digits = 4) {
    if (is.na(x)) return("\u2014")
    formatC(signif(x, digits), format = "fg", flag = "")
  }

  rows_html <- character(0)
  rows_text <- character(0)

  has_ci <- all(c("CI_lower", "CI_upper") %in% names(summ))

  for (i in seq_len(nrow(summ))) {
    row <- summ[i, ]
    param_raw <- as.character(row$Parameter)

    # Get display label
    if (!is.null(label_map) && param_raw %in% names(label_map)) {
      param_label <- label_map[[param_raw]]
    } else {
      param_label <- param_raw
    }

    n_val     <- row$N
    mean_val  <- fmt(row$Mean)
    sd_val    <- fmt(row$SD)
    mean_sd   <- paste0(mean_val, " \u00B1 ", sd_val)
    median_val <- fmt(row$Median)
    range_val  <- paste0(fmt(row$Min), "\u2013", fmt(row$Max))

    # Geometric mean (GSD) if available
    geo_col <- if ("Geo_Mean" %in% names(row)) row$Geo_Mean else NA
    gsd_col <- if ("GSD" %in% names(row)) row$GSD else NA
    if (!is.na(geo_col)) {
      if (!is.na(gsd_col)) {
        geo_str <- paste0(fmt(geo_col), " (", fmt(gsd_col), ")")
      } else {
        geo_str <- fmt(geo_col)
      }
    } else {
      geo_str <- "\u2014"
    }

    # HTML row
    rows_html <- c(rows_html, paste0(
      "<tr>",
      "<td style='text-align:left; padding:4px 8px;'>", htmltools::htmlEscape(param_label), "</td>",
      "<td style='text-align:center; padding:4px 8px;'>", htmltools::htmlEscape(mean_sd), "</td>",
      "<td style='text-align:center; padding:4px 8px;'>", htmltools::htmlEscape(median_val), "</td>",
      "<td style='text-align:center; padding:4px 8px;'>", htmltools::htmlEscape(range_val), "</td>",
      "<td style='text-align:center; padding:4px 8px;'>", htmltools::htmlEscape(geo_str), "</td>",
      "</tr>"
    ))

    # Text row (tab-delimited)
    rows_text <- c(rows_text, paste(param_label, mean_sd, median_val, range_val, geo_str, sep = "\t"))
  }

  # Column headers
  col_headers <- c("Parameter", "Mean \u00B1 SD", "Median", "Range", "Geometric Mean (GSD)")

  # HTML table
  header_html <- paste0(
    "<tr>",
    paste0("<th style='text-align:center; padding:6px 8px; border-bottom:2px solid #333;'>",
           htmltools::htmlEscape(col_headers), "</th>", collapse = ""),
    "</tr>"
  )
  # Left-align the Parameter header
  header_html <- sub("text-align:center;(.*?)>Parameter", "text-align:left;\\1>Parameter", header_html)

  html_table <- paste0(
    "<table style='border-collapse:collapse; width:100%; font-size:13px; margin:10px 0;",
    " border-top:2px solid #333; border-bottom:2px solid #333;'>",
    "<thead>", header_html, "</thead>",
    "<tbody>", paste(rows_html, collapse = ""), "</tbody>",
    "</table>"
  )

  # Legend
  legend_parts <- c(
    paste0("Data are presented as mean \u00B1 standard deviation (SD), median, range, ",
           "and geometric mean with geometric standard deviation (GSD) in parentheses.")
  )
  if (!is_iv) {
    legend_parts <- c(legend_parts,
      "/F indicates parameters are apparent values not corrected for bioavailability.")
  }
  legend_parts <- c(legend_parts,
    paste0("n = ", summ$N[1], " ", study_info$species, "."),
    paste0("\u03BB\u1D63, terminal elimination rate constant; ",
           "AUC, area under the concentration-time curve; ",
           "C\u2098\u2090\u2093, maximum concentration; ",
           "CL, clearance; MRT, mean residence time; ",
           "t\u00BD, elimination half-life; ",
           "T\u2098\u2090\u2093, time to maximum concentration; ",
           "Vd, volume of distribution; ",
           "Vss, volume of distribution at steady state."))
  legend <- paste(legend_parts, collapse = " ")

  # Plain-text table
  text_header <- paste(col_headers, collapse = "\t")
  text_table <- paste(c(text_header, rows_text), collapse = "\n")

  return(list(
    title  = title,
    legend = legend,
    html   = html_table,
    text   = text_table
  ))
}

#' Main function: generate manuscript sections
#'
#' @param api_key Character Claude API key
#' @param study_info List of study information
#' @param analysis_type Character analysis type
#' @param nca_results Data frame (or NULL)
#' @param nca_summary Data frame (or NULL)
#' @param comp_results List (or NULL)
#' @return List with 'content' or 'error'
generate_manuscript <- function(api_key, study_info, analysis_type,
                                nca_results = NULL, nca_summary = NULL,
                                comp_results = NULL, comp_summary_stats = NULL) {
  # Build results text
  results_text <- format_results_for_prompt(
    analysis_type, nca_results, nca_summary, comp_results
  )

  # Build prompts
  sys_prompt <- manuscript_system_prompt()
  user_prompt <- build_manuscript_prompt(study_info, analysis_type, results_text)

  # Call API
  result <- call_claude_api(api_key, sys_prompt, user_prompt)

  # Build manuscript table from actual data
  result$table <- build_manuscript_table(
    study_info, analysis_type,
    nca_summary = nca_summary,
    comp_summary_stats = comp_summary_stats
  )

  return(result)
}
