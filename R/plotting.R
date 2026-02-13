# Plotting Module for PK Analysis

library(ggplot2)

#' Create mean +/- SEM concentration-time plot
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param log_y Logical, use log-transformed y-axis
#' @param species Character, species name for plot title
#' @param drug_name Character, drug name for plot title
#' @param time_unit Character, time unit for axis label
#' @param conc_unit Character, concentration unit for axis label
#' @param color_mode Character, "color" or "bw"
#' @param show_ci Logical, show SEM ribbon/error bars
#' @return ggplot object
plot_mean_conc_time <- function(df, log_y = FALSE, species = "", drug_name = "",
                                 time_unit = "h", conc_unit = "ng/mL",
                                 color_mode = "color", show_ci = TRUE,
                                 point_size = 2.5, line_width = 1,
                                 errorbar_width = 0.5, show_grid = TRUE) {
  # Compute mean and SEM at each time point
  summary_df <- aggregate(Conc ~ Time, data = df, FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE),
      sem = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))),
      n = sum(!is.na(x)))
  })
  summary_df <- data.frame(
    Time = summary_df$Time,
    Mean = summary_df$Conc[, "mean"],
    SEM = summary_df$Conc[, "sem"],
    N = summary_df$Conc[, "n"]
  )
  summary_df$Lower <- summary_df$Mean - summary_df$SEM
  summary_df$Upper <- summary_df$Mean + summary_df$SEM

  # For log scale, floor lower bound at small positive value
  if (log_y) {
    min_pos <- min(summary_df$Mean[summary_df$Mean > 0], na.rm = TRUE) * 0.01
    summary_df$Lower <- pmax(summary_df$Lower, min_pos)
    summary_df$Mean[summary_df$Mean <= 0] <- NA
  }

  title_parts <- c()
  if (nzchar(drug_name)) title_parts <- c(title_parts, drug_name)
  if (nzchar(species)) title_parts <- c(title_parts, paste0("(", species, ")"))
  plot_title <- paste("Concentration-Time Profile",
                      if (length(title_parts) > 0) paste("-", paste(title_parts, collapse = " ")) else "")

  is_bw <- (color_mode == "bw")
  line_col <- if (is_bw) "black" else "#2166AC"
  fill_col <- if (is_bw) "grey60" else "#2166AC"

  p <- ggplot(summary_df, aes(x = Time, y = Mean))

  if (show_ci) {
    p <- p +
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = fill_col)
  }

  p <- p +
    geom_errorbar(aes(ymin = Lower, ymax = Upper),
                  width = max(diff(range(summary_df$Time))) * 0.015,
                  color = line_col, linewidth = errorbar_width) +
    geom_line(color = line_col, linewidth = line_width) +
    geom_point(color = line_col, size = point_size) +
    labs(
      title = plot_title,
      x = paste0("Time (", time_unit, ")"),
      y = paste0("Concentration (", conc_unit, ")"),
      caption = paste("Mean \u00B1 SEM, n =", max(summary_df$N))
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_line(color = "grey90"),
      plot.caption = element_text(hjust = 0, size = 10, color = "grey50")
    )

  if (!show_grid) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  if (log_y) {
    p <- p + scale_y_log10(
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
      labs(y = paste0("Concentration (", conc_unit, ") - Log Scale"))
  }

  return(p)
}

#' Create individual concentration-time plots (spaghetti plot)
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param log_y Logical, use log-transformed y-axis
#' @param time_unit Character, time unit
#' @param conc_unit Character, concentration unit
#' @param color_mode Character, "color" or "bw"
#' @return ggplot object
plot_individual_conc_time <- function(df, log_y = FALSE,
                                       time_unit = "h", conc_unit = "ng/mL",
                                       color_mode = "color",
                                       point_size = 1.5, line_width = 0.7,
                                       show_grid = TRUE) {
  is_bw <- (color_mode == "bw")
  n_subj <- length(unique(df$ID))

  if (is_bw) {
    # Use different linetypes and shapes for BW mode
    p <- ggplot(df, aes(x = Time, y = Conc, group = ID, linetype = ID, shape = ID)) +
      geom_line(linewidth = line_width, alpha = 0.8) +
      geom_point(size = point_size, alpha = 0.9) +
      labs(linetype = "Subject", shape = "Subject")
    if (n_subj <= 6) {
      p <- p + scale_linetype_manual(values = rep(c("solid", "dashed", "dotted",
                                                      "dotdash", "longdash", "twodash"), length.out = n_subj))
    }
  } else {
    p <- ggplot(df, aes(x = Time, y = Conc, group = ID, color = ID)) +
      geom_line(linewidth = line_width, alpha = 0.7) +
      geom_point(size = point_size, alpha = 0.8) +
      labs(color = "Subject")
  }

  p <- p +
    labs(
      title = "Individual Concentration-Time Profiles",
      x = paste0("Time (", time_unit, ")"),
      y = paste0("Concentration (", conc_unit, ")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  if (!show_grid) {
    p <- p + theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  }

  if (log_y) {
    p <- p + scale_y_log10(
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
      labs(y = paste0("Concentration (", conc_unit, ") - Log Scale"))
  }

  # Limit legend entries for large datasets
  if (n_subj > 20) {
    p <- p + theme(legend.position = "none") +
      labs(caption = paste(n_subj, "subjects (legend hidden)"))
  }

  return(p)
}

#' Plot compartmental model fit with observed data
#'
#' @param df Observed data with ID, Time, Conc
#' @param predictions List with 'population' and 'individual' prediction data frames
#' @param log_y Logical, log-transformed y-axis
#' @param model_label Character, model description for title
#' @param time_unit Character, time unit
#' @param conc_unit Character, concentration unit
#' @return ggplot object
plot_model_fit <- function(df, predictions, log_y = FALSE, model_label = "",
                            time_unit = "h", conc_unit = "ng/mL") {
  pop_pred <- predictions$population
  indiv_pred <- predictions$individual

  # Observed summary
  obs_summary <- aggregate(Conc ~ Time, data = df, FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE), sem = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
  })
  obs_summary <- data.frame(
    Time = obs_summary$Time,
    Mean = obs_summary$Conc[, "mean"],
    SEM = obs_summary$Conc[, "sem"]
  )

  p <- ggplot() +
    # Population prediction line
    geom_line(data = pop_pred, aes(x = Time, y = Pred_pop),
              color = "#D6604D", linewidth = 1.2, linetype = "solid") +
    # Individual predictions (thin lines)
    geom_line(data = indiv_pred, aes(x = Time, y = Pred_indiv, group = ID),
              color = "grey60", linewidth = 0.4, alpha = 0.6) +
    # Observed data points (mean +/- SEM)
    geom_point(data = obs_summary, aes(x = Time, y = Mean),
               color = "#2166AC", size = 3) +
    geom_errorbar(data = obs_summary, aes(x = Time, ymin = Mean - SEM, ymax = Mean + SEM),
                  color = "#2166AC", width = max(diff(range(obs_summary$Time))) * 0.015,
                  linewidth = 0.5) +
    labs(
      title = paste("Model Fit:", model_label),
      x = paste0("Time (", time_unit, ")"),
      y = paste0("Concentration (", conc_unit, ")")
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  if (log_y) {
    min_pos <- min(c(obs_summary$Mean[obs_summary$Mean > 0],
                     pop_pred$Pred_pop[pop_pred$Pred_pop > 0]), na.rm = TRUE) * 0.01
    p <- p + scale_y_log10(
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
      labs(y = paste0("Concentration (", conc_unit, ") - Log Scale"))
  }

  return(p)
}

#' Diagnostic plots for compartmental model
#'
#' @param fit nlme fit object
#' @param df Original data
#' @return List of ggplot objects
plot_diagnostics <- function(fit, df) {
  # Get fitted values and residuals (handle both nlme and gnls)
  fitted_vals <- fitted(fit)
  resid_vals <- tryCatch(
    residuals(fit, type = "normalized"),
    error = function(e) residuals(fit) / sd(residuals(fit))
  )
  observed <- df$Conc[!is.na(df$Conc)]

  # Trim to same length
  n <- min(length(observed), length(fitted_vals))
  diag_df <- data.frame(
    Observed = observed[1:n],
    Fitted = fitted_vals[1:n],
    Residuals = resid_vals[1:n]
  )

  # Observed vs Predicted
  p1 <- ggplot(diag_df, aes(x = Fitted, y = Observed)) +
    geom_point(alpha = 0.6, color = "#2166AC") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Observed vs. Population Predicted",
         x = "Population Predicted", y = "Observed") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    coord_equal()

  # Residuals vs Predicted
  p2 <- ggplot(diag_df, aes(x = Fitted, y = Residuals)) +
    geom_point(alpha = 0.6, color = "#2166AC") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Normalized Residuals vs. Predicted",
         x = "Population Predicted", y = "Normalized Residuals") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # QQ plot of residuals
  p3 <- ggplot(diag_df, aes(sample = Residuals)) +
    stat_qq(color = "#2166AC", alpha = 0.6) +
    stat_qq_line(color = "red", linetype = "dashed") +
    labs(title = "Q-Q Plot of Normalized Residuals",
         x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  # Histogram of residuals
  p4 <- ggplot(diag_df, aes(x = Residuals)) +
    geom_histogram(bins = 20, fill = "#2166AC", color = "white", alpha = 0.7) +
    labs(title = "Distribution of Normalized Residuals",
         x = "Normalized Residuals", y = "Count") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  return(list(obs_vs_pred = p1, resid_vs_pred = p2, qq_plot = p3, resid_hist = p4))
}
