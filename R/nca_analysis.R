# Non-Compartmental Analysis (NCA) Module

#' Compute AUC using the linear-log trapezoidal rule
#'
#' Uses linear trapezoidal for ascending concentrations and
#' log trapezoidal for descending concentrations.
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @return Numeric AUC value
auc_lin_log <- function(time, conc) {
  n <- length(time)
  if (n < 2) return(NA_real_)

  auc <- 0
  for (i in 2:n) {
    dt <- time[i] - time[i - 1]
    c1 <- conc[i - 1]
    c2 <- conc[i]

    if (dt <= 0) next
    if (c1 <= 0 || c2 <= 0) {
      # Linear trapezoidal if any conc is zero or negative
      auc <- auc + dt * (c1 + c2) / 2
    } else if (c2 >= c1) {
      # Ascending: linear trapezoidal
      auc <- auc + dt * (c1 + c2) / 2
    } else {
      # Descending: log trapezoidal
      auc <- auc + dt * (c1 - c2) / log(c1 / c2)
    }
  }
  return(auc)
}

#' Estimate terminal elimination rate constant (lambda_z)
#'
#' Uses log-linear regression on the terminal phase (last 3+ points
#' in the descending portion). Selects the best fit (highest R²)
#' using at least 3 terminal points.
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @return List with lambda_z, t_half, r_squared, n_points
estimate_lambda_z <- function(time, conc) {
  n <- length(time)
  if (n < 3) {
    return(list(lambda_z = NA, t_half = NA, r_squared = NA, n_points = 0))
  }

  # Find Cmax index and work with post-Cmax data
  idx_cmax <- which.max(conc)
  if (idx_cmax >= n - 1) {
    # Not enough post-Cmax points
    return(list(lambda_z = NA, t_half = NA, r_squared = NA, n_points = 0))
  }

  # Use only positive concentrations after Cmax
  post_idx <- (idx_cmax):n
  post_time <- time[post_idx]
  post_conc <- conc[post_idx]
  valid <- post_conc > 0
  post_time <- post_time[valid]
  post_conc <- post_conc[valid]
  n_post <- length(post_time)

  if (n_post < 3) {
    return(list(lambda_z = NA, t_half = NA, r_squared = NA, n_points = 0))
  }

  log_conc <- log(post_conc)

  # Try different numbers of terminal points (minimum 3)
  best_r2 <- -Inf
  best_result <- NULL

  for (np in 3:n_post) {
    idx <- (n_post - np + 1):n_post
    t_sub <- post_time[idx]
    lc_sub <- log_conc[idx]

    fit <- lm(lc_sub ~ t_sub)
    slope <- coef(fit)[2]
    r2 <- summary(fit)$r.squared

    # lambda_z must be positive (slope must be negative)
    if (!is.na(slope) && slope < 0 && r2 > best_r2) {
      best_r2 <- r2
      best_result <- list(
        lambda_z = -slope,
        t_half = log(2) / (-slope),
        r_squared = r2,
        n_points = np,
        intercept = coef(fit)[1]
      )
    }
  }

  if (is.null(best_result)) {
    return(list(lambda_z = NA, t_half = NA, r_squared = NA, n_points = 0))
  }
  return(best_result)
}

#' Run non-compartmental analysis for a single subject
#'
#' @param time Numeric vector of time points
#' @param conc Numeric vector of concentrations
#' @param dose Numeric dose administered
#' @param route Character: "IV" or "EV" (extravascular)
#' @return Named list of NCA parameters
nca_single_subject <- function(time, conc, dose, route) {
  # Sort by time

  ord <- order(time)
  time <- time[ord]
  conc <- conc[ord]

  # Basic parameters
  cmax <- max(conc, na.rm = TRUE)
  tmax <- time[which.max(conc)]
  c_last <- conc[length(conc)]
  t_last <- time[length(time)]

  # AUC_last (to last measurable concentration)
  # Find last positive concentration
  last_pos <- max(which(conc > 0))
  auc_last <- auc_lin_log(time[1:last_pos], conc[1:last_pos])

  # Terminal phase
  term <- estimate_lambda_z(time, conc)
  lambda_z <- term$lambda_z
  t_half <- term$t_half

  # AUC_inf
  if (!is.na(lambda_z) && conc[last_pos] > 0) {
    auc_extrap <- conc[last_pos] / lambda_z
    auc_inf <- auc_last + auc_extrap
    pct_extrap <- (auc_extrap / auc_inf) * 100
  } else {
    auc_inf <- NA
    pct_extrap <- NA
  }

  # AUMC_last (area under first moment curve)
  aumc_last <- auc_lin_log(time[1:last_pos], time[1:last_pos] * conc[1:last_pos])

  # MRT (mean residence time)
  if (!is.na(auc_inf) && auc_inf > 0 && !is.na(lambda_z)) {
    aumc_inf <- aumc_last + (t_last * c_last / lambda_z) + (c_last / lambda_z^2)
    mrt <- aumc_inf / auc_inf
    if (toupper(route) == "IV") {
      mrt <- mrt  # MRT_iv
    }
  } else {
    aumc_inf <- NA
    mrt <- NA
  }

  # Clearance and Volume
  if (!is.na(auc_inf) && auc_inf > 0) {
    cl <- dose / auc_inf
    if (toupper(route) != "IV") {
      cl <- cl  # CL/F for extravascular
    }
  } else {
    cl <- NA
  }

  if (!is.na(cl) && !is.na(lambda_z) && lambda_z > 0) {
    vd <- cl / lambda_z  # Vd or Vd/F
  } else {
    vd <- NA
  }

  # Vss
  if (!is.na(cl) && !is.na(mrt)) {
    vss <- cl * mrt
  } else {
    vss <- NA
  }

  params <- list(
    Cmax = cmax,
    Tmax = tmax,
    AUC_last = auc_last,
    AUC_inf = auc_inf,
    AUC_extrap_pct = pct_extrap,
    Lambda_z = lambda_z,
    Half_life = t_half,
    Lambda_z_R2 = term$r_squared,
    Lambda_z_npoints = term$n_points,
    MRT = mrt,
    CL = cl,
    Vd = vd,
    Vss = vss
  )

  return(params)
}

#' Run NCA on population data (all subjects)
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param dose Numeric dose
#' @param route Character route of administration
#' @return Data frame with NCA parameters per subject, plus summary row
run_nca <- function(df, dose, route) {
  subjects <- unique(df$ID)
  results <- list()

  for (subj in subjects) {
    subj_data <- df[df$ID == subj, ]
    params <- nca_single_subject(subj_data$Time, subj_data$Conc, dose, route)
    params$ID <- as.character(subj)
    results[[as.character(subj)]] <- params
  }

  # Combine into data frame
  param_names <- c("Cmax", "Tmax", "AUC_last", "AUC_inf", "AUC_extrap_pct",
                    "Lambda_z", "Half_life", "Lambda_z_R2", "Lambda_z_npoints",
                    "MRT", "CL", "Vd", "Vss")

  result_df <- data.frame(
    ID = sapply(results, function(x) x$ID),
    stringsAsFactors = FALSE
  )

  for (pname in param_names) {
    result_df[[pname]] <- sapply(results, function(x) {
      val <- x[[pname]]
      if (is.null(val) || length(val) == 0) NA_real_ else as.numeric(val)
    })
  }
  rownames(result_df) <- NULL

  return(result_df)
}

#' Compute summary statistics for NCA results
#'
#' @param nca_df Data frame from run_nca
#' @return Data frame with Mean, SD, SEM, Median, Min, Max per parameter
nca_summary <- function(nca_df) {
  param_cols <- setdiff(names(nca_df), c("ID", "Lambda_z_npoints"))

  summary_list <- list()
  for (pc in param_cols) {
    vals <- nca_df[[pc]]
    vals <- vals[!is.na(vals)]
    n <- length(vals)
    if (n > 0) {
      summary_list[[pc]] <- data.frame(
        Parameter = pc,
        N = n,
        Mean = mean(vals),
        SD = sd(vals),
        SEM = if (n > 1) sd(vals) / sqrt(n) else NA,
        Median = median(vals),
        Min = min(vals),
        Max = max(vals),
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, summary_list)
}

#' Get NCA parameter labels with units
#'
#' @param route Character route of administration
#' @param time_unit Character time unit
#' @param conc_unit Character concentration unit
#' @return Named vector of labels
nca_param_labels <- function(route = "IV", time_unit = "h", conc_unit = "ng/mL") {
  cl_label <- if (toupper(route) == "IV") "CL" else "CL/F"
  vd_label <- if (toupper(route) == "IV") "Vd" else "Vd/F"
  vss_label <- if (toupper(route) == "IV") "Vss" else "Vss/F"

  c(
    Cmax = paste0("Cmax (", conc_unit, ")"),
    Tmax = paste0("Tmax (", time_unit, ")"),
    AUC_last = paste0("AUC_last (", conc_unit, "\u00B7", time_unit, ")"),
    AUC_inf = paste0("AUC_inf (", conc_unit, "\u00B7", time_unit, ")"),
    AUC_extrap_pct = "AUC extrapolated (%)",
    Lambda_z = paste0("\u03BBz (1/", time_unit, ")"),
    Half_life = paste0("t\u00BD (", time_unit, ")"),
    Lambda_z_R2 = "Terminal R\u00B2",
    MRT = paste0("MRT (", time_unit, ")"),
    CL = paste0(cl_label, " (L/", time_unit, ")"),
    Vd = paste0(vd_label, " (L)"),
    Vss = paste0(vss_label, " (L)")
  )
}
