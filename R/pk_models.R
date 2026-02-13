# Compartmental PK Models - Non-Linear Mixed Effects Modeling Module

library(nlme)

# ============================================================
# Model Functions
# ============================================================

#' One-compartment IV bolus model function
#' C(t) = (Dose/V) * exp(-k * t)
#' Parameters: lV (log-Volume), lk (log-elimination rate constant)
one_comp_iv <- function(time, dose, lV, lk) {
  V <- exp(lV)
  k <- exp(lk)
  (dose / V) * exp(-k * time)
}

#' One-compartment extravascular (first-order absorption) model
#' C(t) = (Dose * ka / (V * (ka - k))) * (exp(-k*t) - exp(-ka*t))
#' Parameters: lV, lk, lka (log-absorption rate constant)
one_comp_ev <- function(time, dose, lV, lk, lka) {
  V <- exp(lV)
  k <- exp(lk)
  ka <- exp(lka)
  if (abs(ka - k) < 1e-10) {
    # Flip-flop: use limit form
    (dose / V) * k * time * exp(-k * time)
  } else {
    (dose * ka / (V * (ka - k))) * (exp(-k * time) - exp(-ka * time))
  }
}

#' Two-compartment IV bolus model
#' C(t) = A*exp(-alpha*t) + B*exp(-beta*t)
#' Parameterized as: CL, V1, Q, V2
two_comp_iv <- function(time, dose, lCL, lV1, lQ, lV2) {
  CL <- exp(lCL)
  V1 <- exp(lV1)
  Q  <- exp(lQ)
  V2 <- exp(lV2)

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  beta  <- 0.5 * ((k10 + k12 + k21) - sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  alpha <- 0.5 * ((k10 + k12 + k21) + sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))

  A <- (dose / V1) * (alpha - k21) / (alpha - beta)
  B <- (dose / V1) * (k21 - beta) / (alpha - beta)

  A * exp(-alpha * time) + B * exp(-beta * time)
}

#' Two-compartment extravascular model
two_comp_ev <- function(time, dose, lCL, lV1, lQ, lV2, lka) {
  CL <- exp(lCL)
  V1 <- exp(lV1)
  Q  <- exp(lQ)
  V2 <- exp(lV2)
  ka <- exp(lka)

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2

  beta  <- 0.5 * ((k10 + k12 + k21) - sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))
  alpha <- 0.5 * ((k10 + k12 + k21) + sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21))

  A_coeff <- (dose * ka / V1) * (k21 - alpha) / ((ka - alpha) * (beta - alpha))
  B_coeff <- (dose * ka / V1) * (k21 - beta) / ((ka - beta) * (alpha - beta))
  C_coeff <- (dose * ka / V1) * (k21 - ka) / ((alpha - ka) * (beta - ka))

  A_coeff * exp(-alpha * time) + B_coeff * exp(-beta * time) + C_coeff * exp(-ka * time)
}

#' Three-compartment IV bolus model
#' C(t) = A*exp(-alpha*t) + B*exp(-beta*t) + C*exp(-gamma*t)
three_comp_iv <- function(time, dose, lCL, lV1, lQ2, lV2, lQ3, lV3) {
  CL <- exp(lCL)
  V1 <- exp(lV1)
  Q2 <- exp(lQ2)
  V2 <- exp(lV2)
  Q3 <- exp(lQ3)
  V3 <- exp(lV3)

  k10 <- CL / V1
  k12 <- Q2 / V1
  k21 <- Q2 / V2
  k13 <- Q3 / V1
  k31 <- Q3 / V3

  a0 <- k10 * k21 * k31
  a1 <- k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
  a2 <- k10 + k12 + k13 + k21 + k31

  p <- a1 - (a2^2) / 3
  q <- (2 * a2^3 / 27) - (a1 * a2 / 3) + a0
  r1 <- sqrt(-p^3 / 27)
  phi <- acos(-q / (2 * r1)) / 3

  r1_cbrt <- r1^(1/3)
  alpha <- -(cos(phi) * 2 * r1_cbrt - a2 / 3)
  beta  <- -(cos(phi + 2 * pi / 3) * 2 * r1_cbrt - a2 / 3)
  gamma <- -(cos(phi + 4 * pi / 3) * 2 * r1_cbrt - a2 / 3)

  # Sort so alpha > beta > gamma
  rates <- sort(c(alpha, beta, gamma), decreasing = TRUE)
  alpha <- rates[1]
  beta <- rates[2]
  gamma <- rates[3]

  A <- (dose / V1) * ((k21 - alpha) * (k31 - alpha)) / ((beta - alpha) * (gamma - alpha))
  B <- (dose / V1) * ((k21 - beta) * (k31 - beta)) / ((alpha - beta) * (gamma - beta))
  C <- (dose / V1) * ((k21 - gamma) * (k31 - gamma)) / ((alpha - gamma) * (beta - gamma))

  A * exp(-alpha * time) + B * exp(-beta * time) + C * exp(-gamma * time)
}

#' Three-compartment extravascular model
three_comp_ev <- function(time, dose, lCL, lV1, lQ2, lV2, lQ3, lV3, lka) {
  CL <- exp(lCL)
  V1 <- exp(lV1)
  Q2 <- exp(lQ2)
  V2 <- exp(lV2)
  Q3 <- exp(lQ3)
  V3 <- exp(lV3)
  ka <- exp(lka)

  k10 <- CL / V1
  k12 <- Q2 / V1
  k21 <- Q2 / V2
  k13 <- Q3 / V1
  k31 <- Q3 / V3

  a0 <- k10 * k21 * k31
  a1 <- k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
  a2 <- k10 + k12 + k13 + k21 + k31

  p <- a1 - (a2^2) / 3
  q <- (2 * a2^3 / 27) - (a1 * a2 / 3) + a0
  r1 <- sqrt(-p^3 / 27)
  phi <- acos(-q / (2 * r1)) / 3

  r1_cbrt <- r1^(1/3)
  alpha <- -(cos(phi) * 2 * r1_cbrt - a2 / 3)
  beta  <- -(cos(phi + 2 * pi / 3) * 2 * r1_cbrt - a2 / 3)
  gamma <- -(cos(phi + 4 * pi / 3) * 2 * r1_cbrt - a2 / 3)

  rates <- sort(c(alpha, beta, gamma), decreasing = TRUE)
  alpha <- rates[1]
  beta <- rates[2]
  gamma <- rates[3]

  # Coefficients including absorption
  denom_a <- (beta - alpha) * (gamma - alpha) * (ka - alpha)
  denom_b <- (alpha - beta) * (gamma - beta) * (ka - beta)
  denom_g <- (alpha - gamma) * (beta - gamma) * (ka - gamma)
  denom_k <- (alpha - ka) * (beta - ka) * (gamma - ka)

  coeff <- dose * ka / V1

  A <- coeff * ((k21 - alpha) * (k31 - alpha)) / denom_a
  B <- coeff * ((k21 - beta) * (k31 - beta)) / denom_b
  C <- coeff * ((k21 - gamma) * (k31 - gamma)) / denom_g
  D <- coeff * ((k21 - ka) * (k31 - ka)) / denom_k

  A * exp(-alpha * time) + B * exp(-beta * time) + C * exp(-gamma * time) + D * exp(-ka * time)
}


# ============================================================
# Initial Parameter Estimation
# ============================================================

#' Estimate initial parameters from data for a given model
#' Uses method of residuals (curve stripping) for multi-compartment IV models
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param dose Numeric dose
#' @param model_type Character: "1comp", "2comp", "3comp"
#' @param route Character: "IV" or "EV"
#' @return Named list of starting parameter values (on log scale)
estimate_initial_params <- function(df, dose, model_type, route) {
  mean_data <- aggregate(Conc ~ Time, data = df, FUN = mean)
  mean_data <- mean_data[order(mean_data$Time), ]

  time <- mean_data$Time
  conc <- mean_data$Conc

  cmax <- max(conc)
  tmax_idx <- which.max(conc)
  is_iv <- toupper(route) == "IV"

  # --- Curve stripping for IV 2-comp and 3-comp ---
  if (is_iv && model_type %in% c("2comp", "3comp")) {
    params <- tryCatch(curve_strip_iv(time, conc, dose, model_type), error = function(e) NULL)
    if (!is.null(params)) {
      if (model_type == "2comp" && !is_iv) {
        ka_est <- max(exp(params$lCL) / exp(params$lV1) * 3, 0.5)
        params$lka <- log(ka_est)
      }
      return(params)
    }
  }

  # --- Fallback: simple estimation ---
  if (is_iv && conc[1] > 0) {
    V_est <- dose / conc[1]
  } else {
    V_est <- dose / cmax
  }
  V_est <- max(V_est, 0.01)

  # Terminal slope from last portion of data
  n_pts <- length(conc)
  n_term <- max(3, floor(n_pts / 2))
  term_conc <- conc[(n_pts - n_term + 1):n_pts]
  term_time <- time[(n_pts - n_term + 1):n_pts]
  pos_idx <- term_conc > 0
  if (sum(pos_idx) >= 2) {
    log_conc <- log(term_conc[pos_idx])
    t_sub <- term_time[pos_idx]
    fit <- lm(log_conc ~ t_sub)
    k_est <- max(-coef(fit)[2], 0.01)
  } else {
    k_est <- 0.1
  }

  CL_est <- k_est * V_est

  if (model_type == "1comp") {
    params <- list(lV = log(V_est), lk = log(k_est))
    if (!is_iv) {
      ka_est <- max(k_est * 3, 0.5)
      params$lka <- log(ka_est)
    }
  } else if (model_type == "2comp") {
    params <- list(
      lCL = log(CL_est),
      lV1 = log(V_est * 0.4),
      lQ  = log(CL_est * 0.5),
      lV2 = log(V_est * 0.6)
    )
    if (!is_iv) {
      ka_est <- max(k_est * 3, 0.5)
      params$lka <- log(ka_est)
    }
  } else if (model_type == "3comp") {
    params <- list(
      lCL = log(CL_est),
      lV1 = log(V_est * 0.3),
      lQ2 = log(CL_est * 0.4),
      lV2 = log(V_est * 0.4),
      lQ3 = log(CL_est * 0.1),
      lV3 = log(V_est * 0.3)
    )
    if (!is_iv) {
      ka_est <- max(k_est * 3, 0.5)
      params$lka <- log(ka_est)
    }
  }

  return(params)
}

#' Curve stripping (method of residuals) for IV multi-compartment models
#' Peels exponential terms from the terminal phase back to get alpha, beta, A, B
curve_strip_iv <- function(time, conc, dose, model_type) {
  pos <- conc > 0
  time <- time[pos]
  conc <- conc[pos]
  n <- length(time)

  # Terminal phase: fit log-linear to last ~half of data
  n_term <- max(3, floor(n * 0.5))
  idx_term <- (n - n_term + 1):n
  fit_term <- lm(log(conc[idx_term]) ~ time[idx_term])
  beta <- max(-coef(fit_term)[2], 1e-6)
  B_intercept <- exp(coef(fit_term)[1])

  # Residuals: subtract terminal component from all data
  resid_conc <- conc - B_intercept * exp(-beta * time)

  if (model_type == "2comp") {
    # Fit the residual (distribution phase)
    pos_resid <- resid_conc > 0
    if (sum(pos_resid) >= 2) {
      fit_dist <- lm(log(resid_conc[pos_resid]) ~ time[pos_resid])
      alpha <- max(-coef(fit_dist)[2], beta + 0.01)
      A_intercept <- exp(coef(fit_dist)[1])
    } else {
      alpha <- beta * 5
      A_intercept <- conc[1] - B_intercept
    }
    A_intercept <- max(A_intercept, 0.01)

    # Back-calculate CL, V1, Q, V2 from A, B, alpha, beta
    C0 <- A_intercept + B_intercept
    V1 <- dose / C0
    k21 <- (A_intercept * beta + B_intercept * alpha) / C0
    k10 <- alpha * beta / k21
    k12 <- alpha + beta - k21 - k10
    k12 <- max(k12, 0.001)
    k21 <- max(k21, 0.001)
    CL <- k10 * V1
    Q  <- k12 * V1
    V2 <- Q / k21

    return(list(
      lCL = log(max(CL, 1e-6)),
      lV1 = log(max(V1, 1e-6)),
      lQ  = log(max(Q, 1e-6)),
      lV2 = log(max(V2, 1e-6))
    ))
  } else {
    # 3-comp: peel a second residual
    resid2 <- resid_conc
    pos_r2 <- resid2 > 0
    if (sum(pos_r2) >= 3) {
      n_mid <- max(2, floor(sum(pos_r2) * 0.5))
      idx_mid <- which(pos_r2)
      idx_mid <- idx_mid[(length(idx_mid) - n_mid + 1):length(idx_mid)]
      fit_mid <- lm(log(resid2[idx_mid]) ~ time[idx_mid])
      alpha2 <- max(-coef(fit_mid)[2], beta + 0.1)
      B2_int <- exp(coef(fit_mid)[1])
    } else {
      alpha2 <- beta * 3
      B2_int <- conc[1] * 0.3
    }

    resid3 <- resid2 - B2_int * exp(-alpha2 * time)
    pos_r3 <- resid3 > 0
    if (sum(pos_r3) >= 2) {
      fit_fast <- lm(log(resid3[pos_r3]) ~ time[pos_r3])
      alpha1 <- max(-coef(fit_fast)[2], alpha2 + 0.1)
      A1_int <- exp(coef(fit_fast)[1])
    } else {
      alpha1 <- alpha2 * 5
      A1_int <- conc[1] * 0.4
    }

    C0 <- A1_int + B2_int + B_intercept
    V1 <- dose / C0
    CL <- V1 * alpha1 * alpha2 * beta / (alpha2 * beta + alpha1 * beta + alpha1 * alpha2 -
           alpha1 * alpha1 - beta * beta - alpha2 * alpha2 + V1)
    CL <- max(abs(CL), 1e-6)

    return(list(
      lCL = log(CL),
      lV1 = log(max(V1, 1e-6)),
      lQ2 = log(max(CL * 0.5, 1e-6)),
      lV2 = log(max(V1 * 0.8, 1e-6)),
      lQ3 = log(max(CL * 0.15, 1e-6)),
      lV3 = log(max(V1 * 0.5, 1e-6))
    ))
  }
}


# ============================================================
# NLME Fitting
# ============================================================

#' Fit a compartmental model using nlme with robust multi-level fallback
#'
#' Strategy:
#' 1. Try per-subject nls fits to get better starting values
#' 2. Try nlme with progressively simpler random-effects structures
#' 3. Try nlme with varPower() error model (PK data often has proportional error)
#' 4. Fall back to gnls (no random effects) as last resort
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param dose Numeric dose
#' @param model_type Character: "1comp", "2comp", "3comp"
#' @param route Character: "IV" or "EV"
#' @return List with 'fit' (nlme object or NULL), 'summary', 'params', 'error'
fit_compartmental_model <- function(df, dose, model_type, route) {
  route_upper <- toupper(route)
  is_iv <- route_upper == "IV"

  # Prepare data
  pk_data <- df[, c("ID", "Time", "Conc")]
  pk_data$Dose <- dose

  # Get initial estimates (uses curve stripping for IV multi-compartment)
  init_params <- estimate_initial_params(df, dose, model_type, route)

  # Try to refine starting values with per-subject nls fits
  refined <- tryCatch(
    refine_inits_with_nls(pk_data, init_params, model_type, is_iv),
    error = function(e) NULL
  )
  if (!is.null(refined)) {
    init_params <- refined
  }

  # Build list of fitting attempts (progressively more relaxed)
  fit <- NULL
  last_error <- ""
  used_gnls <- FALSE

  # Define control levels: progressively more relaxed
  ctrl_strict <- nlmeControl(maxIter = 200, pnlsTol = 0.01, msMaxIter = 200,
                              tolerance = 1e-5, returnObject = TRUE)
  ctrl_medium <- nlmeControl(maxIter = 300, pnlsTol = 0.1, msMaxIter = 300,
                              tolerance = 1e-4, returnObject = TRUE)
  ctrl_relaxed <- nlmeControl(maxIter = 500, pnlsTol = 1.0, msMaxIter = 500,
                               tolerance = 1e-3, returnObject = TRUE)
  ctrl_very_relaxed <- nlmeControl(maxIter = 500, pnlsTol = 2.0, msMaxIter = 500,
                                    tolerance = 1e-2, returnObject = TRUE, opt = "nlm")

  # --- Phase 1: standard nlme attempts ---
  attempts <- build_fitting_attempts(model_type, is_iv,
                                      ctrl_strict, ctrl_medium, ctrl_relaxed, ctrl_very_relaxed)

  for (attempt in attempts) {
    fit <- tryCatch({
      do_nlme_fit(pk_data, init_params, model_type, is_iv, attempt)
    }, error = function(e) {
      last_error <<- e$message
      NULL
    })
    if (!is.null(fit)) break
  }

  # --- Phase 2: nlme with varPower() error model ---
  if (is.null(fit)) {
    for (re_type in c("re_cl_only", "re_v1_only", "re_v_only", "re_k_only")) {
      # skip irrelevant re types for the model
      if (model_type == "1comp" && re_type %in% c("re_cl_only", "re_v1_only")) next
      if (model_type != "1comp" && re_type %in% c("re_v_only", "re_k_only")) next

      fit <- tryCatch({
        do_nlme_fit_varpower(pk_data, init_params, model_type, is_iv, re_type, ctrl_relaxed)
      }, error = function(e) {
        last_error <<- e$message
        NULL
      })
      if (!is.null(fit)) break
    }
  }

  # --- Phase 3: gnls fallback (no random effects, population-level only) ---
  if (is.null(fit)) {
    fit <- tryCatch({
      do_gnls_fit(pk_data, init_params, model_type, is_iv)
    }, error = function(e) {
      last_error <<- e$message
      NULL
    })
    if (!is.null(fit)) used_gnls <- TRUE
  }

  if (is.null(fit)) {
    return(list(
      fit = NULL, params = NULL, summary = NULL, predictions = NULL,
      error = paste("Model fitting failed after all attempts:", last_error,
                    "\nTry a simpler model or check your data quality.")
    ))
  }

  # Extract parameters
  tryCatch({
    if (used_gnls) {
      params <- extract_gnls_params(fit, model_type, route, dose, pk_data)
    } else {
      params <- extract_compartmental_params(fit, model_type, route, dose)
    }
    preds <- if (used_gnls) {
      generate_gnls_predictions(fit, pk_data, model_type, route, dose)
    } else {
      generate_predictions(fit, pk_data, model_type, route, dose)
    }
    return(list(
      fit = fit,
      params = params$individual,
      summary = params$population,
      predictions = preds,
      error = NULL
    ))
  }, error = function(e) {
    return(list(
      fit = NULL, params = NULL, summary = NULL, predictions = NULL,
      error = paste("Parameter extraction failed:", e$message,
                    "\nThe model may have converged to invalid estimates. Try a simpler model.")
    ))
  })
}

#' Refine starting values by fitting nls to each subject individually
refine_inits_with_nls <- function(data, inits, model_type, is_iv) {
  subjects <- unique(data$ID)
  all_params <- list()

  pred_fn <- get_pred_function_nls(model_type, is_iv)
  start_list <- as.list(unlist(inits))

  for (subj in subjects) {
    sdata <- data[data$ID == subj, ]
    subj_fit <- tryCatch({
      nls(Conc ~ pred_fn(params, Time, Dose),
          data = sdata,
          start = list(params = unlist(start_list)),
          algorithm = "port",
          lower = rep(-20, length(start_list)),
          upper = rep(20, length(start_list)),
          control = nls.control(maxiter = 100, warnOnly = TRUE))
    }, error = function(e) NULL)

    if (!is.null(subj_fit)) {
      all_params[[length(all_params) + 1]] <- coef(subj_fit)
    }
  }

  if (length(all_params) < 2) return(NULL)

  # Use median of per-subject estimates
  param_mat <- do.call(rbind, all_params)
  median_params <- apply(param_mat, 2, median)
  as.list(setNames(median_params, names(inits)))
}

#' Get a prediction function wrapper for nls fitting
get_pred_function_nls <- function(model_type, is_iv) {
  if (model_type == "1comp" && is_iv) {
    function(p, Time, Dose) (Dose / exp(p[1])) * exp(-exp(p[2]) * Time)
  } else if (model_type == "1comp" && !is_iv) {
    function(p, Time, Dose) .pred_1comp_ev(p[1], p[2], p[3], Time, Dose)
  } else if (model_type == "2comp" && is_iv) {
    function(p, Time, Dose) .pred_2comp_iv(p[1], p[2], p[3], p[4], Time, Dose)
  } else if (model_type == "2comp" && !is_iv) {
    function(p, Time, Dose) .pred_2comp_ev(p[1], p[2], p[3], p[4], p[5], Time, Dose)
  } else if (model_type == "3comp" && is_iv) {
    function(p, Time, Dose) .pred_3comp_iv(p[1], p[2], p[3], p[4], p[5], p[6], Time, Dose)
  } else {
    function(p, Time, Dose) .pred_3comp_ev(p[1], p[2], p[3], p[4], p[5], p[6], p[7], Time, Dose)
  }
}

#' Build list of fitting attempts with progressively simpler structures
build_fitting_attempts <- function(model_type, is_iv, ctrl_s, ctrl_m, ctrl_r, ctrl_vr) {
  if (model_type == "1comp") {
    return(list(
      list(re = "full_diag", ctrl = ctrl_s),
      list(re = "full_sym",  ctrl = ctrl_m),
      list(re = "re_v_only", ctrl = ctrl_r),
      list(re = "re_k_only", ctrl = ctrl_r)
    ))
  } else if (model_type == "2comp") {
    return(list(
      list(re = "full_diag",   ctrl = ctrl_s),     # pdDiag on CL+V1
      list(re = "full_sym",    ctrl = ctrl_m),      # pdSymm on CL+V1
      list(re = "re_cl_only",  ctrl = ctrl_m),      # CL only
      list(re = "re_v1_only",  ctrl = ctrl_m),      # V1 only
      list(re = "re_cl_only",  ctrl = ctrl_r),      # CL only, relaxed
      list(re = "re_v1_only",  ctrl = ctrl_r),      # V1 only, relaxed
      list(re = "re_cl_only",  ctrl = ctrl_vr)      # CL only, very relaxed
    ))
  } else {
    # 3-comp
    return(list(
      list(re = "full_diag",   ctrl = ctrl_s),
      list(re = "full_sym",    ctrl = ctrl_m),
      list(re = "re_cl_only",  ctrl = ctrl_m),
      list(re = "re_v1_only",  ctrl = ctrl_m),
      list(re = "re_cl_only",  ctrl = ctrl_r),
      list(re = "re_v1_only",  ctrl = ctrl_r),
      list(re = "re_cl_only",  ctrl = ctrl_vr)
    ))
  }
}

#' Execute a single nlme fitting attempt
do_nlme_fit <- function(data, inits, model_type, is_iv, attempt) {
  re_type <- attempt$re
  ctrl <- attempt$ctrl

  if (model_type == "1comp") {
    if (is_iv) {
      model_formula <- Conc ~ (Dose / exp(lV)) * exp(-exp(lk) * Time)
      fixed_form <- lV + lk ~ 1
      starts <- c(lV = inits$lV, lk = inits$lk)
      random_form <- switch(re_type,
        "full_diag" = lV + lk ~ 1 | ID,
        "full_sym"  = lV + lk ~ 1 | ID,
        "re_v_only" = lV ~ 1 | ID,
        "re_k_only" = lk ~ 1 | ID
      )
    } else {
      model_formula <- Conc ~ .pred_1comp_ev(lV, lk, lka, Time, Dose)
      fixed_form <- lV + lk + lka ~ 1
      starts <- c(lV = inits$lV, lk = inits$lk, lka = inits$lka)
      random_form <- switch(re_type,
        "full_diag" = lV + lk ~ 1 | ID,
        "full_sym"  = lV + lk ~ 1 | ID,
        "re_v_only" = lV ~ 1 | ID,
        "re_k_only" = lk ~ 1 | ID
      )
    }
    use_diag <- re_type == "full_diag"
  } else if (model_type == "2comp") {
    if (is_iv) {
      model_formula <- Conc ~ .pred_2comp_iv(lCL, lV1, lQ, lV2, Time, Dose)
      fixed_form <- lCL + lV1 + lQ + lV2 ~ 1
      starts <- c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ, lV2 = inits$lV2)
    } else {
      model_formula <- Conc ~ .pred_2comp_ev(lCL, lV1, lQ, lV2, lka, Time, Dose)
      fixed_form <- lCL + lV1 + lQ + lV2 + lka ~ 1
      starts <- c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ,
                   lV2 = inits$lV2, lka = inits$lka)
    }
    random_form <- switch(re_type,
      "full_diag"  = lCL + lV1 ~ 1 | ID,
      "full_sym"   = lCL + lV1 ~ 1 | ID,
      "re_cl_only" = lCL ~ 1 | ID,
      "re_v1_only" = lV1 ~ 1 | ID
    )
    use_diag <- re_type == "full_diag"
  } else {
    # 3-comp
    if (is_iv) {
      model_formula <- Conc ~ .pred_3comp_iv(lCL, lV1, lQ2, lV2, lQ3, lV3, Time, Dose)
      fixed_form <- lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 ~ 1
      starts <- c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
                   lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3)
    } else {
      model_formula <- Conc ~ .pred_3comp_ev(lCL, lV1, lQ2, lV2, lQ3, lV3, lka, Time, Dose)
      fixed_form <- lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 + lka ~ 1
      starts <- c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
                   lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3, lka = inits$lka)
    }
    random_form <- switch(re_type,
      "full_diag"  = lCL + lV1 ~ 1 | ID,
      "full_sym"   = lCL + lV1 ~ 1 | ID,
      "re_cl_only" = lCL ~ 1 | ID,
      "re_v1_only" = lV1 ~ 1 | ID
    )
    use_diag <- re_type == "full_diag"
  }

  if (use_diag) {
    nlme(model_formula,
         data = data,
         fixed = fixed_form,
         random = pdDiag(random_form),
         start = starts,
         control = ctrl,
         na.action = na.omit)
  } else {
    nlme(model_formula,
         data = data,
         fixed = fixed_form,
         random = random_form,
         start = starts,
         control = ctrl,
         na.action = na.omit)
  }
}

#' Execute nlme fit with varPower error model (better for heteroscedastic PK data)
do_nlme_fit_varpower <- function(data, inits, model_type, is_iv, re_type, ctrl) {
  specs <- get_model_specs(inits, model_type, is_iv)

  random_form <- switch(re_type,
    "re_cl_only" = lCL ~ 1 | ID,
    "re_v1_only" = lV1 ~ 1 | ID,
    "re_v_only"  = lV ~ 1 | ID,
    "re_k_only"  = lk ~ 1 | ID
  )

  nlme(specs$formula,
       data = data,
       fixed = specs$fixed,
       random = pdDiag(random_form),
       start = specs$starts,
       weights = varPower(form = ~fitted(.)),
       control = ctrl,
       na.action = na.omit)
}

#' Fit using gnls (no random effects - population-level only, last resort)
do_gnls_fit <- function(data, inits, model_type, is_iv) {
  specs <- get_model_specs(inits, model_type, is_iv)

  ctrl <- gnlsControl(maxIter = 500, nlsTol = 0.01, tolerance = 1e-4,
                       returnObject = TRUE)

  # Try with default error, then varPower
  fit <- tryCatch(
    gnls(specs$formula, data = data, start = specs$starts,
         control = ctrl, na.action = na.omit),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    fit <- gnls(specs$formula, data = data, start = specs$starts,
                weights = varPower(form = ~fitted(.)),
                control = ctrl, na.action = na.omit)
  }

  return(fit)
}

#' Helper to get model formula, fixed form, and start values
get_model_specs <- function(inits, model_type, is_iv) {
  if (model_type == "1comp" && is_iv) {
    list(formula = Conc ~ (Dose / exp(lV)) * exp(-exp(lk) * Time),
         fixed = lV + lk ~ 1,
         starts = c(lV = inits$lV, lk = inits$lk))
  } else if (model_type == "1comp" && !is_iv) {
    list(formula = Conc ~ .pred_1comp_ev(lV, lk, lka, Time, Dose),
         fixed = lV + lk + lka ~ 1,
         starts = c(lV = inits$lV, lk = inits$lk, lka = inits$lka))
  } else if (model_type == "2comp" && is_iv) {
    list(formula = Conc ~ .pred_2comp_iv(lCL, lV1, lQ, lV2, Time, Dose),
         fixed = lCL + lV1 + lQ + lV2 ~ 1,
         starts = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ, lV2 = inits$lV2))
  } else if (model_type == "2comp" && !is_iv) {
    list(formula = Conc ~ .pred_2comp_ev(lCL, lV1, lQ, lV2, lka, Time, Dose),
         fixed = lCL + lV1 + lQ + lV2 + lka ~ 1,
         starts = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ,
                    lV2 = inits$lV2, lka = inits$lka))
  } else if (model_type == "3comp" && is_iv) {
    list(formula = Conc ~ .pred_3comp_iv(lCL, lV1, lQ2, lV2, lQ3, lV3, Time, Dose),
         fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 ~ 1,
         starts = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
                    lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3))
  } else {
    list(formula = Conc ~ .pred_3comp_ev(lCL, lV1, lQ2, lV2, lQ3, lV3, lka, Time, Dose),
         fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 + lka ~ 1,
         starts = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
                    lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3, lka = inits$lka))
  }
}

#' Extract parameters from gnls fit (no random effects, so no individual estimates)
extract_gnls_params <- function(fit, model_type, route, dose, data) {
  fe <- coef(fit)
  is_iv <- toupper(route) == "IV"

  vcov_mat <- tryCatch(vcov(fit), error = function(e) NULL)
  df_resid <- tryCatch(nrow(data) - length(fe), error = function(e) NULL)

  if (model_type == "1comp") {
    pop_params <- extract_1comp_params(fe, dose, is_iv, vcov_mat, df_resid)
  } else if (model_type == "2comp") {
    pop_params <- extract_2comp_params(fe, dose, is_iv, vcov_mat, df_resid)
  } else {
    pop_params <- extract_3comp_params(fe, dose, is_iv, vcov_mat, df_resid)
  }

  pop_params$AIC <- AIC(fit)
  pop_params$BIC <- BIC(fit)
  pop_params$logLik <- as.numeric(logLik(fit))

  # For gnls, individual params = population params for each subject (wide format)
  subjects <- unique(data$ID)
  extract_fn <- switch(model_type,
    "1comp" = extract_1comp_params,
    "2comp" = extract_2comp_params,
    "3comp" = extract_3comp_params
  )
  indiv_list <- list()
  for (subj in subjects) {
    ip <- extract_fn(fe, dose, is_iv)
    row_df <- as.data.frame(as.list(setNames(ip$Estimate, ip$Parameter)),
                            check.names = FALSE, stringsAsFactors = FALSE)
    row_df$ID <- as.character(subj)
    indiv_list[[length(indiv_list) + 1]] <- row_df
  }
  indiv_df <- do.call(rbind, indiv_list)

  return(list(population = pop_params, individual = indiv_df))
}

#' Generate predictions from gnls fit
generate_gnls_predictions <- function(fit, data, model_type, route, dose) {
  time_range <- range(data$Time)
  time_seq <- seq(time_range[1], time_range[2], length.out = 200)

  fe <- coef(fit)
  is_iv <- toupper(route) == "IV"

  pop_pred <- predict_from_params(fe, time_seq, dose, model_type, is_iv)
  pred_df <- data.frame(Time = time_seq, Pred_pop = pop_pred)

  # For gnls, individual predictions = population predictions
  subjects <- unique(data$ID)
  indiv_preds <- list()
  for (subj in subjects) {
    indiv_preds[[as.character(subj)]] <- data.frame(
      Time = time_seq,
      Pred_indiv = pop_pred,
      ID = as.character(subj),
      stringsAsFactors = FALSE
    )
  }
  indiv_df <- do.call(rbind, indiv_preds)

  return(list(population = pred_df, individual = indiv_df))
}

#' Prediction functions for nlme formulas (must be at module level)

.pred_1comp_ev <- function(lV, lk, lka, Time, Dose) {
  v <- exp(lV); k <- exp(lk); ka <- exp(lka)
  dka <- ka - k
  dka <- ifelse(abs(dka) < 1e-10, 1e-10, dka)
  (Dose * ka / (v * dka)) * (exp(-k * Time) - exp(-ka * Time))
}

.pred_2comp_iv <- function(lCL, lV1, lQ, lV2, Time, Dose) {
  cl <- exp(lCL); v1 <- exp(lV1); q <- exp(lQ); v2 <- exp(lV2)
  k10 <- cl / v1; k12 <- q / v1; k21 <- q / v2
  ss <- k10 + k12 + k21
  disc <- ss^2 - 4 * k10 * k21
  disc <- ifelse(disc < 0, 1e-20, disc)
  dd <- sqrt(disc)
  alpha <- 0.5 * (ss + dd); beta <- 0.5 * (ss - dd)
  dab <- alpha - beta
  dab <- ifelse(abs(dab) < 1e-20, 1e-20, dab)
  A <- (Dose / v1) * (alpha - k21) / dab
  B <- (Dose / v1) * (k21 - beta) / dab
  A * exp(-alpha * Time) + B * exp(-beta * Time)
}

.pred_2comp_ev <- function(lCL, lV1, lQ, lV2, lka, Time, Dose) {
  cl <- exp(lCL); v1 <- exp(lV1); q <- exp(lQ); v2 <- exp(lV2); ka <- exp(lka)
  k10 <- cl / v1; k12 <- q / v1; k21 <- q / v2
  ss <- k10 + k12 + k21
  disc <- ss^2 - 4 * k10 * k21
  disc <- ifelse(disc < 0, 1e-20, disc)
  dd <- sqrt(disc)
  alpha <- 0.5 * (ss + dd); beta <- 0.5 * (ss - dd)
  d_ka_a <- ka - alpha; d_ka_a <- ifelse(abs(d_ka_a) < 1e-20, 1e-20, d_ka_a)
  d_ka_b <- ka - beta;  d_ka_b <- ifelse(abs(d_ka_b) < 1e-20, 1e-20, d_ka_b)
  d_ba   <- beta - alpha; d_ba <- ifelse(abs(d_ba) < 1e-20, 1e-20, d_ba)
  d_ab   <- -d_ba
  d_ak   <- alpha - ka; d_ak <- ifelse(abs(d_ak) < 1e-20, 1e-20, d_ak)
  d_bk   <- beta - ka;  d_bk <- ifelse(abs(d_bk) < 1e-20, 1e-20, d_bk)
  coeff <- Dose * ka / v1
  A_c <- coeff * (k21 - alpha) / (d_ka_a * d_ba)
  B_c <- coeff * (k21 - beta)  / (d_ka_b * d_ab)
  C_c <- coeff * (k21 - ka)    / (d_ak * d_bk)
  A_c * exp(-alpha * Time) + B_c * exp(-beta * Time) + C_c * exp(-ka * Time)
}

.pred_3comp_iv <- function(lCL, lV1, lQ2, lV2, lQ3, lV3, Time, Dose) {
  cl <- exp(lCL); v1 <- exp(lV1); q2 <- exp(lQ2); v2 <- exp(lV2)
  q3 <- exp(lQ3); v3 <- exp(lV3)
  k10 <- cl / v1; k12 <- q2 / v1; k21 <- q2 / v2
  k13 <- q3 / v1; k31 <- q3 / v3
  a0 <- k10 * k21 * k31
  a1 <- k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
  a2 <- k10 + k12 + k13 + k21 + k31
  pp <- a1 - (a2^2) / 3
  qq <- (2 * a2^3 / 27) - (a1 * a2 / 3) + a0
  r1 <- (-(pp^3) / 27)^0.5
  phi <- acos(pmin(pmax(-qq / (2 * r1 + 1e-20), -1), 1)) / 3
  r1c <- r1^(1/3)
  rr1 <- -(cos(phi) * 2 * r1c - a2 / 3)
  rr2 <- -(cos(phi + 2 * pi / 3) * 2 * r1c - a2 / 3)
  rr3 <- -(cos(phi + 4 * pi / 3) * 2 * r1c - a2 / 3)
  alpha <- pmax(rr1, pmax(rr2, rr3))
  gamma <- pmin(rr1, pmin(rr2, rr3))
  beta <- rr1 + rr2 + rr3 - alpha - gamma
  d_ba <- beta - alpha; d_ba <- ifelse(abs(d_ba) < 1e-20, 1e-20, d_ba)
  d_ga <- gamma - alpha; d_ga <- ifelse(abs(d_ga) < 1e-20, 1e-20, d_ga)
  d_ab <- -d_ba
  d_gb <- gamma - beta; d_gb <- ifelse(abs(d_gb) < 1e-20, 1e-20, d_gb)
  d_ag <- -d_ga; d_bg <- -d_gb
  A <- (Dose / v1) * ((k21 - alpha) * (k31 - alpha)) / (d_ba * d_ga)
  B <- (Dose / v1) * ((k21 - beta) * (k31 - beta)) / (d_ab * d_gb)
  Cc <- (Dose / v1) * ((k21 - gamma) * (k31 - gamma)) / (d_ag * d_bg)
  A * exp(-alpha * Time) + B * exp(-beta * Time) + Cc * exp(-gamma * Time)
}

.pred_3comp_ev <- function(lCL, lV1, lQ2, lV2, lQ3, lV3, lka, Time, Dose) {
  cl <- exp(lCL); v1 <- exp(lV1); q2 <- exp(lQ2); v2 <- exp(lV2)
  q3 <- exp(lQ3); v3 <- exp(lV3); ka <- exp(lka)
  k10 <- cl / v1; k12 <- q2 / v1; k21 <- q2 / v2
  k13 <- q3 / v1; k31 <- q3 / v3
  a0 <- k10 * k21 * k31
  a1 <- k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
  a2 <- k10 + k12 + k13 + k21 + k31
  pp <- a1 - (a2^2) / 3
  qq <- (2 * a2^3 / 27) - (a1 * a2 / 3) + a0
  r1 <- (-(pp^3) / 27)^0.5
  phi <- acos(pmin(pmax(-qq / (2 * r1 + 1e-20), -1), 1)) / 3
  r1c <- r1^(1/3)
  rr1 <- -(cos(phi) * 2 * r1c - a2 / 3)
  rr2 <- -(cos(phi + 2 * pi / 3) * 2 * r1c - a2 / 3)
  rr3 <- -(cos(phi + 4 * pi / 3) * 2 * r1c - a2 / 3)
  alpha <- pmax(rr1, pmax(rr2, rr3))
  gamma <- pmin(rr1, pmin(rr2, rr3))
  beta <- rr1 + rr2 + rr3 - alpha - gamma
  safe_div <- function(x) ifelse(abs(x) < 1e-20, 1e-20, x)
  coeff <- Dose * ka / v1
  A_c <- coeff * ((k21 - alpha) * (k31 - alpha)) / (safe_div(beta - alpha) * safe_div(gamma - alpha) * safe_div(ka - alpha))
  B_c <- coeff * ((k21 - beta) * (k31 - beta)) / (safe_div(alpha - beta) * safe_div(gamma - beta) * safe_div(ka - beta))
  C_c <- coeff * ((k21 - gamma) * (k31 - gamma)) / (safe_div(alpha - gamma) * safe_div(beta - gamma) * safe_div(ka - gamma))
  D_c <- coeff * ((k21 - ka) * (k31 - ka)) / (safe_div(alpha - ka) * safe_div(beta - ka) * safe_div(gamma - ka))
  A_c * exp(-alpha * Time) + B_c * exp(-beta * Time) + C_c * exp(-gamma * Time) + D_c * exp(-ka * Time)
}


# ============================================================
# Parameter Extraction
# ============================================================

#' Extract PK parameters from fitted nlme model
#'
#' @param fit nlme fit object
#' @param model_type Character model type
#' @param route Character route
#' @param dose Numeric dose
#' @return List with 'population' and 'individual' data frames
extract_compartmental_params <- function(fit, model_type, route, dose) {
  fe <- fixef(fit)
  re <- ranef(fit)
  is_iv <- toupper(route) == "IV"

  # Extract vcov and residual df for SE/CI computation
  vcov_mat <- tryCatch(vcov(fit), error = function(e) NULL)
  df_resid <- tryCatch({
    tt <- summary(fit)$tTable
    if (!is.null(tt) && "DF" %in% colnames(tt)) tt[1, "DF"] else nrow(fit$data) - length(fe)
  }, error = function(e) nrow(fit$data) - length(fe))

  # Select extraction function based on model type
  extract_fn <- switch(model_type,
    "1comp" = extract_1comp_params,
    "2comp" = extract_2comp_params,
    "3comp" = extract_3comp_params
  )

  # Population parameters (with SE/CI)
  pop_params <- extract_fn(fe, dose, is_iv, vcov_mat, df_resid)

  # Individual parameters in wide format (one row per subject, parameter names as columns)
  indiv_list <- list()
  for (i in 1:nrow(re)) {
    subj_id <- rownames(re)[i]
    indiv_fe <- fe
    for (nm in names(re)) {
      indiv_fe[nm] <- indiv_fe[nm] + re[i, nm]
    }
    ip <- extract_fn(indiv_fe, dose, is_iv)
    row_df <- as.data.frame(as.list(setNames(ip$Estimate, ip$Parameter)),
                            check.names = FALSE, stringsAsFactors = FALSE)
    row_df$ID <- subj_id
    indiv_list[[i]] <- row_df
  }
  indiv_df <- do.call(rbind, indiv_list)

  # Add model diagnostics
  pop_params$AIC <- AIC(fit)
  pop_params$BIC <- BIC(fit)
  pop_params$logLik <- as.numeric(logLik(fit))

  return(list(population = pop_params, individual = indiv_df))
}

#' Compute SE for a derived parameter using the numerical delta method
#' @param fe Named vector of fixed effects (log-scale)
#' @param vcov_mat Variance-covariance matrix of fixed effects
#' @param fun Function(fe) -> scalar derived parameter value
#' @return Standard error on the natural scale (NA if computation fails)
delta_method_se <- function(fe, vcov_mat, fun) {
  tryCatch({
    eps <- 1e-5
    grad <- numeric(length(fe))
    f0 <- fun(fe)
    for (j in seq_along(fe)) {
      fe_plus <- fe
      fe_plus[j] <- fe_plus[j] + eps
      grad[j] <- (fun(fe_plus) - f0) / eps
    }
    var_est <- as.numeric(t(grad) %*% vcov_mat %*% grad)
    if (var_est >= 0) sqrt(var_est) else NA_real_
  }, error = function(e) NA_real_)
}

#' Add SE and 95% CI columns to a parameter data frame using the delta method
#' @param params Data frame with Parameter and Estimate columns
#' @param fe Named vector of fixed effects (log-scale)
#' @param vcov_mat Variance-covariance matrix of fixed effects
#' @param param_funs Named list of functions mapping fe -> derived value (names match params$Parameter)
#' @param df_resid Residual degrees of freedom for t-distribution
#' @return params with SE, CI_lower, CI_upper columns appended
add_se_ci_to_params <- function(params, fe, vcov_mat, param_funs, df_resid) {
  se_vec <- numeric(nrow(params))
  for (i in seq_len(nrow(params))) {
    pname <- params$Parameter[i]
    if (pname %in% names(param_funs)) {
      se_vec[i] <- delta_method_se(fe, vcov_mat, param_funs[[pname]])
    } else {
      se_vec[i] <- NA_real_
    }
  }
  t_crit <- tryCatch(qt(0.975, df = df_resid), error = function(e) 1.96)
  params$SE <- se_vec
  params$CI_lower <- params$Estimate - t_crit * se_vec
  params$CI_upper <- params$Estimate + t_crit * se_vec
  return(params)
}

extract_1comp_params <- function(fe, dose, is_iv, vcov_mat = NULL, df_resid = NULL) {
  V <- exp(fe["lV"])
  k <- exp(fe["lk"])
  CL <- V * k
  t_half <- log(2) / k
  Vd <- V

  params <- data.frame(
    Parameter = c(if (is_iv) "V (L)" else "V/F (L)",
                  if (is_iv) "CL (L/h)" else "CL/F (L/h)",
                  "k (1/h)", "Half-life (h)"),
    Estimate = c(V, CL, k, t_half),
    stringsAsFactors = FALSE
  )

  if (!is_iv && "lka" %in% names(fe)) {
    ka <- exp(fe["lka"])
    params <- rbind(params, data.frame(
      Parameter = "ka (1/h)", Estimate = ka, stringsAsFactors = FALSE
    ))
  }

  rownames(params) <- NULL

  # Add SE and 95% CI if vcov available
  if (!is.null(vcov_mat) && !is.null(df_resid)) {
    v_name <- if (is_iv) "V (L)" else "V/F (L)"
    cl_name <- if (is_iv) "CL (L/h)" else "CL/F (L/h)"
    pfuns <- list()
    pfuns[[v_name]]         <- function(x) exp(x["lV"])
    pfuns[[cl_name]]        <- function(x) exp(x["lV"]) * exp(x["lk"])
    pfuns[["k (1/h)"]]      <- function(x) exp(x["lk"])
    pfuns[["Half-life (h)"]] <- function(x) log(2) / exp(x["lk"])
    if (!is_iv && "lka" %in% names(fe)) {
      pfuns[["ka (1/h)"]] <- function(x) exp(x["lka"])
    }
    params <- add_se_ci_to_params(params, fe, vcov_mat, pfuns, df_resid)
  }

  return(params)
}

extract_2comp_params <- function(fe, dose, is_iv, vcov_mat = NULL, df_resid = NULL) {
  CL <- exp(fe["lCL"])
  V1 <- exp(fe["lV1"])
  Q  <- exp(fe["lQ"])
  V2 <- exp(fe["lV2"])

  k10 <- CL / V1
  k12 <- Q / V1
  k21 <- Q / V2
  Vss <- V1 + V2

  ss <- k10 + k12 + k21
  dd <- sqrt(ss^2 - 4 * k10 * k21)
  alpha <- 0.5 * (ss + dd)
  beta  <- 0.5 * (ss - dd)
  t_half_alpha <- log(2) / alpha
  t_half_beta  <- log(2) / beta

  pref <- if (is_iv) "" else "/F"
  params <- data.frame(
    Parameter = c(paste0("CL", pref, " (L/h)"),
                  paste0("V1", pref, " (L)"),
                  paste0("Q", pref, " (L/h)"),
                  paste0("V2", pref, " (L)"),
                  paste0("Vss", pref, " (L)"),
                  "Alpha (1/h)", "Beta (1/h)",
                  "t1/2 alpha (h)", "t1/2 beta (h)"),
    Estimate = c(CL, V1, Q, V2, Vss, alpha, beta, t_half_alpha, t_half_beta),
    stringsAsFactors = FALSE
  )

  if (!is_iv && "lka" %in% names(fe)) {
    ka <- exp(fe["lka"])
    params <- rbind(params, data.frame(
      Parameter = "ka (1/h)", Estimate = ka, stringsAsFactors = FALSE
    ))
  }

  rownames(params) <- NULL

  # Add SE and 95% CI if vcov available
  if (!is.null(vcov_mat) && !is.null(df_resid)) {
    # Helper: compute alpha/beta from log-scale params
    .alpha_beta <- function(x) {
      cl <- exp(x["lCL"]); v1 <- exp(x["lV1"]); q <- exp(x["lQ"]); v2 <- exp(x["lV2"])
      k10 <- cl/v1; k12 <- q/v1; k21 <- q/v2
      s <- k10 + k12 + k21
      d <- sqrt(s^2 - 4 * k10 * k21)
      c(alpha = 0.5*(s+d), beta = 0.5*(s-d))
    }

    pfuns <- list()
    pfuns[[paste0("CL", pref, " (L/h)")]]  <- function(x) exp(x["lCL"])
    pfuns[[paste0("V1", pref, " (L)")]]     <- function(x) exp(x["lV1"])
    pfuns[[paste0("Q", pref, " (L/h)")]]    <- function(x) exp(x["lQ"])
    pfuns[[paste0("V2", pref, " (L)")]]     <- function(x) exp(x["lV2"])
    pfuns[[paste0("Vss", pref, " (L)")]]    <- function(x) exp(x["lV1"]) + exp(x["lV2"])
    pfuns[["Alpha (1/h)"]]       <- function(x) .alpha_beta(x)["alpha"]
    pfuns[["Beta (1/h)"]]        <- function(x) .alpha_beta(x)["beta"]
    pfuns[["t1/2 alpha (h)"]]    <- function(x) log(2) / .alpha_beta(x)["alpha"]
    pfuns[["t1/2 beta (h)"]]     <- function(x) log(2) / .alpha_beta(x)["beta"]
    if (!is_iv && "lka" %in% names(fe)) {
      pfuns[["ka (1/h)"]] <- function(x) exp(x["lka"])
    }
    params <- add_se_ci_to_params(params, fe, vcov_mat, pfuns, df_resid)
  }

  return(params)
}

extract_3comp_params <- function(fe, dose, is_iv, vcov_mat = NULL, df_resid = NULL) {
  CL <- exp(fe["lCL"])
  V1 <- exp(fe["lV1"])
  Q2 <- exp(fe["lQ2"])
  V2 <- exp(fe["lV2"])
  Q3 <- exp(fe["lQ3"])
  V3 <- exp(fe["lV3"])

  Vss <- V1 + V2 + V3
  k10 <- CL / V1

  pref <- if (is_iv) "" else "/F"
  params <- data.frame(
    Parameter = c(paste0("CL", pref, " (L/h)"),
                  paste0("V1", pref, " (L)"),
                  paste0("Q2", pref, " (L/h)"),
                  paste0("V2", pref, " (L)"),
                  paste0("Q3", pref, " (L/h)"),
                  paste0("V3", pref, " (L)"),
                  paste0("Vss", pref, " (L)")),
    Estimate = c(CL, V1, Q2, V2, Q3, V3, Vss),
    stringsAsFactors = FALSE
  )

  if (!is_iv && "lka" %in% names(fe)) {
    ka <- exp(fe["lka"])
    params <- rbind(params, data.frame(
      Parameter = "ka (1/h)", Estimate = ka, stringsAsFactors = FALSE
    ))
  }

  rownames(params) <- NULL

  # Add SE and 95% CI if vcov available
  if (!is.null(vcov_mat) && !is.null(df_resid)) {
    pfuns <- list()
    pfuns[[paste0("CL", pref, " (L/h)")]]  <- function(x) exp(x["lCL"])
    pfuns[[paste0("V1", pref, " (L)")]]     <- function(x) exp(x["lV1"])
    pfuns[[paste0("Q2", pref, " (L/h)")]]   <- function(x) exp(x["lQ2"])
    pfuns[[paste0("V2", pref, " (L)")]]     <- function(x) exp(x["lV2"])
    pfuns[[paste0("Q3", pref, " (L/h)")]]   <- function(x) exp(x["lQ3"])
    pfuns[[paste0("V3", pref, " (L)")]]     <- function(x) exp(x["lV3"])
    pfuns[[paste0("Vss", pref, " (L)")]]    <- function(x) exp(x["lV1"]) + exp(x["lV2"]) + exp(x["lV3"])
    if (!is_iv && "lka" %in% names(fe)) {
      pfuns[["ka (1/h)"]] <- function(x) exp(x["lka"])
    }
    params <- add_se_ci_to_params(params, fe, vcov_mat, pfuns, df_resid)
  }

  return(params)
}


# ============================================================
# Predictions
# ============================================================

#' Generate model predictions for plotting
#'
#' @param fit nlme fit object
#' @param data Original data
#' @param model_type Character model type
#' @param route Character route
#' @param dose Numeric dose
#' @return Data frame with Time, Pred_pop, ID, Pred_indiv columns
generate_predictions <- function(fit, data, model_type, route, dose) {
  time_range <- range(data$Time)
  time_seq <- seq(time_range[1], time_range[2], length.out = 200)

  # Population predictions
  new_data <- data.frame(
    Time = time_seq,
    Dose = dose,
    ID = data$ID[1]  # placeholder for population prediction
  )

  # Population-level prediction using fixed effects only
  fe <- fixef(fit)
  is_iv <- toupper(route) == "IV"

  pop_pred <- predict_from_params(fe, time_seq, dose, model_type, is_iv)

  pred_df <- data.frame(Time = time_seq, Pred_pop = pop_pred)

  # Individual predictions
  re <- ranef(fit)
  indiv_preds <- list()
  for (i in 1:nrow(re)) {
    subj_id <- rownames(re)[i]
    indiv_fe <- fe
    for (nm in names(re)) {
      indiv_fe[nm] <- indiv_fe[nm] + re[i, nm]
    }
    indiv_pred <- predict_from_params(indiv_fe, time_seq, dose, model_type, is_iv)
    indiv_preds[[subj_id]] <- data.frame(
      Time = time_seq,
      Pred_indiv = indiv_pred,
      ID = subj_id,
      stringsAsFactors = FALSE
    )
  }

  indiv_df <- do.call(rbind, indiv_preds)

  return(list(population = pred_df, individual = indiv_df))
}

#' Compute predictions from parameter vector
predict_from_params <- function(params, time, dose, model_type, is_iv) {
  if (model_type == "1comp") {
    V <- exp(params["lV"])
    k <- exp(params["lk"])
    if (is_iv) {
      pred <- (dose / V) * exp(-k * time)
    } else {
      ka <- exp(params["lka"])
      if (abs(ka - k) < 1e-10) {
        pred <- (dose / V) * k * time * exp(-k * time)
      } else {
        pred <- (dose * ka / (V * (ka - k))) * (exp(-k * time) - exp(-ka * time))
      }
    }
  } else if (model_type == "2comp") {
    CL <- exp(params["lCL"]); V1 <- exp(params["lV1"])
    Q <- exp(params["lQ"]); V2 <- exp(params["lV2"])
    k10 <- CL / V1; k12 <- Q / V1; k21 <- Q / V2
    ss <- k10 + k12 + k21
    dd <- sqrt(ss^2 - 4 * k10 * k21)
    alpha <- 0.5 * (ss + dd); beta <- 0.5 * (ss - dd)

    if (is_iv) {
      A <- (dose / V1) * (alpha - k21) / (alpha - beta)
      B <- (dose / V1) * (k21 - beta) / (alpha - beta)
      pred <- A * exp(-alpha * time) + B * exp(-beta * time)
    } else {
      ka <- exp(params["lka"])
      A_c <- (dose * ka / V1) * (k21 - alpha) / ((ka - alpha) * (beta - alpha))
      B_c <- (dose * ka / V1) * (k21 - beta) / ((ka - beta) * (alpha - beta))
      C_c <- (dose * ka / V1) * (k21 - ka) / ((alpha - ka) * (beta - ka))
      pred <- A_c * exp(-alpha * time) + B_c * exp(-beta * time) + C_c * exp(-ka * time)
    }
  } else if (model_type == "3comp") {
    CL <- exp(params["lCL"]); V1 <- exp(params["lV1"])
    Q2 <- exp(params["lQ2"]); V2 <- exp(params["lV2"])
    Q3 <- exp(params["lQ3"]); V3 <- exp(params["lV3"])
    k10 <- CL / V1; k12 <- Q2 / V1; k21 <- Q2 / V2
    k13 <- Q3 / V1; k31 <- Q3 / V3
    a0 <- k10 * k21 * k31
    a1 <- k10 * k31 + k21 * k31 + k21 * k13 + k10 * k21 + k31 * k12
    a2 <- k10 + k12 + k13 + k21 + k31
    pp <- a1 - (a2^2) / 3
    qq <- (2 * a2^3 / 27) - (a1 * a2 / 3) + a0
    r1 <- (-(pp^3) / 27)^0.5
    phi <- acos(-qq / (2 * r1)) / 3
    r1c <- r1^(1/3)
    rr <- c(
      -(cos(phi) * 2 * r1c - a2 / 3),
      -(cos(phi + 2 * pi / 3) * 2 * r1c - a2 / 3),
      -(cos(phi + 4 * pi / 3) * 2 * r1c - a2 / 3)
    )
    rr <- sort(rr, decreasing = TRUE)
    alpha <- rr[1]; beta <- rr[2]; gamma <- rr[3]

    if (is_iv) {
      A <- (dose / V1) * ((k21 - alpha) * (k31 - alpha)) / ((beta - alpha) * (gamma - alpha))
      B <- (dose / V1) * ((k21 - beta) * (k31 - beta)) / ((alpha - beta) * (gamma - beta))
      C <- (dose / V1) * ((k21 - gamma) * (k31 - gamma)) / ((alpha - gamma) * (beta - gamma))
      pred <- A * exp(-alpha * time) + B * exp(-beta * time) + C * exp(-gamma * time)
    } else {
      ka <- exp(params["lka"])
      A_c <- (dose * ka / V1) * ((k21 - alpha) * (k31 - alpha)) / ((beta - alpha) * (gamma - alpha) * (ka - alpha))
      B_c <- (dose * ka / V1) * ((k21 - beta) * (k31 - beta)) / ((alpha - beta) * (gamma - beta) * (ka - beta))
      C_c <- (dose * ka / V1) * ((k21 - gamma) * (k31 - gamma)) / ((alpha - gamma) * (beta - gamma) * (ka - gamma))
      D_c <- (dose * ka / V1) * ((k21 - ka) * (k31 - ka)) / ((alpha - ka) * (beta - ka) * (gamma - ka))
      pred <- A_c * exp(-alpha * time) + B_c * exp(-beta * time) + C_c * exp(-gamma * time) + D_c * exp(-ka * time)
    }
  }

  # Ensure non-negative
  pred[pred < 0] <- 0
  return(pred)
}
