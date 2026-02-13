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
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param dose Numeric dose
#' @param model_type Character: "1comp", "2comp", "3comp"
#' @param route Character: "IV" or "EV"
#' @return Named list of starting parameter values (on log scale)
estimate_initial_params <- function(df, dose, model_type, route) {
  # Aggregate data for initial estimates
  mean_data <- aggregate(Conc ~ Time, data = df, FUN = mean)
  mean_data <- mean_data[order(mean_data$Time), ]

  time <- mean_data$Time
  conc <- mean_data$Conc

  cmax <- max(conc)
  tmax_idx <- which.max(conc)

  # Rough V estimate
  if (toupper(route) == "IV" && conc[1] > 0) {
    V_est <- dose / conc[1]
  } else {
    V_est <- dose / cmax
  }
  V_est <- max(V_est, 0.01)

  # Rough terminal slope
  post_peak <- conc[tmax_idx:length(conc)]
  post_time <- time[tmax_idx:length(time)]
  pos_idx <- post_peak > 0
  if (sum(pos_idx) >= 2) {
    log_conc <- log(post_peak[pos_idx])
    t_sub <- post_time[pos_idx]
    fit <- lm(log_conc ~ t_sub)
    k_est <- max(-coef(fit)[2], 0.01)
  } else {
    k_est <- 0.1
  }

  CL_est <- k_est * V_est

  if (model_type == "1comp") {
    params <- list(lV = log(V_est), lk = log(k_est))
    if (toupper(route) != "IV") {
      # Estimate ka from ascending phase
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
    if (toupper(route) != "IV") {
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
    if (toupper(route) != "IV") {
      ka_est <- max(k_est * 3, 0.5)
      params$lka <- log(ka_est)
    }
  }

  return(params)
}


# ============================================================
# NLME Fitting
# ============================================================

#' Fit a compartmental model using nlme
#'
#' @param df Data frame with ID, Time, Conc columns
#' @param dose Numeric dose
#' @param model_type Character: "1comp", "2comp", "3comp"
#' @param route Character: "IV" or "EV"
#' @return List with 'fit' (nlme object or NULL), 'summary', 'params', 'error'
fit_compartmental_model <- function(df, dose, model_type, route) {
  route_upper <- toupper(route)
  is_iv <- route_upper == "IV"

  # Prepare data as groupedData
  pk_data <- df[, c("ID", "Time", "Conc")]
  pk_data$Dose <- dose

  # Get initial estimates
  init_params <- estimate_initial_params(df, dose, model_type, route)

  # Build the model - try full random effects first, fall back to reduced
  fit <- NULL
  fit_error <- NULL

  # Attempt 1: Full random effects (lCL + lV1 ~ 1 | ID)
  fit <- tryCatch({
    if (model_type == "1comp" && is_iv) {
      fit_1comp_iv(pk_data, init_params)
    } else if (model_type == "1comp" && !is_iv) {
      fit_1comp_ev(pk_data, init_params)
    } else if (model_type == "2comp" && is_iv) {
      fit_2comp_iv(pk_data, init_params)
    } else if (model_type == "2comp" && !is_iv) {
      fit_2comp_ev(pk_data, init_params)
    } else if (model_type == "3comp" && is_iv) {
      fit_3comp_iv(pk_data, init_params)
    } else if (model_type == "3comp" && !is_iv) {
      fit_3comp_ev(pk_data, init_params)
    }
  }, error = function(e) {
    fit_error <<- e$message
    NULL
  })

  # Attempt 2: Reduced random effects (lCL ~ 1 | ID) for 2comp/3comp
  if (is.null(fit) && model_type %in% c("2comp", "3comp")) {
    message("Full model failed, trying reduced random effects...")
    fit <- tryCatch({
      if (model_type == "2comp" && is_iv) {
        fit_2comp_iv_reduced(pk_data, init_params)
      } else if (model_type == "2comp" && !is_iv) {
        fit_2comp_ev_reduced(pk_data, init_params)
      } else if (model_type == "3comp" && is_iv) {
        fit_3comp_iv_reduced(pk_data, init_params)
      } else if (model_type == "3comp" && !is_iv) {
        fit_3comp_ev_reduced(pk_data, init_params)
      }
    }, error = function(e) {
      fit_error <<- e$message
      NULL
    })
  }

  if (is.null(fit)) {
    return(list(
      fit = NULL, params = NULL, summary = NULL, predictions = NULL,
      error = paste("Model fitting failed:", fit_error,
                    "\nTry a simpler model or check your data quality.")
    ))
  }

  # Extract parameters
  tryCatch({
    params <- extract_compartmental_params(fit, model_type, route, dose)
    return(list(
      fit = fit,
      params = params$individual,
      summary = params$population,
      predictions = generate_predictions(fit, pk_data, model_type, route, dose),
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

#' Fit one-compartment IV model
fit_1comp_iv <- function(data, inits) {
  nlme(
    Conc ~ (Dose / exp(lV)) * exp(-exp(lk) * Time),
    data = data,
    fixed = lV + lk ~ 1,
    random = lV + lk ~ 1 | ID,
    start = c(lV = inits$lV, lk = inits$lk),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' One-compartment extravascular prediction function
.pred_1comp_ev <- function(lV, lk, lka, Time, Dose) {
  v <- exp(lV); k <- exp(lk); ka <- exp(lka)
  dka <- ka - k
  dka <- ifelse(abs(dka) < 1e-10, 1e-10, dka)
  (Dose * ka / (v * dka)) * (exp(-k * Time) - exp(-ka * Time))
}

#' Fit one-compartment extravascular model
fit_1comp_ev <- function(data, inits) {
  nlme(
    Conc ~ .pred_1comp_ev(lV, lk, lka, Time, Dose),
    data = data,
    fixed = lV + lk + lka ~ 1,
    random = lV + lk ~ 1 | ID,
    start = c(lV = inits$lV, lk = inits$lk, lka = inits$lka),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Two-compartment IV prediction function
#' Defined at module level for reliable nlme formula evaluation
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

#' Fit two-compartment IV model
fit_2comp_iv <- function(data, inits) {
  nlme(
    Conc ~ .pred_2comp_iv(lCL, lV1, lQ, lV2, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ + lV2 ~ 1,
    random = lCL + lV1 ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ, lV2 = inits$lV2),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Fit two-compartment IV model with reduced random effects (fallback)
fit_2comp_iv_reduced <- function(data, inits) {
  nlme(
    Conc ~ .pred_2comp_iv(lCL, lV1, lQ, lV2, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ + lV2 ~ 1,
    random = lCL ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ, lV2 = inits$lV2),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.5, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Two-compartment extravascular prediction function
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

#' Fit two-compartment extravascular model
fit_2comp_ev <- function(data, inits) {
  nlme(
    Conc ~ .pred_2comp_ev(lCL, lV1, lQ, lV2, lka, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ + lV2 + lka ~ 1,
    random = lCL + lV1 ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ,
              lV2 = inits$lV2, lka = inits$lka),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Fit two-compartment extravascular model with reduced random effects (fallback)
fit_2comp_ev_reduced <- function(data, inits) {
  nlme(
    Conc ~ .pred_2comp_ev(lCL, lV1, lQ, lV2, lka, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ + lV2 + lka ~ 1,
    random = lCL ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ = inits$lQ,
              lV2 = inits$lV2, lka = inits$lka),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.5, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Three-compartment IV prediction function
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

#' Fit three-compartment IV model
fit_3comp_iv <- function(data, inits) {
  nlme(
    Conc ~ .pred_3comp_iv(lCL, lV1, lQ2, lV2, lQ3, lV3, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 ~ 1,
    random = lCL + lV1 ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
              lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Fit three-compartment IV model with reduced random effects (fallback)
fit_3comp_iv_reduced <- function(data, inits) {
  nlme(
    Conc ~ .pred_3comp_iv(lCL, lV1, lQ2, lV2, lQ3, lV3, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 ~ 1,
    random = lCL ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
              lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.5, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Three-compartment extravascular prediction function
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

#' Fit three-compartment extravascular model
fit_3comp_ev <- function(data, inits) {
  nlme(
    Conc ~ .pred_3comp_ev(lCL, lV1, lQ2, lV2, lQ3, lV3, lka, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 + lka ~ 1,
    random = lCL + lV1 ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
              lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3, lka = inits$lka),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.1, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
}

#' Fit three-compartment extravascular model with reduced random effects (fallback)
fit_3comp_ev_reduced <- function(data, inits) {
  nlme(
    Conc ~ .pred_3comp_ev(lCL, lV1, lQ2, lV2, lQ3, lV3, lka, Time, Dose),
    data = data,
    fixed = lCL + lV1 + lQ2 + lV2 + lQ3 + lV3 + lka ~ 1,
    random = lCL ~ 1 | ID,
    start = c(lCL = inits$lCL, lV1 = inits$lV1, lQ2 = inits$lQ2,
              lV2 = inits$lV2, lQ3 = inits$lQ3, lV3 = inits$lV3, lka = inits$lka),
    control = nlmeControl(maxIter = 200, pnlsTol = 0.5, msMaxIter = 200,
                          returnObject = TRUE),
    na.action = na.omit
  )
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

  if (model_type == "1comp") {
    pop_params <- extract_1comp_params(fe, dose, is_iv)
    # Individual parameters
    indiv_list <- list()
    for (i in 1:nrow(re)) {
      subj_id <- rownames(re)[i]
      indiv_fe <- fe
      for (nm in names(re)) {
        indiv_fe[nm] <- indiv_fe[nm] + re[i, nm]
      }
      ip <- extract_1comp_params(indiv_fe, dose, is_iv)
      ip$ID <- subj_id
      indiv_list[[i]] <- ip
    }
    indiv_df <- do.call(rbind, lapply(indiv_list, as.data.frame, stringsAsFactors = FALSE))
  } else if (model_type == "2comp") {
    pop_params <- extract_2comp_params(fe, dose, is_iv)
    indiv_list <- list()
    for (i in 1:nrow(re)) {
      subj_id <- rownames(re)[i]
      indiv_fe <- fe
      for (nm in names(re)) {
        indiv_fe[nm] <- indiv_fe[nm] + re[i, nm]
      }
      ip <- extract_2comp_params(indiv_fe, dose, is_iv)
      ip$ID <- subj_id
      indiv_list[[i]] <- ip
    }
    indiv_df <- do.call(rbind, lapply(indiv_list, as.data.frame, stringsAsFactors = FALSE))
  } else if (model_type == "3comp") {
    pop_params <- extract_3comp_params(fe, dose, is_iv)
    indiv_list <- list()
    for (i in 1:nrow(re)) {
      subj_id <- rownames(re)[i]
      indiv_fe <- fe
      for (nm in names(re)) {
        indiv_fe[nm] <- indiv_fe[nm] + re[i, nm]
      }
      ip <- extract_3comp_params(indiv_fe, dose, is_iv)
      ip$ID <- subj_id
      indiv_list[[i]] <- ip
    }
    indiv_df <- do.call(rbind, lapply(indiv_list, as.data.frame, stringsAsFactors = FALSE))
  }

  # Add model diagnostics
  pop_params$AIC <- AIC(fit)
  pop_params$BIC <- BIC(fit)
  pop_params$logLik <- as.numeric(logLik(fit))

  return(list(population = pop_params, individual = indiv_df))
}

extract_1comp_params <- function(fe, dose, is_iv) {
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
  return(params)
}

extract_2comp_params <- function(fe, dose, is_iv) {
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
  return(params)
}

extract_3comp_params <- function(fe, dose, is_iv) {
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
