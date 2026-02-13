# Data Validation and Preprocessing Module

library(readxl)
library(tools)

#' Validate and load PK data from uploaded file
#'
#' Expects wide format: first column is Time, each subsequent column
#' contains concentrations for one animal/subject. Column headers become
#' subject IDs (e.g., "Animal1", "Dog_01", "1").
#' Optionally a "Drug" or "Compound" column can be present for multi-compound data.
#'
#' @param file_path Path to uploaded file
#' @param file_name Original file name (for extension detection)
#' @return List with 'data' (data.frame in long format or NULL) and 'error' (character or NULL)
load_pk_data <- function(file_path, file_name) {
  ext <- tolower(file_ext(file_name))

  tryCatch({
    if (ext == "csv") {
      df <- read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE)
    } else if (ext %in% c("xls", "xlsx")) {
      df <- as.data.frame(readxl::read_excel(file_path), check.names = FALSE)
    } else {
      return(list(data = NULL, error = "Unsupported file format. Please upload a CSV or XLS/XLSX file."))
    }

    if (ncol(df) < 2) {
      return(list(data = NULL, error = "File must have at least 2 columns: Time and one or more animal concentration columns."))
    }

    validation <- validate_and_reshape(df)
    if (!is.null(validation$error)) {
      return(validation)
    }

    df_long <- standardize_pk_data(validation$data)
    return(list(data = df_long, error = NULL))

  }, error = function(e) {
    return(list(data = NULL, error = paste("Error reading file:", e$message)))
  })
}

#' Identify the Time column, optional Drug column, and reshape wide-format data to long format
#'
#' @param df Data frame in wide format (Time + optional Drug + animal columns)
#' @return List with 'data' (long-format data.frame) and 'error' (character or NULL)
validate_and_reshape <- function(df) {
  # Find the Time column (case-insensitive)
  time_aliases <- c("time", "time_h", "time_hr", "time_min", "timepoint", "time.h",
                     "time.hr", "time.min", "time (h)", "time (min)", "time (hr)")
  cols_lower <- tolower(names(df))

  time_col_idx <- NULL
  for (alias in time_aliases) {
    idx <- which(cols_lower == alias)
    if (length(idx) > 0) {
      time_col_idx <- idx[1]
      break
    }
  }

  # If no alias matched, assume the first column is Time
  if (is.null(time_col_idx)) {
    time_col_idx <- 1
    message("No 'Time' column header found; using first column as Time.")
  }

  # Detect Drug/Compound column (case-insensitive)
  drug_aliases <- c("drug", "compound", "treatment", "group")
  drug_col_idx <- NULL
  for (alias in drug_aliases) {
    idx <- which(cols_lower == alias)
    if (length(idx) > 0) {
      drug_col_idx <- idx[1]
      break
    }
  }

  time_col_name <- names(df)[time_col_idx]
  exclude_idx <- time_col_idx
  drug_col_name <- NULL
  if (!is.null(drug_col_idx)) {
    drug_col_name <- names(df)[drug_col_idx]
    exclude_idx <- c(exclude_idx, drug_col_idx)
  }

  other_col_names <- names(df)[-exclude_idx]

  if (length(other_col_names) < 1) {
    return(list(data = NULL, error = "No concentration columns found. File needs Time plus at least one animal column."))
  }

  # Check that Time column is numeric
  time_vals <- suppressWarnings(as.numeric(df[[time_col_name]]))
  if (all(is.na(time_vals))) {
    return(list(data = NULL,
                error = paste0("The Time column ('", time_col_name, "') does not contain numeric values.")))
  }

  # Keep only columns that contain numeric data (skip ID, text columns like Average, SD if non-numeric)
  conc_col_names <- c()
  skipped <- c()
  for (cn in other_col_names) {
    # Treat "0%" and similar percentage strings as numeric by cleaning first
    cleaned_vals <- gsub("%$", "", trimws(as.character(df[[cn]])))
    test_vals <- suppressWarnings(as.numeric(cleaned_vals))
    if (all(is.na(test_vals))) {
      skipped <- c(skipped, cn)
    } else {
      conc_col_names <- c(conc_col_names, cn)
    }
  }

  # Skip columns named Average/Mean/SD/SEM (summary columns)
  summary_aliases <- c("average", "avg", "mean", "sd", "sem", "se", "stdev", "cv", "cv%")
  keep_cols <- c()
  for (cn in conc_col_names) {
    if (tolower(trimws(cn)) %in% summary_aliases) {
      skipped <- c(skipped, cn)
    } else {
      keep_cols <- c(keep_cols, cn)
    }
  }
  conc_col_names <- keep_cols

  if (length(skipped) > 0) {
    message(paste("Skipped non-numeric or summary columns:", paste(skipped, collapse = ", ")))
  }

  if (length(conc_col_names) < 1) {
    return(list(data = NULL, error = "No numeric concentration columns found after the Time column."))
  }

  # Get drug values if drug column exists
  drug_vals <- NULL
  if (!is.null(drug_col_name)) {
    drug_vals <- trimws(as.character(df[[drug_col_name]]))
  }

  # Reshape from wide to long
  long_list <- list()
  for (cn in conc_col_names) {
    # Clean percentage signs and convert to numeric
    cleaned_vals <- gsub("%$", "", trimws(as.character(df[[cn]])))
    conc_numeric <- suppressWarnings(as.numeric(cleaned_vals))

    subj_df <- data.frame(
      ID = cn,
      Time = time_vals,
      Conc = conc_numeric,
      stringsAsFactors = FALSE
    )

    if (!is.null(drug_vals)) {
      subj_df$Drug <- drug_vals
    }

    long_list[[cn]] <- subj_df
  }

  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL

  return(list(data = long_df, error = NULL))
}

#' Standardize PK data types and handle missing values
#'
#' @param df Data frame with ID, Time, Conc (and optionally Drug) columns (long format)
#' @return Cleaned data frame
standardize_pk_data <- function(df) {
  df$ID <- as.factor(df$ID)
  df$Time <- as.numeric(df$Time)
  df$Conc <- as.numeric(df$Conc)

  if ("Drug" %in% names(df)) {
    df$Drug <- as.factor(df$Drug)
  }

  # Remove rows with NA in critical columns
  n_before <- nrow(df)
  df <- df[!is.na(df$Time) & !is.na(df$Conc), ]
  n_after <- nrow(df)

  if (n_before > n_after) {
    message(paste("Removed", n_before - n_after, "rows with missing Time or Conc values."))
  }

  # Sort by Drug (if present), then ID, then Time
  if ("Drug" %in% names(df)) {
    df <- df[order(df$Drug, df$ID, df$Time), ]
  } else {
    df <- df[order(df$ID, df$Time), ]
  }
  rownames(df) <- NULL

  return(df)
}

#' Get summary statistics for uploaded data
#'
#' @param df Validated PK data frame (long format)
#' @return Named list of summary info
pk_data_summary <- function(df) {
  has_drug <- "Drug" %in% names(df)
  summary_info <- list(
    n_subjects = length(unique(df$ID)),
    n_observations = nrow(df),
    time_range = range(df$Time, na.rm = TRUE),
    conc_range = range(df$Conc, na.rm = TRUE),
    subjects = sort(unique(as.character(df$ID))),
    timepoints = sort(unique(df$Time)),
    has_drug = has_drug
  )

  if (has_drug) {
    summary_info$drugs <- sort(unique(as.character(df$Drug)))
    summary_info$n_drugs <- length(summary_info$drugs)
  }

  return(summary_info)
}
