# Data Validation and Preprocessing Module

library(readxl)
library(tools)

#' Validate and load PK data from uploaded file
#'
#' Expects wide format: first column is Time, each subsequent column
#' contains concentrations for one animal/subject. Column headers become
#' subject IDs (e.g., "Animal1", "Dog_01", "1").
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

#' Identify the Time column and reshape wide-format data to long format
#'
#' @param df Data frame in wide format (Time + animal columns)
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

  time_col_name <- names(df)[time_col_idx]
  conc_col_names <- names(df)[-time_col_idx]

  if (length(conc_col_names) < 1) {
    return(list(data = NULL, error = "No concentration columns found. File needs Time plus at least one animal column."))
  }

  # Check that Time column is numeric
  time_vals <- suppressWarnings(as.numeric(df[[time_col_name]]))
  if (all(is.na(time_vals))) {
    return(list(data = NULL,
                error = paste0("The Time column ('", time_col_name, "') does not contain numeric values.")))
  }

  # Check that concentration columns are numeric
  for (cn in conc_col_names) {
    test_vals <- suppressWarnings(as.numeric(df[[cn]]))
    if (all(is.na(test_vals))) {
      return(list(data = NULL,
                  error = paste0("Column '", cn, "' does not contain numeric concentration values.")))
    }
  }

  # Reshape from wide to long
  long_list <- list()
  for (cn in conc_col_names) {
    subj_df <- data.frame(
      ID = cn,
      Time = time_vals,
      Conc = suppressWarnings(as.numeric(df[[cn]])),
      stringsAsFactors = FALSE
    )
    long_list[[cn]] <- subj_df
  }

  long_df <- do.call(rbind, long_list)
  rownames(long_df) <- NULL

  return(list(data = long_df, error = NULL))
}

#' Standardize PK data types and handle missing values
#'
#' @param df Data frame with ID, Time, Conc columns (long format)
#' @return Cleaned data frame
standardize_pk_data <- function(df) {
  df$ID <- as.factor(df$ID)
  df$Time <- as.numeric(df$Time)
  df$Conc <- as.numeric(df$Conc)

  # Remove rows with NA in critical columns
  n_before <- nrow(df)
  df <- df[!is.na(df$Time) & !is.na(df$Conc), ]
  n_after <- nrow(df)

  if (n_before > n_after) {
    message(paste("Removed", n_before - n_after, "rows with missing Time or Conc values."))
  }

  # Sort by ID then Time
  df <- df[order(df$ID, df$Time), ]
  rownames(df) <- NULL

  return(df)
}

#' Get summary statistics for uploaded data
#'
#' @param df Validated PK data frame (long format)
#' @return Named list of summary info
pk_data_summary <- function(df) {
  list(
    n_subjects = length(unique(df$ID)),
    n_observations = nrow(df),
    time_range = range(df$Time, na.rm = TRUE),
    conc_range = range(df$Conc, na.rm = TRUE),
    subjects = sort(unique(as.character(df$ID))),
    timepoints = sort(unique(df$Time))
  )
}
