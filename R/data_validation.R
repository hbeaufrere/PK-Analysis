# Data Validation and Preprocessing Module

library(readxl)
library(tools)

#' Validate and load PK data from uploaded file
#'
#' @param file_path Path to uploaded file
#' @param file_name Original file name (for extension detection)
#' @return List with 'data' (data.frame or NULL) and 'error' (character or NULL)
load_pk_data <- function(file_path, file_name) {
  ext <- tolower(file_ext(file_name))

  tryCatch({
    if (ext == "csv") {
      df <- read.csv(file_path, stringsAsFactors = FALSE)
    } else if (ext %in% c("xls", "xlsx")) {
      df <- as.data.frame(readxl::read_excel(file_path))
    } else {
      return(list(data = NULL, error = "Unsupported file format. Please upload a CSV or XLS/XLSX file."))
    }

    validation <- validate_pk_columns(df)
    if (!is.null(validation$error)) {
      return(validation)
    }

    df <- standardize_pk_data(validation$data)
    return(list(data = df, error = NULL))

  }, error = function(e) {
    return(list(data = NULL, error = paste("Error reading file:", e$message)))
  })
}

#' Check that required columns exist (case-insensitive matching)
#'
#' @param df Data frame to validate
#' @return List with 'data' and 'error'
validate_pk_columns <- function(df) {
  cols_lower <- tolower(names(df))
  required <- c("id", "time", "concentration")
  # Also accept common variants
  aliases <- list(
    id = c("id", "subject", "subj", "subject_id", "animal"),
    time = c("time", "time_h", "time_hr", "time_min", "timepoint"),
    concentration = c("concentration", "conc", "concentrations", "dv", "cp", "plasma_conc")
  )

  matched <- list()
  for (req in names(aliases)) {
    found <- FALSE
    for (alias in aliases[[req]]) {
      idx <- which(cols_lower == alias)
      if (length(idx) > 0) {
        matched[[req]] <- names(df)[idx[1]]
        found <- TRUE
        break
      }
    }
    if (!found) {
      return(list(
        data = NULL,
        error = paste0(
          "Required column '", req, "' not found. ",
          "Expected one of: ", paste(aliases[[req]], collapse = ", "), ". ",
          "Found columns: ", paste(names(df), collapse = ", ")
        )
      ))
    }
  }

  # Rename to standard names
  names(df)[names(df) == matched$id] <- "ID"
  names(df)[names(df) == matched$time] <- "Time"
  names(df)[names(df) == matched$concentration] <- "Conc"

  return(list(data = df, error = NULL))
}

#' Standardize PK data types and handle missing values
#'
#' @param df Data frame with ID, Time, Conc columns
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
#' @param df Validated PK data frame
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
