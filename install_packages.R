# Install required packages for PK Analysis App
# Run this script once before launching the app

required_packages <- c(
  "shiny",       # Web application framework
  "DT",          # Interactive data tables
  "ggplot2",     # Plotting
  "nlme",        # Non-linear mixed effects modeling
  "readxl",      # Reading Excel files
  "scales",      # Axis scale formatting
  "tools"        # File extension utilities (base R, usually included)
)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing", pkg, "..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message(paste(pkg, "is already installed."))
  }
}

invisible(lapply(required_packages, install_if_missing))

message("\nAll packages installed. You can now run the app with: shiny::runApp()")
