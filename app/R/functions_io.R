# ============================================================
# functions_io.R
# Input/output functions for the NGC-Ucar Toolkit
# ============================================================


read_ngc_toolkit_dataset <- function(filepath) {
  
  # ----------------------------------------------------------
  # 1. Basic file checks
  # ----------------------------------------------------------
  
  if (is.null(filepath) || !file.exists(filepath)) {
    stop("Dataset file not found.")
  }
  
  lines <- readLines(
    filepath,
    warn = FALSE,
    encoding = "UTF-8"
  )
  
  if (length(lines) < 4) {
    stop(
      paste(
        "Invalid NGC-Ucar Toolkit dataset.",
        "The file must contain at least four lines:",
        "1) dataset identifier,",
        "2) units,",
        "3) header sigma3,sigma1,",
        "4+) experimental data."
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 2. Read metadata
  # ----------------------------------------------------------
  
  dataset_id <- trimws(lines[1])
  units      <- trimws(lines[2])
  
  if (dataset_id == "") {
    stop("Invalid dataset: line 1 (dataset identifier) is empty.")
  }
  
  if (units == "") {
    stop("Invalid dataset: line 2 (units) is empty.")
  }
  
  
  # ----------------------------------------------------------
  # 3. Validate header
  # ----------------------------------------------------------
  
  header <- gsub("\\s+", "", lines[3])
  
  if (header != "sigma3,sigma1") {
    stop(
      paste0(
        "Invalid dataset header. ",
        "Line 3 must be exactly: sigma3,sigma1"
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 4. Read experimental data
  # ----------------------------------------------------------
  
  data <- tryCatch(
    {
      read.csv(
        filepath,
        skip = 2,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        strip.white = TRUE
      )
    },
    error = function(e) {
      stop(
        paste0(
          "The experimental data could not be read: ",
          e$message
        )
      )
    }
  )
  
  
  # ----------------------------------------------------------
  # 5. Validate table structure
  # ----------------------------------------------------------
  
  if (ncol(data) != 2) {
    stop(
      paste0(
        "Invalid dataset: exactly two data columns are required ",
        "(sigma3 and sigma1)."
      )
    )
  }
  
  if (!identical(names(data), c("sigma3", "sigma1"))) {
    stop(
      paste0(
        "Invalid column names. ",
        "The columns must be named sigma3 and sigma1."
      )
    )
  }
  
  if (nrow(data) == 0) {
    stop("Invalid dataset: no experimental data were found.")
  }
  
  
  # ----------------------------------------------------------
  # 6. Validate numeric values
  # ----------------------------------------------------------
  
  if (!is.numeric(data$sigma3) || !is.numeric(data$sigma1)) {
    stop(
      "Invalid dataset: sigma3 and sigma1 must contain numeric values only."
    )
  }
  
  if (anyNA(data$sigma3) || anyNA(data$sigma1)) {
    stop(
      "Invalid dataset: missing or non-numeric values were detected."
    )
  }
  
  if (any(!is.finite(data$sigma3)) ||
      any(!is.finite(data$sigma1))) {
    stop(
      "Invalid dataset: all stress values must be finite numbers."
    )
  }
  
  
  # ----------------------------------------------------------
  # 7. Mechanical consistency
  # ----------------------------------------------------------
  
  invalid_rows <- which(data$sigma1 < data$sigma3)
  
  if (length(invalid_rows) > 0) {
    stop(
      paste0(
        "Invalid dataset: sigma1 must be greater than or equal to sigma3. ",
        "Problem detected in row(s): ",
        paste(invalid_rows, collapse = ", "),
        "."
      )
    )
  }
  
  
  # ----------------------------------------------------------
  # 8. Return validated dataset
  # ----------------------------------------------------------
  
  return(
    list(
      dataset_id = dataset_id,
      units = units,
      n_tests = nrow(data),
      data = data
    )
  )
}