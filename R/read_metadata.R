#' Read Metadata from a File
#'
#' This function reads metadata from a specified file, validates the file
#' format, and extracts columns that match a predefined set of column names.
#'
#' @param file_path A character string specifying the path to the metadata file.
#'   The file must be in CSV or Excel format.
#' @param header Logical. If `TRUE`, the first row of the file is used as column
#'   names. Defaults to `TRUE`.
#' @param sep A character string specifying the column separator for CSV files.
#'   Defaults to `","`.
#' @param string_as_factors Logical. If `TRUE`, character columns in the data
#'   frame are converted to factors. Defaults to `FALSE`.
#' @param predefined_columns A character vector of column names to match in the
#'   file. Defaults to a predefined set of column names:
#'   \code{c("Plate_ID", "Well_ID", "Row", "Column", "Species", "Cell_type", "Model_type",
#'   "Time", "Unit", "Treatment_1", "Concentration_1", "Unit_1", "Sample_type", "Barcode", "Project")}.
#' @import readxl
#' @importFrom utils read.csv
#' @return A data frame containing only the matched columns if at least one
#'   predefined column is found. Returns `NULL` and prints an error message if
#'   no predefined columns are found.
#' @export
#'
#' @examples
#' # Example CSV file
#' file_path <- system.file("/extdata/PMMSq033_metadata.csv", package = "macpie")
#' result <- read_metadata(file_path)
#'
#' @details The function first checks if the file exists and validates its
#'   format (CSV or Excel). It then attempts to read the file and match its
#'   column names to a predefined set of column names. Supported file formats
#'   include:
#' \itemize{
#'   \item CSV files with `.csv` extension.
#'   \item Excel files with `.xls` or `.xlsx` extensions (requires the \code{readxl} package).
#' }

read_metadata <- function(file_path, header = TRUE, sep = ",", string_as_factors = FALSE,
                          predefined_columns = c("Plate_ID", "Well_ID", "Row", "Column", "Species",
                                                 "Treatment_1", "Concentration_1", 
                                                 "Sample_type", "Barcode")) {
  #check file exists
  if (!file.exists(file_path)) {
    stop("The file does not exist. Please provide a valid file path.")
  }

  # Read the metadata file based on its format
  read_file <- function(file_path) {
    tryCatch({
      if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
        data <- read.csv(file_path)
      } else if (grepl("\\.xls[x]?$", file_path, ignore.case = TRUE)) {
        data <- readxl::read_excel(file_path)
      } else {
        stop("Unsupported file format. Please select a CSV or Excel file.")
      }
      return(data)
    }, error = function(e) {
      stop("Error reading the file: ", e$message)
    })
  }

  data <- read_file(file_path) %>%
    mutate(across(
      where(is.numeric),
      ~ if (n_distinct(.) < 10) as.factor(.) else .  # or < 20 if you want
    ))

  # Match predefined columns with user columns
  # Strict string matching required here 
  match_columns <- function(data, predefined_columns) {
    user_columns <- colnames(data)
    matched_columns <- intersect(predefined_columns, user_columns)
    return(matched_columns)
  }

  # Match columns
  matched_columns <- match_columns(data, predefined_columns)

  #fix return outside a function error
  check_columns <- function(data, matched_columns) {
    if (length(matched_columns) == length(predefined_columns)) {
      cleaned_data <- data %>%
        arrange("Barcode")
      return(cleaned_data)
    } else {
      cat("Error: Some required columns were missing. Check your file.\n")
      cat("Missing columns:\n", paste(setdiff(predefined_columns, matched_columns), collapse = ", "), "\n")
      return(NULL)
    }
  }
  output <- check_columns(data = data, matched_columns = matched_columns)
  return(output)
}
