
#' User interactive function to find the metadata file
#' @author Mark Li
#' @description This is a rigid function compared to findmetadata_liberal
#' @description The function ask user to select a metadatafile
#' @description It checks the read permission of the file and its format
#' @description Then it will read the user file and check column headers agaisnt the predefined columns
#' @description It will also provide a quick summary of matched columns identified
#' @description User should ideally have 14 matched columns
#' @importFrom dplyr readr readxl htt2
#' @param predefined_columns this is the most important and rigid input for the entire package
#' @param predefined_columns for findmetadata function, we are string matching the 14 columns from user input metadatafile
#' @return it returns a dataframe of metadafile
#'
#' @examples meta_dt <- findmetadata()


findmetadata <- function(predefined_columns = c("Plate_ID", "Well_ID", "Row", "Column", "Species",
                                                "Cell_type", "Model_type", "Time", "Unit",
                                                "Treatment_1", "Concentration_1", "Unit_1",
                                                "Sample_type", "Barcode")) {
  # User select a metadata file
  select_metadatafile <- function() {
    cat("Please select a metadata file in CSV or Excel format...\n")
    file_path <- file.choose()
    cat("You selected:", file_path, "\n")

    # Check if the file is readable
    if (!file.exists(file_path)) {
      stop("Error: The selected file does not exist. Please try again.")
    }
    if (!file.access(file_path, 4) == 0) {
      stop("Error: The selected file is not readable. Please check file permissions.")
    }
    return(file_path)
  }

  # Read the metadata file based on its format
  read_file <- function(file_path) {
    tryCatch({
      if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
        data <- readr::read_csv(file_path)
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

  # Match predefined columns with user columns
  # Strict string matching required here as requested by Xin

  match_columns <- function(data, predefined_columns) {
    user_columns <- colnames(data)
    matched_columns <- intersect(predefined_columns, user_columns)
    return(matched_columns)
  }


  file_path <- select_metadatafile()
  data <- read_file(file_path)

  # Print predefined column names
  cat("Predefined column names:\n", paste(predefined_columns, collapse = ", "), "\n")

  # Match columns
  matched_columns <- match_columns(data, predefined_columns)

  if (length(matched_columns) > 0) {
    cat("Matched columns:\n", paste(matched_columns, collapse = ", "), "\n")
    cat("Number of matched columns:\n", length(matched_columns), "\n")
    cleaned_data <- data[, matched_columns, drop = FALSE]
    return(cleaned_data)
  } else {
    cat("Error: No columns matching predefined criteria found. Check your file.\n")
    return(NULL)
  }
}


