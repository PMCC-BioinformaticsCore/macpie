
#' User interactive function to find the metadata file with freedom
#' Yet to be decided by the dev team discussion
#' @author Mark Li
#' @description This is a liberal function compared to findmetadata
#' @description The function ask user to select a metadatafile
#' @description It checks the read permission of the file and its format
#' @description Then it will read the user file and check column headers
#' @description Then it will ask user to identify what column of the input file matches the predefined column one by one
#' @description ML can also modify the matching process to use a number vector
#' @description User should ideally have 14 matched columns
#' @param predefined_columns this is the most important and rigid input for the entire package
#' @return it returns a dataframe of metadafile
#'
#' @examples ML is too lazy to write
#' @examples testmetadata <- findmetadata_liberal()


findmetadata_liberal <- function(predefined_columns = c("Plate_ID", "Well_ID", "Row", "Column", "Species",
                                                "Cell_type", "Model_type", "Time", "Unit",
                                                "Treatment_1", "Concentration_1", "Unit_1",
                                                "Sample_type", "Barcode")) {
  # Step 1: Select a metadata file
  select_metadatafile <- function() {
    cat("Please select a metadata file in CSV or Excel format...\n")
    file_path <- file.choose()
    cat("You selected:", file_path, "\n")
    return(file_path)
  }

  # Step 2: Read the metadata file based on its format
  read_file <- function(file_path) {
    if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
      data <- readr::read_csv(file_path)
    } else if (grepl("\\.xls[x]?$", file_path, ignore.case = TRUE)) {
      data <- readxl::read_excel(file_path)
    } else {
      stop("Unsupported file format. Please select a CSV or Excel file.")
    }
    return(data)
  }

  # Step 3: Match column names interactively
  match_columns <- function(data, predefined_columns) {
    # Display predefined columns
    cat("Predefined columns:\n")
    print(predefined_columns)

    # Extract user column names
    user_columns <- colnames(data)
    cat("User-provided columns:\n")
    print(user_columns)

    # Interactive matching
    matched_columns <- character(length(predefined_columns))
    for (i in seq_along(predefined_columns)) {
      cat("\nMatching for predefined column:", predefined_columns[i], "\n")
      print(DT::datatable(data.frame(User_Columns = user_columns),
                          selection = "single", rownames = FALSE))
      selected_index <- as.integer(readline(prompt = "Enter the row number of the matching column: "))
      if (!is.na(selected_index) && selected_index > 0 && selected_index <= length(user_columns)) {
        matched_columns[i] <- user_columns[selected_index]
        cat("Matched predefined column", predefined_columns[i], "to user column", user_columns[selected_index], "\n")
      } else {
        cat("No match selected for", predefined_columns[i], "\n")
      }
    }

    # Return matching results
    return(data.frame(Predefined_Column = predefined_columns, Matched_Column = matched_columns))
  }

  # Main execution flow
  file_path <- select_metadatafile()
  data <- read_file(file_path)

  # Match columns
  matched_columns <- match_columns(data, predefined_columns)

  # Print summary
  cat("\nMatching Summary:\n")
  print(matched_columns)

  # Return results
  return(list(file_path = file_path, matched_columns = matched_columns))
}



