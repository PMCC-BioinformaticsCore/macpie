#'Check and clean the metadata file generated from findmetadata function
#'
#'This function validates a metadata dataset, cleans specific fields, checks for
#'issues, and generates a summary table grouped by \code{Plate_ID}.
#'@author Mark Li
#'@description It checks the column completeness and removes special characters
#'  such as spaces and comma except from Treatment_1
#'@description It replaces special characters such as spaces and comma with
#'  underscore
#'@description It generates a quick summary of the different treatment and
#'  sample groups in the metadata
#'@param metadata this is the output from findmetadata function
#'@import dplyr
#'@return A list containing: \item{cleaned_metadata}{The cleaned metadata data
#'  frame after validation and modifications.}
#' \item{summary_table}{A summary table grouped by \code{Plate_ID}, showing the count of unique values
#' in selected columns.}
#'
#'@examples
#' \dontrun{
#' # Example CSV file
#' file_path <- system.file("/extdata/PMMSq033_metadata.csv", package = "macpie")
#' metadata <- read_metadata(file_path)
#' metadt_qc <- validate_metadata(metadata)
#' }
#'
#'@export
validate_metadata <- function(metadata) {
  if (is.null(metadata) || nrow(as.data.frame(metadata)) == 0) {
    stop("Error: Input metadata is empty or wrong format. Ensure the metadata is correctly provided.")
  }

  # Check Well_ID pattern
  check_well_id <- function(well_ids) {
    pattern <- "^[A-P](0[1-9]|1[0-9]|2[0-4])$"
    valid <- grepl(pattern, well_ids)
    return(valid)
  }

  # Initialize an empty list to store validation results
  issues <- list()

  # 1. Check for empty rows and summarize
  empty_rows <- apply(metadata, 2, function(column) sum(is.na(column) | column == ""))

  for (col_name in names(empty_rows)) {
    if (empty_rows[col_name] > 0) {
      issues[[col_name]] <- paste("Contains", empty_rows[col_name], "empty rows")
    }
  }

  # 2. Validate specific columns
  metadata$Plate_ID <- as.character(metadata$Plate_ID)
  metadata$Row <- as.character(metadata$Row)
  metadata$Column <- as.character(metadata$Column)
  metadata$Well_ID <- as.character(metadata$Well_ID)
  metadata$Time <- as.character(metadata$Time)
  metadata$Concentration_1 <- as.character(metadata$Concentration_1)

  # Check Plate_ID: Replace spaces and special characters with underscores
  metadata$Plate_ID <- gsub("[^[:alnum:]_]", "_", metadata$Plate_ID)

  # Validate individual columns
  validate_columns <- function(metadata) {

    # Validate Row
    if (any(!metadata$Row %in% LETTERS[1:16])) {
      issues[["Row"]] <- "Contains values outside A to P."
    }

    # Validate Column
    if (any(!grepl("^([1-9]|1[0-9]|2[0-4])$", metadata$Column))) {
      issues[["Column"]] <- "Contains values outside 1 to 24."
    }

    # Validate Well_ID
    if (any(!check_well_id(metadata$Well_ID))) {
      issues[["Well_ID"]] <- "Contains invalid Well_ID values."
    }

    # Validate numeric columns
    numeric_columns <- c("Time", "Concentration_1")
    for (col in numeric_columns) {
      if (any(!grepl("^\\d+(\\.\\d+)?$", metadata[[col]]))) {
        issues[[col]] <- paste("Contains non-numeric values in", col)
      }
    }

    return(issues)
  }

  issues <- append(issues, validate_columns(metadata))


  # Check all other columns: Should be character strings without special characters
  other_columns <- setdiff(names(metadata),
                           c("Plate_ID", "Row", "Column", "Well_ID", "Time", "Concentration_1", "Treatment_1"))
  for (col_name in other_columns) {
    if (any(grepl("[^[:alnum:]_ ]", metadata[[col_name]]))) {
      issues[[col_name]] <- "Contains special characters."
    }
    # Remove spaces in all other columns
    metadata[[col_name]] <- gsub("\\s", "", metadata[[col_name]])
  }

  # Special handling for Treatment_1
  if ("Treatment_1" %in% names(metadata)) {
    metadata$Treatment_1 <- gsub("[,\\s]", "_", metadata$Treatment_1)
  }

  # Print issues
  if (length(issues) > 0) {
    cat("\nValidation Issues:\n")
    for (issue_col in names(issues)) {
      cat(issue_col, ": ", issues[[issue_col]], "\n", sep = "")
    }
  } else {
    cat("\nNo validation issues found. Metadata is clean.\n")
  }

  # 3. Generate a summary table grouped by Plate_ID
  summary_columns <- c("Species", "Cell_type", "Model_type", "Time", "Unit",
                       "Treatment_1", "Concentration_1", "Unit_1", "Sample_type")

  cat("\nGenerating summary table grouped by Plate_ID...\n")

  summary_table <- metadata %>%
    dplyr::group_by(Plate_ID) %>%
    dplyr::summarize(across(all_of(summary_columns), ~ length(unique(.)), .names = "count_{col}")) %>%
    dplyr::ungroup()

  print(as.data.frame(summary_table))
}
