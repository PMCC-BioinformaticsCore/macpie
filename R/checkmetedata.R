
#' Check and clean the metadata file generated from findmetadata function
#' @author Mark Li
#' @description As said in title
#' @description It checks the column completeness and removes special characters such as spaces and comma except from Treatment_1
#' @description It replaces special characters such as spaces and comma with underscore
#' @description It generates a quick summary of the different treatment and sample groups in the metadata
#' @param metadata this is the output from findmetadata function
#' @return it prints out a summary
#'
#' @examples metadt_qc <- validate_metadata(meta_dt)
#'
validate_metadata <- function(metadata) {
  if (is.null(metadata)) {
    stop("Error: Input metadata is NULL. Ensure the metadata is correctly provided.")
  }

  # Helper function to check Well_ID pattern
  check_well_id <- function(well_ids) {
    pattern <- "^[A-P](0[1-9]|1[0-9]|2[0-4])$"
    valid <- grepl(pattern, well_ids)
    return(valid)
  }

  # Initialize an empty list to store validation results
  issues <- list()

  # 1. Check for empty rows and summarize
  empty_rows <- apply(metadata, 2, function(column) sum(is.na(column) | column == ""))
  total_rows <- nrow(metadata)

  cat("\nColumn-wise row presence summary:\n")
  for (col_name in names(empty_rows)) {
    non_empty <- total_rows - empty_rows[col_name]
    cat(col_name, ": ", non_empty, " non-empty rows out of ", total_rows, "\n", sep = "")
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

  # Check Row: Should be A to P
  if (any(!metadata$Row %in% LETTERS[1:16])) {
    issues[["Row"]] <- "Contains values outside A to P."
  }

  # Check Column: Should be 1 to 24
  if (any(!grepl("^([1-9]|1[0-9]|2[0-4])$", metadata$Column))) {
    issues[["Column"]] <- "Contains values outside 1 to 24."
  }

  # Check Well_ID: Should follow A01, A02, ..., P24
  if (any(!check_well_id(metadata$Well_ID))) {
    issues[["Well_ID"]] <- "Contains invalid Well_ID values."
  }

  # Check Time: Should be numeric
  if (any(!grepl("^\\d+(\\.\\d+)?$", metadata$Time))) {
    issues[["Time"]] <- "Contains non-numeric values."
  }

  # Check Concentration_1: Should be numeric
  if (any(!grepl("^\\d+(\\.\\d+)?$", metadata$Concentration_1))) {
    issues[["Concentration_1"]] <- "Contains non-numeric values."
  }

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

  print(summary_table)

  # Return cleaned metadata and summary table
  return(list(cleaned_metadata = metadata, summary_table = summary_table))
}



