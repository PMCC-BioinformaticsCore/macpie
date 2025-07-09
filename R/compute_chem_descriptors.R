#' Compute chemical descriptors from SMILES
#'
#' This function parses SMILES strings and computes chemical descriptors using rcdk.
#' It stores cleaned, non-redundant descriptors in `tools$chem_descriptors`.
#'
#' @param data A tidyseurat object with a `smiles` column and `Treatment_1` column.
#' @param r_squared R squared value, default of 0.6
#' @param descriptors Specify any descriptors of interest from rcdk 
#' @returns The same tidyseurat object with a new entry in `tools$chem_descriptors`.
#' @importFrom dplyr select where bind_rows mutate
#' @importFrom stringr str_replace_all str_trim str_to_title
#' @importFrom stats cor
#' @export
#'
#' @examples
#' \dontrun{
#' mock_data <- tibble::tibble(
#'   Treatment_1 = c("Aspirin", "Caffeine", "NonExistentCompound_123")
#' )
#' result <- compute_smiles(mock_data,compound_column = "Treatment_1" )
#' data <- compute_chem_descriptors(result, descriptors =
#' "org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor")
#' }
compute_chem_descriptors <- function(data,
                                     r_squared = 0.6,
                                     descriptors = NULL) {
  if (inherits(data, "tbl_df")) {
    if (!"smiles" %in% colnames(data)) {
      stop("The input must contain a `smiles` column. Run compute_smiles() first.")
    }
    if (!"Treatment_1" %in% colnames(data)) {
      stop("The input tibble must contain a column named 'Treatment_1'")
    }
  } else {
    if (!"smiles" %in% colnames(data@meta.data)) {
      stop("The input must contain a `smiles` column. Run compute_smiles() first.")
    }
    if (!"Treatment_1" %in% colnames(data@meta.data)) {
      stop("The input object must contain a column named 'Treatment_1'")
    }
  }
  
  if (!inherits(r_squared, "numeric")) {
    stop("The r value must be numeric")
  }
  
  if (!is.null(descriptors) && !is.character(descriptors)) {
    stop("`descriptors` must be a character vector (or NULL to use the defaults)")
  }
  
  # Prepare compound names
  smiles_list <- data %>%
    select("Treatment_1", "smiles") %>%
    distinct() %>%
    mutate(clean_compound_name = str_to_title(
      str_trim(
        str_replace_all(.data$Treatment_1, "_", " ")
        )
      )
    )
  
  smiles_list <- smiles_list %>% filter(!is.na(.data$smiles))
  
  # Parse SMILES and name molecules
  safe_parse <- function(s) {
    tryCatch(rcdk::parse.smiles(s)[[1]], error = function(e) NULL)
  }
  
  mol <- lapply(smiles_list$smiles, safe_parse)
  names(mol) <- smiles_list$clean_compound_name
  mol <- Filter(Negate(is.null), mol)
  
  # Filter safe descriptors (exclude 3D/charge-based)
  safe_descs <- setdiff(
    rcdk::get.desc.names("all"),
    grep("charge|peoe|ATS", rcdk::get.desc.names("all"), value = TRUE)
  )
  
  # Select the descriptors of interest 
  if (!is.null(descriptors)){
    safe_descs <- intersect(safe_descs, descriptors)
  }
  
  # Evaluate descriptors
  descriptor_df <- bind_rows(lapply(mol, function(m) {
    if (!is.null(m)) {
      suppressWarnings(
        tryCatch(rcdk::eval.desc(m, safe_descs), error = function(e) NULL)
      )
    } else {
      NULL
    }
  }), .id = "clean_compound_name")
  
  # Clean descriptors: remove constant or NA columns
  descriptor_df_clean <- descriptor_df %>%
    select(where(~ all(!is.na(.)))) %>%
    select(where(~ n_distinct(.) > 1))
  
  # Remove highly correlated columns (R² > 0.6)
  if (ncol(descriptor_df_clean) > 1) {
    corr_mat <- cor(descriptor_df_clean %>% select(-"clean_compound_name"), use = "pairwise.complete.obs")^2
    upper_tri <- which(upper.tri(corr_mat) & corr_mat > r_squared, arr.ind = TRUE)
    cols_to_remove <- unique(colnames(corr_mat)[upper_tri[, 2]])
    
    descriptor_df_clean <- descriptor_df_clean %>%
      select(-all_of(cols_to_remove))
  }

  descriptor_df_clean$Treatment_1 <- smiles_list$Treatment_1
  
  # Store in @tools
  if (inherits(data, "tbl_df")) {
    data <- data %>%
      left_join(., descriptor_df_clean, join_by("Treatment_1"))
  } else {
    data@tools$chem_descriptors <- descriptor_df_clean
  }
  return(data)
}
