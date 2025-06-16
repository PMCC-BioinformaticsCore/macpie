#' Compute chemical descriptors from SMILES
#'
#' This function parses SMILES strings and computes chemical descriptors using rcdk.
#' It stores cleaned, non-redundant descriptors in `tools$chem_descriptors`.
#'
#' @param data A tidyseurat object with a `smiles` column and `Treatment_1` column.
#' @param r_squared R squared value, default of 0.6
#' @returns The same tidyseurat object with a new entry in `tools$chem_descriptors`.
#' @importFrom rcdk parse.smiles get.desc.names eval.desc
#' @importFrom dplyr select where bind_rows mutate
#' @importFrom stringr str_replace_all str_trim str_to_title
#' @importFrom stats cor
#' @export
#'
#' @examples
#' mock_data <- tibble::tibble(
#' Treatment_1 = c("Aspirin", "Caffeine", "NonExistentCompound_123")
#' )
#' result <- compute_smiles(mock_data)
#' data <- compute_chem_descriptors(result)
compute_chem_descriptors <- function(data, r_squared = 0.6) {
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
  
  # Prepare compound names
  data$clean_compound_name <- data %>%
    select("Treatment_1") %>%
    pull() %>%
    str_replace_all("_", " ") %>%
    str_trim() %>%
    str_to_title()
  
  smiles_list <- data %>% select("clean_compound_name", "Treatment_1", "smiles") %>% distinct()
  
  # Parse SMILES and name molecules
  mol <- parse.smiles(smiles_list$smiles[!is.na(smiles_list$smiles)])
  names(mol) <- smiles_list$clean_compound_name[!is.na(smiles_list$smiles)]
  
  # Filter safe descriptors (exclude 3D/charge-based)
  safe_descs <- setdiff(
    get.desc.names("all"),
    grep("charge|peoe|ATS", get.desc.names("all"), value = TRUE)
  )
  
  # Evaluate descriptors
  descriptor_df <- bind_rows(lapply(mol, function(m) {
    if (!is.null(m)) {
      suppressWarnings(
        tryCatch(eval.desc(m, safe_descs), error = function(e) NULL)
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
  corr_mat <- cor(descriptor_df_clean %>% select(-"clean_compound_name"), use = "pairwise.complete.obs")^2
  upper_tri <- which(upper.tri(corr_mat) & corr_mat > r_squared, arr.ind = TRUE)
  cols_to_remove <- unique(colnames(corr_mat)[upper_tri[, 2]])
  
  descriptor_df_clean <- descriptor_df_clean %>%
    select(-all_of(cols_to_remove))
  
  descriptor_df_clean$Treatment_1 <- smiles_list$Treatment_1[!is.na(smiles_list$smiles)]
  
  # Store in @tools
  if (inherits(data, "tbl_df")) {
    data <- data %>%
      left_join(., descriptor_df_clean, join_by("Treatment_1"))
  } else {
    data@tools$chem_descriptors <- descriptor_df_clean
  }
  return(data)
}
