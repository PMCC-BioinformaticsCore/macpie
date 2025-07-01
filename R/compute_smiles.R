#' Annotate tidyseurat object with SMILES
#'
#' This function retrieves isomeric SMILES from PubChem based on compound names
#' in the `Treatment_1` column of a tidyseurat object.
#'
#' @param data A tidyseurat object
#' @param compound_column Column with the generic name of a compound 
#' @importFrom webchem get_cid pc_prop
#' @importFrom stringr str_replace_all str_trim str_to_title
#' @importFrom dplyr mutate left_join rename select pull
#' @returns A tidyseurat object with a `smiles` column added to the metadata.
#' @export
#'
#' @examples
#' mock_data <- tibble::tibble(
#' Treatment_1 = c("Aspirin", "Caffeine", "NonExistentCompound_123")
#' )
#' result <- compute_smiles(mock_data, compound_column = Treatment_1)

compute_smiles <- function(data, compound_column) {
  if (!compound_column %in% colnames(data@meta.data)) {
    stop("The compound_column does not exist in metadata.")
  } 
  
  # Clean compound names
  data$clean_compound_name <- data %>%
    pull(!!sym(compound_column)) %>%
    str_replace_all("_", " ") %>%
    str_trim() 
  
  # Retrieve compound CIDs from PubChem
  cids <- get_cid(unique(data$clean_compound_name), from = "name", match = "first", verbose = FALSE)
  
  # Retrieve SMILES and merge with CIDs
  smiles_df <- pc_prop(cids$cid, properties = "IsomericSMILES") %>%
    mutate("CID" = as.character(.data$CID)) %>%
    left_join(cids, by = c("CID" = "cid")) %>%
    rename(smiles = "IsomericSMILES")
  
  # Merge back into main object
  data <- data %>%
    left_join(., smiles_df, by = c("clean_compound_name" = "query")) %>%
    select(-"clean_compound_name")
  return(data)
}
