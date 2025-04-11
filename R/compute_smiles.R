#' Annotate tidyseurat object with SMILES
#'
#' This function retrieves isomeric SMILES from PubChem based on compound names
#' in the `Treatment_1` column of a tidyseurat object.
#'
#' @param data A tidyseurat object with a column named `Treatment_1`.
#' @importFrom webchem get_cid pc_prop
#' @importFrom stringr str_replace_all str_trim str_to_title
#' @importFrom dplyr mutate left_join rename select pull
#' @returns A tidyseurat object with a `smiles` column added to the metadata.
#' @export
#'
#' @examples
#' mac <- compute_smiles(mac)

compute_smiles <- function(data) {
  if (!"Treatment_1" %in% colnames(data)) {
    stop("The input object must contain a column named 'Treatment_1'")
  }
  
  # Clean compound names
  data$clean_compound_name <- data %>%
    select(Treatment_1) %>%
    pull() %>%
    str_replace_all("_", " ") %>%
    str_trim() 
  
  # Retrieve compound CIDs from PubChem
  cids <- get_cid(unique(data$clean_compound_name), from = "name", match = "first", verbose = FALSE)
  
  # Retrieve SMILES and merge with CIDs
  smiles_df <- pc_prop(cids$cid, properties = "IsomericSMILES") %>%
    mutate(CID = as.character(CID)) %>%
    left_join(cids, by = c("CID" = "cid")) %>%
    rename(smiles = IsomericSMILES)
  
  # Merge back into main object
  data <- data %>%
    left_join(., smiles_df, by = c("clean_compound_name" = "query")) %>%
    select(-clean_compound_name)
  return(data)
}