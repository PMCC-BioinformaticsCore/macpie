#' Collapse Replicates by Differential Expression
#'
#' This function collapses replicate samples in a tidyseurat object based on 
#' differential expression (DE) results. It returns a new tidyseurat object 
#' with a 'DE' assay containing the aggregated values.
#'
#' @param data A tidyseurat object with DE results
#' @param metric_col A character string indicating the column name in each 
#' DE dataframe to use for aggregation. Defaults to "metric".
#'
#' @details
#' This function performs the following steps:
#' - Converts the list of DE comparisons into a gene-by-condition matrix.
#' - Aggregates metadata across replicates by `combined_id`.
#' - Creates a new Seurat object with this matrix as a "DE" assay
#'
#' @return A new tidyseurat object with collapsed metadata 
#' @import tidyverse 
#' @import Seurat 
#' @importFrom Matrix Matrix
#' @importFrom tibble column_to_rownames
#' @import tidyseurat
#' @examples
#' data(mini_mac)
#' mac_collapsed <- aggregate_by_de(mini_mac)
#'
#' @export
aggregate_by_de <- function(data, metric_col = "metric") {
  
  validate_inputs <- function(de_list) {
    if (!inherits(de_list, "list") || length(de_list) == 0) {
      stop("Error: There are no DE comparisons, please run compute_multi_de.")
    }
    if (!inherits(metric_col, "character") || nchar(metric_col) == 0) {
      stop("Error: Parameter metric_col should be acharacter.")
    }
    required_cols <- c("gene", "combined_id", metric_col)
    if (!all(sapply(de_list, function(x) all(required_cols %in% colnames(x))))) {
      stop("Each element of the DE list must contain 'gene', 'combined_id', and the specified metric column.")
    }
  }
  
  # 1. Extract DE list and convert to wide matrix
  de_list <- data@tools$diff_exprs
  validate_inputs(de_list)
  
  de_matrix <- bind_rows(de_list) %>%
    select(gene, combined_id, all_of(metric_col)) %>%
    pivot_wider(names_from = combined_id, values_from = all_of(metric_col)) %>%
    column_to_rownames("gene") %>%
    as.matrix()
  
  # 2. Aggregate metadata by combined_id
  metadata <- data@meta.data %>%
    select(-starts_with("nCount"), -starts_with("nFeature"), -starts_with("percent"))
  
  metadata_collapsed <- metadata %>%
    group_by(combined_id) %>%
    summarise(
      across(
        everything(),
        ~ {
          u <- unique(.x)
          if (length(u) == 1) {
            u
          } else {
            switch(
              typeof(.x),
              character = NA_character_,
              integer   = mean(as.integer(as.character(.x)), na.rm = TRUE),
              double    = mean(.x, na.rm = TRUE),
              logical   = NA,              # Or use NA_logical_
              NA                             # fallback for unknown types
            )
          }
        },
        .names = "{.col}"
      ),
      .groups = "drop"
    ) %>%
    mutate(cell_id = combined_id) %>%
    column_to_rownames("cell_id") %>%
    select(-orig.ident)
  
  # 3. Create new Seurat object with DE matrix as RNA assay
  seurat_new <- CreateSeuratObject(
    counts = Matrix::Matrix(de_matrix, sparse = TRUE),
    data = Matrix::Matrix(de_matrix, sparse = TRUE),
    assay = "DE"
  )
  
  seurat_new <- seurat_new %>%
    inner_join(metadata_collapsed,by = c(".cell" = "combined_id")) 
  
  # 4. Add all tools slots from original object
  for (tool_name in names(data@tools)) {
    seurat_new@tools[[tool_name]] <- data@tools[[tool_name]]
  }
  
  return(seurat_new)
}
