
#' Prepare DE-based umap
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr pivot_wider
#' @importFrom umap umap
#' @returns A tidyseurat object with umap_de data frame in slot tools
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac <- compute_de_umap(mini_mac)

compute_de_umap <- function(data = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(de_list) {
    if (!inherits(de_list, "list") && length(de_list) > 0) {
      stop("Error: argument 'data' must contain a list of DE comparisons in the slot tool.")
    }
  }

  de_list <- data@tools$diff_exprs
  validate_inputs(de_list)

  df <- bind_rows(de_list)
  df_wide <- df %>%
    select("gene", "combined_id", "metric") %>%
    pivot_wider(names_from = "combined_id", values_from = "metric")

  set.seed(1)
  df_umap <- umap(t(df_wide[, -1]))

  # Get UMAP coordinates
  umap_mat <- df_umap$layout
  rownames(umap_mat) <- colnames(df_wide)  # these are the "combined_id"s
  
  # Create DimReduc object
  umap_reduction <- CreateDimReducObject(
    embeddings = umap_mat,
    key = "UMAPde_",
    assay = DefaultAssay(data)
  )
  
  # Add it to the Seurat object under reductions
  data@reductions[["umap_de"]] <- umap_reduction
  return(data)
}
