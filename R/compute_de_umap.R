
#' Prepare DE-based umap
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr pivot_wider
#' @importFrom stats dist
#' @returns A tidyseurat object with umap_de data frame in slot tools
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac <- compute_de_umap(mini_mac)

compute_de_umap <- function(data = NULL) {
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop(
      "compute_de_umap(): the following package is required but not installed: umap",
      "\nPlease install via `install.packages()`.")
  }

  # Helper function to validate input data
  validate_inputs <- function(de_list) {
    if (!inherits(de_list, "list") && length(de_list) > 0) {
      stop("argument 'data' must contain a list of DE comparisons in the slot tool.")
    }
  }

  de_list <- data@tools$diff_exprs
  validate_inputs(de_list)

  df <- bind_rows(de_list)
  df_wide <- df %>%
    select("gene", "combined_id", "metric") %>%
    pivot_wider(names_from = "combined_id", values_from = "metric")

  set.seed(1)
  df_umap <- umap::umap(t(df_wide[, -1]), method = "naive")

  # Get UMAP coordinates
  umap_mat <- df_umap$layout
  colnames(umap_mat) <- c("UMAPde_1", "UMAPde_2")
  
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
