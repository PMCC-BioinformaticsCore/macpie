utils::globalVariables(c("cell", "combined_id", "n_replicates"))
#' Filter genes by expression and grouping
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param group_by Name of the column that defines groups of replicates
#' @param min_counts Minimum number of reads per gene per group
#' @param min_samples Minimum number of samples in a group
#' @returns tidyseurat object
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac <- filter_genes_by_expression(mini_mac,
#'                                        group_by = "combined_id", 
#'                                        min_counts = 10, min_samples = 2)

filter_genes_by_expression <- function(data, 
                                       group_by = "combined_id",
                                       min_counts = 10,
                                       min_samples = 2) {
  # Check if metadata column exists
  if (!group_by %in% colnames(data@meta.data)) {
    stop(paste("Metadata column", group_by, "not found in Seurat object."))
  }
  
  # Extract expression matrix
  expr_mat <- GetAssayData(data, layer = "counts", assay = "RNA")
  
  # Get metadata and filter for nCount_RNA
  meta <- data@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, !!group_by)
  
  # Subset expression to those cells
  expr_mat <- expr_mat[, meta$cell, drop = FALSE]
  
  # Add combined_id info
  combined_ids <- meta[[group_by]]
  names(combined_ids) <- meta$cell
  
  # Melt the matrix to long format
  expr_long <- as.data.frame(as.matrix(expr_mat)) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expr") %>%
    mutate(combined_id = combined_ids[cell])
  
  # Count expression
  expr_summary <- expr_long %>%
    group_by(gene, combined_id) %>%
    summarise(n_replicates = sum(expr > min_counts), .groups = "drop") %>%
    filter(n_replicates >= min_samples)
  
  # Keep only genes passing the filter in any group
  genes_to_keep <- unique(expr_summary$gene)
  
  data <- subset(data, features = genes_to_keep)
  # Filter all assay slots (counts, data, scale.data)
  for (slot in names(data@assays[["RNA"]]@layers)) {
    if (!is.null(data@assays[["RNA"]][slot]) &&
      nrow(data@assays[["RNA"]][slot]) > 0) {
      data@assays[["RNA"]][slot] <- data@assays[["RNA"]][slot][genes_to_keep, , drop = FALSE]
    }
  }
  
  # update both nFeature_RNA and nCount_RNA columns 
  data$nFeature_RNA <- Matrix::colSums(data@assays$RNA$counts>0)
  data$nCount_RNA <- Matrix::colSums(data@assays$RNA$counts)
  return(data)
}
