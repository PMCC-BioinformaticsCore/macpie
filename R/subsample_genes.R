#' Subsample genes (fast helper function for zero-inflation checks)
#' 
#' @description
#' Quickly subsample a specified number of **genes** from a Seurat object and
#' return a smaller Seurat object containing the selected features and all
#' original wells/samples. This is a lightweight convenience wrapper around
#' **`seqgendiff::select_counts()`** and is intended for creating a small
#' working object to run **`check_zeroinflation()`** (or similar
#' diagnostics) rapidly.
#' 
#' @param data   A Seurat object (v4 or v5) with counts in assay "RNA".
#' @param ngene  Integer. Number of genes to keep (must be <= total genes).
#' @param gselect Gene-selection strategy as used by
#'   **`seqgendiff::select_counts()`**. Defaults to `"random"`.
#' @param seed   Integer random seed for reproducibility.
#' @return A Seurat object containing the subsampled genes and all original wells/samples.
#' @export
#' @examples
#' data(mini_mac)
#' subsample_genes(mini_mac, ngene = 50, gselect = "random", seed = 1 )
#' 
subsample_genes <- function(data, 
                            ngene = 100, 
                            gselect = "random", 
                            seed = 1){
  if (utils::packageVersion("Seurat") < "5.0.0"
  ) {
    count_matrix <- GetAssayData(data, assay = "RNA", slot = "counts")
  } else {
    count_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
  }
  count_matrix <- as.matrix(count_matrix)
  if (!inherits(data, "Seurat")) {
    stop("argument 'data' must be a Seurat or TidySeurat object.")
  }
  if (ngene > nrow(count_matrix)){
    stop("ngene should be less than the total number of genes in the dataset.")
  }
  set.seed(seed)
  sub_matrix <- seqgendiff::select_counts(count_matrix, ngene = ngene, gselect = gselect)
  subsample_genes <- rownames(sub_matrix)
  data_sub <- data[subsample_genes, ]
  return(data_sub)
  
}
