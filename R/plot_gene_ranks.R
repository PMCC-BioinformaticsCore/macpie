#' Generate a knee plot
#' A knee plot to show total number of total read counts for each gene in a given
#' treatment group
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param group_by A column that specifies the treatment group in the input data
#' @param samples Treatment group
#' @param scale Boolean statement for log10 transformation on both axes
#' @import dplyr
#' @importFrom ggplot2 ggplot aes geom_point
#' @return A ggplot object with a knee plot
#' @export
#' @examples
#' data("mini_mac")
#' p <- plot_gene_ranks(mini_mac,group_by = "combined_id", 
#' samples = "Staurosporine_10", scale = TRUE)



plot_gene_ranks <- function(data = NULL,
                            group_by = NULL,
                            samples = NULL,
                            scale = TRUE) {
  # Helper function to validate input data
  validate_inputs <- function(data, group_by, samples, scale) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    if (is.null(samples)) {
      stop("Missing the vectors of samples.")
    }
    if (!is.logical(scale) || length(scale) != 1) {
      stop("Error: 'scale' must be a single logical value (TRUE or FALSE).")
    }
  }
  validate_inputs(data, group_by, samples, scale)
  count_table <- as.data.frame(data@assays$RNA$count)
  combined_id_barcodes <- data@meta.data[[group_by]]
  names(combined_id_barcodes) <- rownames(data@meta.data)
  colnames(count_table) <- combined_id_barcodes[colnames(count_table)]
  if (samples != "all") {
    count_table_subset <- as.data.frame(count_table[, colnames(count_table) %in% samples])
    colnames(count_table_subset) <- paste0(samples, seq_len(ncol(count_table_subset)))
    rank_counts_genes <-
      count_table_subset %>%
      mutate(sum = rowSums(across(where(is.numeric)))) %>%
      arrange(-sum)
    rank_counts_genes$rank <- seq_len(nrow(rank_counts_genes))
  } else {
    colnames(count_table) <- paste0(colnames(count_table), "_", seq_len(ncol(count_table)))
    rank_counts_genes <-
      count_table %>%
      mutate(sum = rowSums(across(where(is.numeric)))) %>%
      arrange(-sum)
    rank_counts_genes$rank <- seq_len(nrow(rank_counts_genes))
  }
  if (scale){
    p <- ggplot(rank_counts_genes, aes(x = rank, y = sum)) +
      geom_point(shape = 1) +
      scale_y_log10() +
      scale_x_log10() +
      ylab("Read counts per gene") +
      xlab("Rank") + 
      theme_bw()
  } else {
    p <- ggplot(rank_counts_genes, aes(x = rank, y = sum)) +
      geom_point(shape = 1) +
      ylab("Read counts per gene") +
      xlab("Rank") + 
      theme_bw()
  }
  return(p)
}
