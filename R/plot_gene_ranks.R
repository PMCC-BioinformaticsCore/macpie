#' Generate a knee plot
#' A knee plot to show total number of total read counts for each gene in a given
#' treatment group
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param group_by A column that specifies the treatment group in the input data
#' @param samples Treatment group
#' @import dplyr
#' @import ggplot2
#' @export
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_gene_ranks(mac, group_by = "combined_id", samples = "Staurosporine_10", scale_y_log = TRUE)


plot_gene_ranks <- function(data = NULL,
                            group_by = NULL,
                            samples = NULL,
                            scale_y_log = TRUE) {
  # Helper function to validate input data
  validate_inputs <- function(data, group_by, samples, scale_y_log) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    if (is.null(samples)) {
      stop("Missing the vectors of samples.")
    }
    if (!is.logical(scale_y_log) || length(scale_y_log) != 1) {
      stop("Error: 'scale_y_log' must be a single logical value (TRUE or FALSE).")
    }
  }
  validate_inputs(data, group_by, samples, scale_y_log)
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
  if (scale_y_log){
    p <- ggplot(rank_counts_genes, aes(x = rank, y = log1p(sum))) +
      geom_point(shape = 1) +
      macpie_theme(show_x_title = TRUE, show_y_title = TRUE) 
  } else {
    p <- ggplot(rank_counts_genes, aes(x = rank, y = sum)) +
      geom_point(shape = 1) +
      macpie_theme(show_x_title = TRUE, show_y_title = TRUE)
  }
  return(p)
}
