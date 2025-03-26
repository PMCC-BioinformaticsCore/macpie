utils::globalVariables(c(".", ".data", "gene", "log2FC", "metric"))
#' Generate heatmap of DE genes from multiple treatments
#' This is the function to generate a heatmap of DE genes from running compute_multi_DE
#' that shared by more than one treatment group. There are a few options available
#' to help you to extract shared DE genes.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param group_by A column that specifies the treatment group in the input data
#' @param value To use logCPM, log2FC or metric in the heatmap
#' @param p_value_cutoff Cutoff for adjusted p-value (column p_value_adj), default 0.01
#' @param direction Direction to select up or down regulated genes or in both directions
#' @param n_genes Top n genes to be extracted from each treatment comparison
#' @param control The control group to be included in the final heatmap, usually DMSO_0
#' @param by Extract top n genes by either absolute fold change or by adjusted p-value
#' @param gene_list External list of genes to plot the heatmap on
#' @import pheatmap
#' @import dplyr
#' @importFrom purrr map2_dfr
#' @returns a pheatmap object
#' @export
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_multi_de(mac, group_by = "combined_id", value = "lcpm", p_value_cutoff = 0.01, direction="up", n_genes = 10, control = "DMSO_0", by="fc")


plot_multi_de <- function(data = NULL,
                          group_by = NULL,
                          value = NULL,
                          p_value_cutoff = 0.01,
                          direction = "both",
                          n_genes = 10,
                          control = "DMSO_0",
                          by = "fc",
                          gene_list = NULL) {
  # Helper function to validate input data
  validate_inputs <- function(data, group_by, value, p_value_cutoff, direction, n_genes, control, by, gene_list) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    if (!is.null(value) && !value %in% c("lcpm", "log2FC", "metric"))
      if (!inherits(p_value_cutoff, "numeric")) {
        stop("Error: argument 'p_value_cutoff' must be numeric.")
      }


    if (!is.null(direction) && !direction %in% c("up", "down", "both")) {
      stop("Value of the direction paramater should be up, down or both. Default is both.")
    }
    if (!inherits(n_genes, "numeric")) {
      stop("Error: argument 'n_genes' must be numeric.")
    }
    if (!is.null(control) && !inherits(control, "character")) {
      stop("Error: argument 'control' must be control group in same format as combined_id.")
    }
    if (!is.null(by) && !by %in% c("fc", "adj_p_val")) {
      stop("Value of the by parameter show be either fc or adj_p_val.")
    }
    if (length(data@tools$diff_exprs) == 0) {
      stop("Missing information on DE genes. Run compute_multi_de first.")
    }
    if (!is.null(gene_list) && !inherits(gene_list, "character")) {
      stop("Error: argument 'gene_list' must be a vector of characters.")
    }
  }
  validate_inputs(data, group_by, value, p_value_cutoff, direction, n_genes, control, by, gene_list)
  #extract de information from the object data
  all_de <- data@tools$diff_exprs
  de_df <- unique(map2_dfr(all_de, names(all_de), ~ mutate(.x, source = .y)))


  filtered_de_df <- de_df %>%
    filter(.data$p_value_adj < .env$p_value_cutoff) %>%
    filter(direction == "both" | (direction == "up" & .data$log2FC > 0) | (direction == "down" & .data$log2FC < 0))
  top_genes_per_combined_id <- filtered_de_df %>%
    group_by(.data[[group_by]]) %>% {
      if (by == "fc") {
        slice_max(., order_by = abs(.data$log2FC), n = n_genes)
      } else {
        slice_min(., order_by = .data$p_value_adj, n = n_genes)
      }
    }


  #find genes shared by at least 2 different treatment groups
  common_genes <- top_genes_per_combined_id %>%
    group_by(.data$gene) %>%
    filter(n_distinct(.data[[group_by]]) > 1)
  if (length(common_genes) == 0) {
    stop("No genes are shared by at least 2 different treatment groups.")
  }
  if(is.null(gene_list)){
    features <- unique(common_genes$gene)
  } else {
    features <- unique(gene_list)
  }


  common_genes_treatments <- common_genes %>% dplyr::distinct(.data[[group_by]]) %>% pull() %>% unique()
  #add control group
  common_genes_treatments <- c(common_genes_treatments, control)
  if (value == "lcpm") {
    #calculate cpm on full samples
    lcpm <- cpm(data@assays$RNA$counts, log = TRUE)
    combined_id_barcodes <- data@meta.data[[group_by]]
    names(combined_id_barcodes) <- rownames(data@meta.data)
    colnames(lcpm) <- combined_id_barcodes[colnames(lcpm)]
    sub_lcpm <- lcpm[features, colnames(lcpm) %in% common_genes_treatments]
    #plot lcpm
    p <- pheatmap(sub_lcpm,
                  cexRow = 0.1,
                  cexCol = 0.2,
                  col = macpie_colours$continuous)
  } else if (value == "log2FC") {
    # plot log2fc
    fc_common_genes <- de_df %>% 
      filter(gene %in% features) %>%
      select(gene, log2FC, source) %>%
      pivot_wider(names_from = source, values_from = log2FC, values_fill = 0)
    fc_matrix <- as.matrix(fc_common_genes[, -1])
    rownames(fc_matrix) <- fc_common_genes$gene 
    storage.mode(fc_matrix) <- "numeric"
    p <- pheatmap(fc_matrix,
                  cexRow = 0.1,
                  cexCol = 0.2,
                  col = rev(macpie_colours$continuous),
                  cluster_cols = F)
  } else {
    # plot metric
    fc_common_genes <- de_df %>% 
      filter(gene %in% features) %>%
      select(gene, metric, any_of(group_by)) %>%
      pivot_wider(names_from = any_of(group_by), values_from = metric, values_fill = 0)
    fc_matrix <- as.matrix(fc_common_genes[, -1])
    rownames(fc_matrix) <- fc_common_genes$gene
    storage.mode(fc_matrix) <- "numeric
    p <- pheatmap(fc_matrix,
                  cexRow = 0.1,
                  cexCol = 0.2,
                  col = macpie_colours$continuous)
  }
  return(p)
}

