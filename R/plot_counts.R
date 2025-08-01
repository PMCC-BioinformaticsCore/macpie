utils::globalVariables(c(".data", "Treatment", "Expression", "spikes", "Genes", ".N", "Replicate"))
#' Generate a box plot to show gene expression (CPM)
#' 
#' This is the function to generate a box plot to show CPM levels of DE genes
#' among selected treatment samples and control samples.
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param genes Genes to be plotted
#' @param group_by A column that specifies the treatment group in the input data
#' @param treatment_samples Value in the column "combined_id" representing replicates of treatment samples in the data
#' @param control_samples Value in the column "combined_id"  representing replicates of control samples in the data
#' @param color_by A column that specifies the group coloring
#' @param normalisation One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch To indicate patch factor
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter facet_wrap scale_x_discrete scale_fill_manual labs
#' @import dplyr
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' data(mini_mac)
#' genes <- mini_mac@tools$diff_exprs$Staurosporine_10$gene[1:6]
#' p <- plot_counts(mini_mac, genes = genes, group_by = "combined_id", 
#' treatment_samples = "Staurosporine_10", 
#' control_samples = "DMSO_0",
#' normalisation = "clr")


plot_counts <- function(data = NULL,
                        genes = NULL,
                        group_by = NULL,
                        treatment_samples = NULL,
                        control_samples = NULL,
                        color_by = NULL,
                        normalisation = NULL,
                        batch = 1) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
      "plot_counts(): the following package is required but not installed: data.table",
      "\nPlease install via `install.packages()`.")
  }
  # Helper function to validate input data
  validate_inputs <- function(data, genes, group_by, treatment_samples, control_samples, color_by, normalisation) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (!inherits(genes, "character")) {
      stop("Error: arguemnt 'genes' must be characters.")
    }
    if (length(genes) > 20) {
      stop("Try plotting with less genes.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    if (is.null(treatment_samples) || is.null(control_samples)) {
      stop("Missing the vectors of treatment and control samples.")
    }
    normalisation <- if (is.null(normalisation)) "raw" else normalisation
    color_by <- if (is.null(color_by)) "combined_id" else color_by
  }
  validate_inputs(data, genes, group_by, treatment_samples, control_samples, color_by, normalisation)
  # Fetch and transform count matrix
  count_matrix <-
    switch(
      normalisation,
      raw = data@assays$RNA$counts,
      logNorm = compute_normalised_counts(data, method = "logNorm", batch = batch),
      cpm = compute_normalised_counts(data, method = "cpm", batch = batch),
      clr = compute_normalised_counts(data, method = "clr", batch = batch),
      SCT = compute_normalised_counts(data, method = "SCT", batch = batch),
      DESeq2 = compute_normalised_counts(data, method = "DESeq2", batch = batch),
      edgeR = compute_normalised_counts(data, method = "edgeR", batch = batch),
      limma_voom = compute_normalised_counts(data, method = "limma_voom", batch = batch),
      RUVg = compute_normalised_counts(data, method = "RUVg", batch = batch, spikes = spikes),
      RUVs = compute_normalised_counts(data, method = "RUVs", batch = batch),
      RUVr = compute_normalised_counts(data, method = "RUVr", batch = batch),
      zinb = compute_normalised_counts(data, method = "zinb", batch = batch),
      stop("Unsupported normalization method.")
    )
  
  #calculate expression on full samples
  combined_id_barcodes <- data@meta.data[[group_by]]
  names(combined_id_barcodes) <- rownames(data@meta.data)
  colnames(count_matrix) <- combined_id_barcodes[colnames(count_matrix)]
  treatment_control <- c(treatment_samples, control_samples)
  sub_count_matrix <- count_matrix[genes, colnames(count_matrix) %in% treatment_control, drop = FALSE]
  col_labels <- colnames(sub_count_matrix)
  #colnames(sub_count_matrix) <- make.names(col_labels, unique = TRUE)  # makes them unique (e.g. DMSO_0, DMSO_0.1, ...)
  colnames(sub_count_matrix) <- make.unique(col_labels, sep = "_repXXX")
  dt <- data.table::as.data.table(sub_count_matrix, keep.rownames = "Genes")
  dt_long <- data.table::melt(dt, id.vars = "Genes", variable.name = "ReplicateID", value.name = "Expression")
  data.table::setDT(dt_long)
  
  # make sure just remove numbers after "." at end of string that were for replicate numbers
  data.table::set(dt_long, j = "Treatment", value = gsub("_repXXX\\d*$", "", dt_long$ReplicateID)) 
  # dt_long[, Replicate := sequence(.N), by = .(Genes, Treatment)] 
  # convert the above code using data.table function to dplyr
  dt_long <- dt_long %>%
    dplyr::group_by(Genes, Treatment) %>%
    dplyr::mutate(Replicate = dplyr::row_number()) %>%
    dplyr::ungroup()
  data.table::setDT(dt_long)
  
  data.table::setcolorder(dt_long, c("Genes", "Treatment", "Replicate", "Expression"))
  n_samples <- length(colnames(sub_count_matrix))
  p <- ggplot(dt_long, aes(x = Treatment, y = Expression, group = Treatment)) +
    geom_boxplot(aes(fill = Treatment)) + facet_wrap(~Genes, ncol = 3, scales = "free_y") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_manual(values = c(macpie_colours$high, macpie_colours$low, macpie_colours$discrete[1:n_samples])) +
    labs(y = "Gene Expression") +
    macpie_theme()
  return(p)
}
