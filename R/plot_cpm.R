utils::globalVariables(c(".data", "Treatments", "CPM"))
#' Generate a box plot to show gene expression (CPM)
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
#' @import ggplot2
#' @import dplyr
#'
#' @returns a ggplot2 object
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' treatment_samples <- "Staurosporine_0.1"
#' control_samples <- "DMSO_0"
#' group_by <- "combined_id"
#' top_table_edgeR <- compute_single_de(mac, treatment_samples, control_samples, method = "edgeR")
#' genes <- c("TBRG4", "MRPL52", "DCTPP1", "ZFP36L1", "LSM1", "POLR2G")
#' plot_cpm(mac,genes, group_by, treatment_samples, control_samples)

plot_cpm <- function(data = NULL,
                      genes = NULL,
                      group_by = NULL,
                      treatment_samples = NULL,
                      control_samples = NULL,
                      color_by = NULL) {
  
  # Helper function to validate input data
  validate_inputs <- function(data, genes, group_by, treatment_samples, control_samples, color_by) {
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
    color_by <- if (is.null(color_by)) "combined_id" else color_by
  }
  validate_inputs(data, genes, group_by, treatment_samples, control_samples, color_by)
  #calculate cpm on full samples
  lcpm <- cpm(data@assays$RNA$counts, log = FALSE)
  combined_id_barcodes <- data@meta.data[[group_by]]
  names(combined_id_barcodes) <- rownames(data@meta.data)
  colnames(lcpm) <- combined_id_barcodes[colnames(lcpm)]
  treatment_control <- c(treatment_samples, control_samples)
  sub_lcpm <- lcpm[genes, colnames(lcpm) %in% treatment_control]
  sub_lcpm_m <- melt(sub_lcpm)
  colnames(sub_lcpm_m) <- c("Genes", "Treatments", "CPM")
  
  # Add color_by column
  sub_lcpm_m$color_by <- data@meta.data[match(sub_lcpm_m$Treatments, combined_id_barcodes), color_by]
  n_samples <- length(treatment_control)
  
  #reorder factors
  sub_lcpm_m$Treatments <- factor(sub_lcpm_m$Treatments, levels = c(treatment_samples, control_samples))

  p <- ggplot(sub_lcpm_m, aes(x = Treatments, y = CPM, group = Treatments)) +
    geom_boxplot(aes(fill = Treatments)) + facet_wrap(~Genes, ncol = 3) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    scale_fill_manual(values = macpie_colours$discrete[1:n_samples]) +
    labs(y = "Gene Expression (CPM)") +
    macpie_theme()

  return(p)
}
