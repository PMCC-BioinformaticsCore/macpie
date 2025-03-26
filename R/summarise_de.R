utils::globalVariables(c("p_value_adj","log2FC"))
#' Generate a table to summarise gene numbers from a differential
#' expression test. 
#'
#' @param top_table A data table showing results from compute_single_de
#' @param lfc_threshold  Threshold of log2FC
#' @param padj_threshold threshold of adjusted p value
#'
#' @import dplyr
#'
#' @returns a data table
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' treatment_samples <- "Staurosporine_0.1"
#' control_samples <- "DMSO_0"
#' group_by <- "combined_id"
#' top_table_edgeR <- compute_single_de(mac, treatment_samples, control_samples, method = "edgeR")
#' summarise_de(top_table_edgeR, lfc_threshold = 1, padj_threshold = 0.05)

summarise_de <- function(top_table,
                         lfc_threshold = 1,
                         padj_threshold = 0.01) {
  validate_inputs <- function(top_table, 
                              lfc_threshold,
                              padj_threshold) {
    if (!inherits(top_table, "data.frame")) {
      stop("Error: argument 'top_table' must be a data frame.")
    }
    if (!inherits(lfc_threshold, "numeric")) {
      stop("Error: arguemnt 'lfc_threshold' must be numeric.")
    }
    if (!inherits(padj_threshold, "numeric")) {
      stop("Error: arguemnt 'padj_threshod' must be numeric.")
    }
  }
  validate_inputs(top_table, lfc_threshold, padj_threshold)
  total_genes <- nrow(top_table)
  up <- top_table %>% filter(log2FC >= lfc_threshold & p_value_adj <= padj_threshold) %>% nrow()
  down <-  top_table %>% filter(log2FC <= -lfc_threshold & p_value_adj <= padj_threshold) %>% nrow()
  n_sig <- up + down
  summary_tbl <- tibble(
    Total_genes_tested = total_genes,
    Significantly_upregulated = up,
    Significantly_downregulated = down,
    Total_significant = n_sig,
    Padj_threshold = padj_threshold,
    Log2FC_threshold = lfc_threshold
  )
return(summary_tbl)
}