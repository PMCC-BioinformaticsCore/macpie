utils::globalVariables(c("p_value_adj", "log2FC", "Significantly_upregulated", "Significantly_downregulated"))
#' Generate a table to summarise gene numbers from a differential
#' expression test.
#'
#' @param top_table A data table showing results from compute_single_de
#' @param lfc_threshold  Threshold of log2FC
#' @param padj_threshold Threshold of adjusted p value
#' @param multi To indicate to summarise for single de comparison or multi de comparison
#' @param group_by Name of the column that defines groups of replicates
#' @import dplyr
#'
#' @returns a data table
#' @export
#'
#' @examples
#' data("mini_mac")
#' top_table <- mini_mac@tools$diff_exprs$Staurosporine_10
#' summarise_de(top_table, lfc_threshold = 1, padj_threshold = 0.05)
#' summarise_de(mini_mac, lfc_threshold = 1, 
#' padj_threshold = 0.01, multi=TRUE)

summarise_de <- function(top_table,
                         lfc_threshold = 1,
                         padj_threshold = 0.01,
                         multi = FALSE, 
                         group_by = "combined_id") {
  if (!is.null(multi) && !multi %in% c(FALSE, TRUE)) {
    stop("Error: arugment 'multi' must be FALSE for single DE and TRUE for multi DE.")
  }
  validate_inputs_single <- function(top_table,
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
  validate_inputs_multi <- function(top_table,
                                    lfc_threshold,
                                    padj_threshold) {
    if (!inherits(top_table, "Seurat")) {
      stop("Error: argument 'top_table' must be a Seurat or TidySeurat object.")
    }
    if (!inherits(lfc_threshold, "numeric")) {
      stop("Error: arguemnt 'lfc_threshold' must be numeric.")
    }
    if (!inherits(padj_threshold, "numeric")) {
      stop("Error: arguemnt 'padj_threshod' must be numeric.")
    }
  }
  if (multi==FALSE) {
    validate_inputs_single(top_table, lfc_threshold, padj_threshold)
    total_genes <- nrow(top_table)
    up <- top_table %>% filter(log2FC >= lfc_threshold & p_value_adj <= padj_threshold) %>% nrow()
    down <-  top_table %>% filter(log2FC < -lfc_threshold & p_value_adj <= padj_threshold) %>% nrow()
    n_sig <- up + down
    summary_tbl <- tibble(
      Total_genes_tested = total_genes,
      Significantly_upregulated = up,
      Significantly_downregulated = down,
      Total_significant = n_sig,
      Padj_threshold = padj_threshold,
      Log2FC_threshold = lfc_threshold
    )
  } else {
    validate_inputs_multi(top_table, lfc_threshold, padj_threshold)
    all_de <- top_table@tools$diff_exprs
    de_df <- bind_rows(all_de)
    summary_tbl <- de_df %>% group_by(.data[[group_by]]) %>%
      summarise(
        Total_genes_tested = n(),
        Significantly_upregulated = sum(log2FC >= lfc_threshold & p_value_adj <= padj_threshold, na.rm = TRUE),
        Significantly_downregulated = sum(log2FC < -lfc_threshold & p_value_adj <= padj_threshold, na.rm = TRUE),
        Total_significant = Significantly_upregulated + Significantly_downregulated,
        padj_threshold = padj_threshold,
        Log2FC_threshold = lfc_threshold,
        .groups = "drop"
      )
  }
  return(summary_tbl)
}
