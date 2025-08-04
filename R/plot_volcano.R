#' Volcano plot of differentially expressed genes
#'
#' @param top_table Data frame with columns that contain FC and P-values, default to log2FC and p_value_adj
#' @param x Name of the column with logFC values
#' @param y Name of the column with adjusted p-values
#' @param fdr_cutoff Cutoff for labels to be plotted based on their adjusted p-values
#' @param max.overlaps Maximum number of overlaps of points for package ggrepel
#' @import stringr
#' @importFrom stats setNames
#' @returns A ggplot object
#' @export
#'
#' @examples
#' data("mini_mac")
#' top_table <- mini_mac@tools$diff_exprs$Staurosporine_10
#' plot_volcano(top_table = top_table, max.overlaps = 100)
#'
plot_volcano <- function(top_table, x = "log2FC", y = "p_value_adj", fdr_cutoff = 0.05, max.overlaps = 30) {

  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop(
      "plot_volcano(): the following package is required but not installed: forcats",
      "\nPlease install via `install.packages()`.")
  }
  
  # Helper function to validate input data
  validate_inputs <- function(top_table, x, y, fdr_cutoff, max.overlaps) {
    if (!inherits(as.data.frame(top_table), "data.frame")) {
      stop("Error: 'top_table' must be convertable to a data frame.")
    }
    if (!all(x %in% colnames(top_table) && y %in% colnames(top_table))) {
      stop("Error: 'top_table' does not contain expected column names, check paramaterrs 'x' and 'y'.")
    }
    if (!inherits(fdr_cutoff, "numeric") && fdr_cutoff <= 1) {
      stop("FDR should be numeric and smaller than 1.")
    }
    if (!inherits(max.overlaps, "numeric")) {
      stop("max.overlaps should be numeric.")
    }
  }

  validate_inputs(top_table, x, y, fdr_cutoff, max.overlaps)

  top_table$diff_expressed <- "no"
  top_table$gene_labels <- ""

  top_table <- top_table %>%
    mutate(diff_expressed = case_when(!!rlang::sym(x) >= 1 & !!rlang::sym(y) <= fdr_cutoff ~ "up",
                                      !!rlang::sym(x) <= -1 & !!rlang::sym(y) <= fdr_cutoff ~ "down",
                                      TRUE ~ diff_expressed),
           diff_expressed = forcats::fct_expand(diff_expressed, "up", "down", "no"), # Add all expected levels
           gene_labels = case_when(.data$diff_expressed != "no" ~ .data$gene, TRUE ~ ""),
           colors = case_when(
             .data$diff_expressed == "up"   ~ macpie_colours$high,
             .data$diff_expressed == "down" ~ macpie_colours$low,
             TRUE                           ~ "gray"
           ))

  color_mapping <- unique(top_table[, c("diff_expressed", "colors")])
  named_colors <- setNames(color_mapping$colors, color_mapping$diff_expressed)

  suppressWarnings({
    p <- ggplot(top_table, aes(x = .data$log2FC, y = -log10(.data$p_value_adj),
                               group = .data$diff_expressed, col = .data$diff_expressed,
                               label = gene_labels)) +
      scale_y_continuous(name = expression(-log[10](p-value[adj]))) +
      geom_point() +
      ggrepel::geom_text_repel(min.segment.length = 5, max.overlaps = max.overlaps, show.legend = FALSE) +
      scale_color_manual(values = named_colors) + # Map colors dynamically
      # Boundaries
      geom_vline(xintercept = c(-1, 1), col = "#003366", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "#003366", linetype = 'dashed') +
      # Theme
      macpie_theme()
  })

  return(p)
}
