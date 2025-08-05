utils::globalVariables(c("mac"))
#' Calculate QC metrics
#'
#' To calculate QC metrics such as standard deviation (sd), median absolute deviation (MAD),
#' interquartile range (IQR), and Z score
#' for read counts per condition/group of interest

#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param group_by A metadata column name to group data
#' @param order_by A column name or "median" for median of read counts.
#'
#' @return a list with a data frame with QC metrics and a box plot showing read counts per condition
#'
#' @import Seurat tidyseurat
#' @examples
#' data(mini_mac)
#' compute_qc_metrics(data = mini_mac, group_by = "combined_id", order_by = "median")
#' @export
compute_qc_metrics <- function(data = NULL, group_by = NULL, order_by = NULL) {
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop(
      "compute_qc_metrics(): the following package is required but not installed: forcats",
      "\nPlease install via `install.packages()`.")
  }
  # Helper function to validate input data
  validate_inputs <- function(data, group_by, order_by) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by

    if (!is.null(group_by) && !any(c(group_by) %in% c(colnames(data@meta.data), "median"))) {
      stop("Your parameter group_by is not a metadata column or 'median'.")
    }
    list(data = data, group_by = group_by, order_by = order_by)
  }

  validated <- validate_inputs(data, group_by, order_by)
  group_by <- validated$group_by
  meta <- validated$data@meta.data

  # Group by the specified column and assign replicate numbers
  meta <- meta %>%
    group_by(across(all_of(validated$group_by))) %>%
    mutate(Replicate = row_number())

  #global median and MAD
  global_median <- median(meta$nCount_RNA)
  global_mad <- mad(meta$nCount_RNA)

  stats_summary <- meta %>%
    group_by(across(all_of(validated$group_by))) %>%
    summarise(sd_value = round(sd(.data$nCount_RNA, na.rm = TRUE), 3),
              mad_value = round(mad(.data$nCount_RNA, na.rm = TRUE), 3),
              group_median = round(median(.data$nCount_RNA), 3),
              z_score = round((.data$group_median - global_median) / global_mad, 3),
              IQR = round(stats::IQR(.data$nCount_RNA), 3),
              .groups = "drop")# Calculate qc metrics per group
  if (order_by == "median") {
    meta[[group_by]] <- forcats::fct_reorder(meta[[group_by]], meta$nCount_RNA, median)
  } else if (order_by %in% colnames(stats_summary)) {
    # Order by a group-level summary metric
    ordering <- stats_summary %>%
      arrange(.data[[order_by]]) %>%
      pull(!!sym(group_by))
    meta[[group_by]] <- factor(meta[[group_by]], levels = ordering)
  } else if (order_by %in% colnames(meta)) {
    # Order by a per-cell metadata column (like Row, Column, Time, etc.)
    meta[[group_by]] <- forcats::fct_reorder(meta[[group_by]], meta[[order_by]], median)
  } else {
    stop(glue::glue("`order_by` = '{order_by}' not found in stats_summary or metadata"))
  }
  p <- ggplot(meta, aes(.data[[group_by]], .data$nCount_RNA)) +
    geom_boxplot(outlier.colour = "NA") +
    geom_jitter(width = 0.1, aes(color = .data[[group_by]])) +
    theme_minimal() +
    labs(x = group_by, y = "UMI counts") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position = "none")
  p
  return(list(stats_summary = stats_summary, plot = p))
}
