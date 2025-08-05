#' Create a lollipop chart
#'
#' Plot QC metrics such as standard deviation (sd), median absolute deviation (MAD),
#' interquartile range (IQR), or Z score
#' or read counts per condition/group of interest in a lollipop chart

#' @param stats_summary stats_summary a list
#' @param group_by A metadata column name to group data
#' @param metrics to specify one of the QC metrics from stats_summary calculated from QC_metrics
#'
#' @import ggplot2
#' @importFrom stats reorder mad sd
#' @returns A ggplot object with a lollipop chart of the specified QC metrics.
#' @examples
#' data("mini_mac")
#' qc_stats <- compute_qc_metrics(mini_mac, group_by = "combined_id", order_by = "median")
#' plot_qc_metrics(qc_stats, group_by = "combined_id", metrics = "sd_value")
#' @export


plot_qc_metrics <- function(stats_summary, group_by, metrics) {
  # Helper function to validate input data
  validate_inputs <- function(stats_summary, group_by, metrics) {
    if (!inherits(stats_summary, "list")) {
      stop("argument 'stats_summary' must be calculated QC metrics by running QC_metrics function.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by

    column_names <- stats_summary$stats_summary %>%
      head() %>%
      colnames()
    if (!all(c(metrics) %in% column_names)) {
      stop("Your metrics are not present in the stats_summary.")
    }
    list(stats_summary = stats_summary$stats_summary, group_by = group_by, metrics = metrics)

  }

  validate_inputs <- validate_inputs(stats_summary, group_by, metrics)

  var_stats <- validate_inputs$stats_summary
  group_by <- validate_inputs$group_by
  metrics <- validate_inputs$metrics

  p <- ggplot(var_stats, aes(x = reorder(.data[[group_by]], .data[[metrics]]), y = .data[[metrics]])) +
    geom_segment(aes(xend = .data[[group_by]], y = 0, yend = .data[[metrics]]),
                 color = "grey", size = 1) +
    geom_point(color = macpie_colours$discrete[[1]], size = 2) +
    coord_flip() +
    labs(x = group_by, y = metrics,
         title = '') +
    macpie_theme(show_x_title = TRUE, show_y_title = FALSE, legend_position_ = 'bottom', x_labels_angle = 0)
  p
}
