utils::globalVariables(c("Metric", "Value", "Normalized"))
#' Create a heatmap for multiple QC metrics
#'
#' Plots multiple QC metrics from stats_summary as a heatmap, normalizing values between 0 and 1.
#'
#' @param stats_summary A list containing QC metrics.
#' @param group_by A metadata column name to group data.
#' @param metrics A vector of QC metrics to visualize.
#' @param order_by The metric to use for sorting group_by values.
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @importFrom scales rescale
#' @examples
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' qc_stats <- compute_qc_metrics(mac, group_by = "combined_id", order_by = "median")
#' plot_qc_metrics_heatmap(stats_summary = qc_stats$stats_summary, group_by = "combined_id")
#' 
#' @export

plot_qc_metrics_heatmap <- function(stats_summary = NULL, 
                                    group_by = NULL, 
                                    metrics = NULL, 
                                    order_by = NULL) {
  # Helper function to validate input data
  validate_inputs <- function(stats_summary, group_by, metrics, order_by) {
    if (!inherits(stats_summary, "tbl_df")) {
      stop("Error: argument 'stats_summary' must be calculated QC metrics by running QC_metrics function.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    metrics <- if (is.null(metrics)) c("sd_value","z_score","mad_value","IQR") else metrics
    order_by <- if (is.null(order_by)) "sd_value" else order_by
    
    if (!all(metrics %in% colnames(stats_summary))) {
      stop("Error: One or more specified metrics are not present in stats_summary.")
    }
    
    if (!(order_by %in% metrics)) {
      stop("Error: 'order_by' metric must be one of the selected metrics.")
    }
    
    list(stats_summary = stats_summary, group_by = group_by, metrics = metrics, order_by = order_by)
  }
  
  validated_inputs <- validate_inputs(stats_summary, group_by, metrics, order_by)
  var_stats <- validated_inputs$stats_summary
  group_by <- validated_inputs$group_by
  metrics <- validated_inputs$metrics
  order_by <- validated_inputs$order_by
  # abs z-score
  var_stats$z_score <- abs(var_stats$z_score)
  # Reshape data to long format
  heatmap_data <- var_stats %>%
    select(all_of(c(group_by, metrics))) %>%
    pivot_longer(cols = all_of(metrics), names_to = "Metric", values_to = "Value") %>%
    group_by(Metric) %>%
    mutate(Normalized = scales::rescale(Value, to = c(0, 1))) %>%  # Normalize to 0-1
    ungroup()
  
  # Compute sorting order based on `order_by` metric
  sorting_order <- heatmap_data %>%
    filter(Metric == order_by) %>%
    arrange(desc(Value)) %>%  # Sort descending
    pull(.data[[group_by]])
  
  heatmap_data$Metric <- factor(heatmap_data$Metric, levels = metrics)
  
  # Plot heatmap
  p <- ggplot(heatmap_data, aes(x = Metric, y = factor(.data[[group_by]], levels = sorting_order), fill = Normalized)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = macpie_colours$discrete[7]) +  # Change color scale if needed
    labs(x = "QC Metric", y = group_by, fill = "Normalized Value") +
    theme_minimal() +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4))
  
  return(p)
}
