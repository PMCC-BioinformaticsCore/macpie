#' Plot multi-screen profile from fgsea results
#'
#' @param data A tidyseurat object with `screen_profile` in the @tools slot.
#' @param color_by A string specifying which column in screen_profile to use for color (default: automatically chosen).
#' @param size_by A string specifying which column in screen_profile to use for point size (default: log10(padj)).
#' @param label_angle Angle of x-axis labels (default: 90).
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap scale_color_gradientn scale_color_manual theme_minimal labs
#' @importFrom dplyr mutate arrange
#' @returns A ggplot object
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10")
#' plot_multi_screen_profile(mini_mac, color_by = "Sample_type")
plot_multi_screen_profile <- function(data,
                                      color_by = NULL,
                                      size_by = "logPadj",
                                      label_angle = 90) {
  
  # Validate input
  if (is.null(data@tools$screen_profile)) {
    stop("The input object does not contain 'screen_profile' in the @tools slot.")
  }
  if (is.null(data@tools$screen_profile)) {
    stop("The input object does not contain 'screen_profile' in the @tools slot.")
  }
  color_by <- if (is.null(color_by)) "Sample_type" else color_by
  column_names <- data %>%
    head() %>%
    colnames()
  if (!all(c(color_by) %in% column_names)) {
    stop("Your color arguments is not present in the data.")
  }
  df <- data@tools$screen_profile
  
  # Choose default color_by
  metric <- if (!is.null(df$percent_target_activity)) {
    "percent_target_activity"
  } else {
    "NES"
  }
  
  # Add logPadj if not already present
  if (!"logPadj" %in% colnames(df)) {
    df <- df %>% mutate(logPadj = -log10(padj))
  }
  
  # Order target by effect size
  df <- df %>%
    arrange(desc(metric)) %>%
    mutate(target = forcats::fct_reorder(target_id, !!rlang::sym(metric),.desc = TRUE)) %>%
    left_join(data@meta.data, join_by(target_id == combined_id))
  
  # Continuous or discrete?
  is_continuous <- is.numeric(df[[color_by]])
  
  # Base plot
  p <- ggplot(df, aes(x = target, 
                      y = !!rlang::sym(metric), 
                      color = !!rlang::sym(color_by),
                      label = target)) +
    geom_point_interactive(aes(x = target,
                               y = !!rlang::sym(metric), 
                               color = !!rlang::sym(color_by),
                               tooltip = target,
                               data_id = target)) +
    geom_point(aes(size = .data[[size_by]])) +
    labs(
      x = "Target profile",
      y = metric,
      title = "Profile enrichment"
    ) +
    macpie_theme(x_labels_angle = label_angle, show_x_title = FALSE) +
    theme(
      panel.grid.major.x = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
      panel.grid.minor.y = element_blank()
    )
  
  if (is_continuous) {
    p <- p + scale_color_gradientn(colors = rev(macpie_colours$divergent))
  } else {
    p <- p + scale_color_manual(values = macpie_colours$discrete)
  }
  p
  return(p)
}

