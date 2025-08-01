utils::globalVariables(c("tooltip_text", "group_id"))
#' Plot MAC-seq data on a plate layout
#' @importFrom stats median
#' @importFrom dplyr mutate
#' @importFrom rlang .data sym
#' @importFrom utils head
#' @import tidyseurat
#' @import ggplot2 
#' @import ggiraph
#' @importFrom dplyr select
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param metric A string specifying which column in data will be used to color
#'   a sample. Defaults to "nCount_RNA".
#' @param annotation A string specifying which column in data will be used to
#'   label a sample. Defaults to "Treatment_1".
#' @param midpoint A value to be used in heatmap scale bar, it can be either mean value or median value.
#' @return ggplot object
#'
#' @examples
#' 
#' data("mini_mac")
#' p <- plot_plate_layout(mini_mac, metric = "nCount_RNA",
#' annotation = "Treatment_1", midpoint = "mean")


#' @export

plot_plate_layout <- function(data = NULL, metric = NULL, annotation = NULL, midpoint = NULL) {
  
  
  req_pkgs <- c("gtools", "forcats")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "plot_plate_layout(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }
  
  # Helper function to validate input data
  validate_inputs <- function(data, metric, annotation, midpoint) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    metric <- if (is.null(metric)) "nCount_RNA" else metric
    annotation <- if (is.null(annotation)) "Treatment_1" else annotation
    midpoint <- if (is.null(midpoint)) "mean" else midpoint
    if (!midpoint %in% c("mean","median")){
      stop("Error: argument 'midpoint' must be either mean or median.")
    }
    
    suppressWarnings({column_names <- data %>%
      head() %>%
      colnames()
    })
    if (!all(c(metric, annotation) %in% column_names)) {
      stop("Your column names are not present in the data or metadata.")
    }
    list(data = data, metric = metric, annotation = annotation, midpoint = midpoint)
  }
  
  # Validate inputs
  validated <- validate_inputs(data, metric, annotation, midpoint)
  # validated <- validate_inputs(mac_outliers, "nCount_RNA","combined_id", midpoint = "median")
  metric <- validated$metric
  annotation <- validated$annotation
  data <- validated$data
  midpoint <- validated$midpoint
  
  # Find stats for the data to arrange plotting
  suppressWarnings({
    data <- data %>%
      mutate(Col = as.character(.data$Column)) %>%
      mutate(Col = factor(.data$Col, levels = gtools::mixedsort(unique(.data$Col)))) %>%
      mutate(Row = factor(.data$Row, levels = gtools::mixedsort(unique(.data$Row)))) %>%
      mutate(median_value = median(!!rlang::sym(metric))) %>%
      mutate(mean_value = mean(!!rlang::sym(metric))) %>%
      mutate(
        tooltip_text = paste0("Sample: ", .data[[annotation]], "\nValue: ", .data[[metric]]),
        group_id = .data[[annotation]]
      )
  })
  # Plot the results
  tryCatch({
    if (midpoint == "median"){
      midpoint <- unique(data$median_value)
    } else {
      midpoint <- unique(data$mean_value)
    }
    p <- ggplot(data, aes(
      x = .data$Col,
      y = forcats::fct_rev(forcats::as_factor(.data$Row)),
      fill = !!rlang::sym(metric)
    )) +
      facet_wrap(~Plate_ID, ncol = 1) +
      geom_tile_interactive(aes(tooltip = tooltip_text,
                                data_id = group_id), color = "grey90") +
      scale_x_discrete(position = "bottom") +
      scale_fill_gradient2(
        high = macpie_colours$high,
        mid = "white",
        low = macpie_colours$low,
        midpoint = midpoint,
        name = rlang::as_name(metric)
      ) +
      ylab("Row") +
      macpie_theme(show_x_title = FALSE, show_y_title = FALSE, legend_position_ = 'right') +
      guides(fill = guide_colorbar(barwidth = 0.5, barheight = 10))
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
