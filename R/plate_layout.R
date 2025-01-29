#' Plot MAC-seq data on a plate layout
#' @importFrom stats median
#' @importFrom gtools mixedsort
#' @importFrom forcats fct_reorder
#' @importFrom rlang sym
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom utils head
#' @import tidyseurat
#' @import rlang
#' @import ggplot2
#' @importFrom dplyr select
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param metric A string specifying which column in data will be used to color
#'   a sample. Defaults to "nCount_RNA".
#' @param annotation A string specifying which column in data will be used to
#'   label a sample. Defaults to "Treatment_1".
#' @return ggplot object
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plate_layout(mac,"nCount_RNA","Treatment_1")

plate_layout <- function(data = NULL, metric = NULL, annotation = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, metric, annotation) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    metric <- if (is.null(metric)) "nCount_RNA" else metric
    annotation <- if (is.null(annotation)) "Treatment_1" else annotation

    suppressWarnings({column_names <- data %>%
      head() %>%
      colnames()
    })
    if (!all(c(metric, annotation) %in% column_names)) {
      stop("Your column names are not present in the data or metadata.")
    }
    list(data = data, metric = metric, annotation = annotation)
  }

  # Validate inputs
  validated <- validate_inputs(data, metric, annotation)
  metric <- validated$metric
  annotation <- validated$annotation
  data <- validated$data

  # Find stats for the data to arrange plotting
  suppressWarnings({
    data <- data %>%
      mutate(Col = as.character(.data$Column)) %>%
      mutate(Col = factor(.data$Col, levels = gtools::mixedsort(unique(.data$Col)))) %>%
      mutate(Row = factor(.data$Row, levels = gtools::mixedsort(unique(.data$Row)))) %>%
      mutate(median_value = median(!!rlang::sym(metric))) %>%
      mutate(max_value = max(!!rlang::sym(metric))) %>%
      mutate(min_value = min(!!rlang::sym(metric)))
  })
  # Plot the results

  tryCatch({
    p <- ggplot(data, aes(.data$Col,
                          forcats::fct_rev(forcats::as_factor(.data$Row)),
                          fill = !!rlang::sym(metric))) +
      geom_tile() +
      scale_x_discrete(position = "top") +
      scale_fill_gradient2(low = "#3C5488FF",
                           mid = "white",
                           high = "#E64B35FF",
                           midpoint = unique(data$median_value),
                           name = {{metric}}) +
      theme(
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank()
      ) +
      geom_text(aes(label = !!rlang::sym(annotation)), size = 1.5)
    p
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })

}
