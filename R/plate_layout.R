#' Plot MAC-seq data on a plate layout
#' @importFrom stats median
#' @importFrom gtools mixedsort
#' @importFrom forcats fct_reorder
#' @importFrom rlang sym
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @import tidyseurat
#' @import rlang
#' @import ggplot2
#' @param data A tidyseurat object merged with metadata. Must contain columns "Well_ID", "Row", "Column"
#' @param metric A string specifying which column in data will be used to color a sample.
#' @param annotation A string specifying which column in data will be used to label a sample.
#' @return ggplot object
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plate_layout(mac,"nCount_RNA","Treatment_1")

plate_layout <- function(data, metric, annotation) {
  data <- data %>%
    dplyr::select("Well_ID", "Row", "Column", {{metric}}, {{annotation}}) %>%
    mutate(Col = as.character(.data$Column)) %>%
    mutate(Col = factor(.data$Col, levels = gtools::mixedsort(unique(.data$Col)))) %>%
    mutate(median_value = median(!!rlang::sym(metric))) %>%
    mutate(max_value = max(!!rlang::sym(metric))) %>%
    mutate(min_value = min(!!rlang::sym(metric)))
  p <- ggplot(data, aes(.data$Col,
    forcats::fct_rev(forcats::as_factor(.data$Row)),
    fill = !!rlang::sym(metric),
  )
  ) +
    geom_tile() +
    scale_x_discrete(position = "top") +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = unique(data$median_value),
                         name = {{metric}}) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) +
    geom_text(aes(label = !!rlang::sym(annotation)), size = 2)
  p
}
