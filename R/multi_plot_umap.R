
#' Plot UMAP dimensionality reduction on DE genes
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param group_by A string specifying which column in data will be used to
#'   color the samples.
#' @param label A string specifying which column in data will be used to
#'   label a sample.
#' @param max_overlaps Maximum number of overlaps for ggrepel
#' @importFrom umap umap
#' @importFrom ggiraph geom_point_interactive girafe
#' @import ggrepel
#' @import ggplot2
#' @returns ggplot object
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' mac <- multi_prepare_umap(mac)
#' multi_plot_umap(mac)
multi_plot_umap <- function(data = NULL, group_by = NULL, label = NULL, max_overlaps = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, group_by, label, max_overlaps) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "Sample_type" else group_by
    label <- if (is.null(label)) "combined_id" else label
    max_overlaps <- if (is.null(max_overlaps)) 20 else max_overlaps
    column_names <- data %>%
      head() %>%
      colnames()
    if (!all(c(group_by, label) %in% column_names)) {
      stop("Your column names are not present in the data or metadata.")
    }
    list(data = data, group_by = group_by, label = label, max_overlaps = max_overlaps)
  }

  # Validate inputs
  validated <- validate_inputs(data, group_by, label, max_overlaps)
  group_by <- validated$group_by
  label <- validated$label
  max_overlaps <- validated$max_overlaps
  data <- validated$data

  df_umap_data <- data@tools[["umap_de"]]
  df_umap_data <- df_umap_data %>%
    left_join(., data@meta.data, join_by("combined_id"))

  p <- ggplot(df_umap_data, aes(x = .data$UMAP_1,
                                y = .data$UMAP_2,
                                color = .data$Sample_type,
                                label = !!rlang::sym(label))) +
    geom_point_interactive(aes(x = .data$UMAP_1,
                               y = .data$UMAP_2,
                               color = .data$Sample_type,
                               tooltip = !!rlang::sym(label),
                               data_id = !!rlang::sym(label))) +
    geom_text_repel(aes(label = .data$combined_id),                    # Smart label repulsion
                    size = 3.5,
                    max.overlaps = max_overlaps) +        # Add sample labels
    theme_minimal() +                             # Minimal theme
    scale_color_npg() +
    labs(
      title = "UMAP plot",
      x = "Dimension 1",
      y = "Dimension 2"
    ) +
    theme_minimal()
  p
}
