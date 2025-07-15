
#' Plot UMAP dimensionality reduction on DE genes
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param color_by A string specifying which column in data will be used to
#'   color the samples.
#' @param label A string specifying which column in data will be used to
#'   label a sample.
#' @param max_overlaps Maximum number of overlaps for ggrepel
#' @importFrom ggiraph geom_point_interactive girafe
#' @import ggplot2
#' @returns ggplot object
#' @export
#'
# #' @examples
plot_de_umap <- function(data = NULL, color_by = NULL, label = NULL, max_overlaps = NULL) {
  req_pkgs <- c("umap", "ggrepel")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "plot_de_umap(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }
  # Helper function to validate input data
  validate_inputs <- function(data, color_by, label, max_overlaps) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    color_by <- if (is.null(color_by)) "seurat_clusters" else color_by
    label <- if (is.null(label)) "combined_id" else label
    max_overlaps <- if (is.null(max_overlaps)) 20 else max_overlaps
    column_names <- data %>%
      head() %>%
      colnames()
    if (!all(c(color_by, label) %in% column_names)) {
      stop("Your color or label arguments are not present in the data. Try running FindNeighbors.")
    }
    list(data = data, color_by = color_by, label = label, max_overlaps = max_overlaps)
  }

  # Validate inputs
  validated <- validate_inputs(data, color_by, label, max_overlaps)
  color_by <- validated$color_by
  label <- validated$label
  max_overlaps <- validated$max_overlaps
  data <- validated$data

  cell_coords <- Embeddings(data, reduction = "umap_de") %>%
    as.data.frame() %>%
    rownames_to_column("combined_id") %>%
    left_join(data@meta.data, join_by("combined_id"))
  
  is_continuous <- is.numeric(cell_coords[[color_by]])
  
  # Plot with clusters and labels
  p <- ggplot(cell_coords, aes(x = UMAPde_1, 
                          y = UMAPde_2, 
                          color = !!rlang::sym(color_by),
                          label = !!rlang::sym(label))) +
    geom_point_interactive(aes(x = .data$UMAPde_1,
                               y = .data$UMAPde_2,
                               color = !!rlang::sym(color_by),
                               tooltip = !!rlang::sym(label),
                               data_id = !!rlang::sym(label))) +
    ggrepel::geom_text_repel(aes(label = !!rlang::sym(label)),                    
                    size = 3.5,
                    max.overlaps = max_overlaps) +        
    theme_minimal() +                             
    labs(
      title = "UMAP plot",
      x = "Dimension 1",
      y = "Dimension 2"
    ) +
    theme_minimal()
  
  if (is_continuous) {
    p <- p + scale_color_gradientn(colors = rev(macpie_colours$divergent))
  } else {
    p <- p + scale_color_manual(values = macpie_colours$discrete)
  }
  p
}
