#' Plot MDS dimensionality reduction
#'
#' This function uses limma's plot MDS to visualise groupoing of the data
#' points. The function is used with gene.selection = "common" to increase the
#' speed.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param group_by A string specifying which column in data will be used to
#'   color the samples.
#' @param label A string specifying which column in data will be used to
#'   label a sample.
#' @importFrom limma plotMDS
#' @importFrom utils head
#' @importFrom ggsci scale_color_npg
#' @importFrom ggiraph geom_point_interactive girafe
#' @import edgeR
#' @import ggrepel

#' @return ggplot object
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_mds(mac)
#' @export
plot_mds <- function(data = NULL, group_by = NULL, label = NULL, max_overlaps = NULL) {

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

  dge <- edgeR::DGEList(counts = data@assays$RNA$counts)
  # Extract MDS coordinates without plotting
  #use gene.selection = "commmon" to increase speed
  mds_result <- limma::plotMDS(dge,
                               gene.selection = "common",
                               top = 500, plot = FALSE)

  data$PCA1 <- mds_result$x
  data$PCA2 <- mds_result$y
  # Plot the results
  tryCatch({
    p <- ggplot(data, aes(x = .data$PCA1,
                          y = .data$PCA2,
                          color = .data$Sample_type,
                          label = !!rlang::sym(label))) +
      #geom_point(size = 2) +
      geom_point_interactive(aes(x = .data$PCA1,
                                 y = .data$PCA2,
                                 color = .data$Sample_type,
                                 tooltip = !!rlang::sym(label),
                                 data_id = !!rlang::sym(label))) +
      geom_text_repel(aes(label = .data$combined_id),                    # Smart label repulsion
                      size = 3.5,
                      max.overlaps = max_overlaps) +        # Add sample labels
      theme_minimal() +                             # Minimal theme
      scale_color_npg() +
      labs(
        title = "MDS plot",
        x = "Dimension 1",
        y = "Dimension 2"
      ) +
      theme_minimal()
    p
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
