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
#' @param n_labels An integer specifying number of labels to show, based on PC1 and PC2 extremeness. Default set to 50
#' @param max_overlaps Maximum number of overlaps for ggrepel
#' @importFrom limma plotMDS
#' @importFrom utils head
#' @importFrom ggsci scale_color_npg
#' @importFrom ggiraph geom_point_interactive girafe
#' @import edgeR
#' @import ggrepel
#' @return ggplot object
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' p <- plot_mds(mac)
#' @export
#'
plot_mds <- function(data = NULL, group_by = NULL, label = NULL, max_overlaps = NULL, n_labels = 50) {

  # Helper function to validate input data
  validate_inputs <- function(data, group_by, label, max_overlaps) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "Sample_type" else group_by
    label <- if (is.null(label)) "combined_id" else label
    max_overlaps <- if (is.null(max_overlaps)) 10 else max_overlaps
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

  data <- data %>%
    mutate(extremeness = abs(.data$PCA1) + abs(.data$PCA2))

  data$group <- data %>%
    select(!!rlang::sym(group_by)) %>%
    pull()
  
  data$label <- data %>%
    select(!!rlang::sym(label)) %>%
    pull()
  
  # Select the top 50 most extreme points
  top_n_data <- data@meta.data %>%  # Sort by extremeness in descending order
    arrange(desc(.data$extremeness)) %>% head(n_labels)

  # Plot the results
  tryCatch({
    p <- ggplot(data, aes(x = .data$PCA1, y = .data$PCA2, color = .data$group, label = !!rlang::sym(label))) +
      #geom_point(size = 2) +
      geom_point_interactive(aes(x = .data$PCA1, y = .data$PCA2,
                                 color = .data$group, 
                                 tooltip = !!rlang::sym(label), 
                                 data_id = !!rlang::sym(label))) +
      geom_text_repel(data = top_n_data, aes(label = .data$label), size = 3.5, max.overlaps = max_overlaps, show.legend = F) + # Add sample labels
      scale_color_manual(values = macpie_colours$discrete) +
      labs(title = "MDS plot",
           x = "Dimension 1",
           y = "Dimension 2") +
      macpie_theme()
    p
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
