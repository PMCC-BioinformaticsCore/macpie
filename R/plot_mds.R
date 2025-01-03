#' Plot MDS dimensionality reduction
#'
#' This function uses limma's plot MDS to visualise groupoing of the data
#' points. The function is used with gene.selection = "common" to increase the
#' speed.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param annotation A string specifying which column in data will be used to
#'   label a sample.
#' @importFrom limma plotMDS
#' @importFrom utils head
#' @import edgeR
#' @import ggrepel

#' @return ggplot object
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_mds(mac, "Treatment_1")
#' @export
plot_mds <- function(data = NULL, annotation = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, annotation) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    annotation <- if (is.null(annotation)) "Treatment_1" else annotation

    column_names <- data %>%
      head() %>%
      colnames()
    if (!all(c(annotation) %in% column_names)) {
      stop("Your column names are not present in the data or metadata.")
    }
    list(data = data, annotation = annotation)
  }

  # Validate inputs
  validated <- validate_inputs(data, annotation)
  annotation <- validated$annotation
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
                          label = !!rlang::sym(annotation))) +
      geom_point(size = 2) +
      geom_text_repel(aes(label = .data$Treatment_1),                    # Smart label repulsion
                      size = 3.5,
                      max.overlaps = 20) +        # Add sample labels
      theme_minimal() +                             # Minimal theme
      labs(
        title = "MDS plot",
        x = "Dimension 1",
        y = "Dimension 2"
      )
    p
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
