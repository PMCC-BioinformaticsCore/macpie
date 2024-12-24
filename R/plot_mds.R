#' Plot MDS dimensionality reduction
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns "Well_ID", "Row", "Column"
#' @param annotation A string specifying which column in data will be used to label a sample.
#' @importFrom limma plotMDS
#' @import edgeR
#' @import ggrepel
#' @return ggplot object
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_mds(mac, "Treatment_1")
plot_mds <- function(data, annotation) {

  dge <- edgeR::DGEList(counts = data@assays$RNA$counts)
  mds_result <- limma::plotMDS(dge, top = 500, plot = FALSE) # Extract MDS coordinates without plotting
  data$mds_dim1 <- mds_result$x
  data$mds_dim2 <- mds_result$y

  ggplot(data, aes(x = .data$mds_dim1,
                   y = .data$mds_dim2,
                   color = .data$Sample_type,
                   label = !!rlang::sym(annotation))) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = .data$Treatment_1),                    # Smart label repulsion
                    size = 3.5,
                    max.overlaps = 30) +        # Add sample labels
    theme_minimal() +                             # Minimal theme
    labs(
      title = "MDS Plot Using edgeR",
      x = "Dimension 1",
      y = "Dimension 2"
    )
}
