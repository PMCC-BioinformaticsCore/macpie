#' Create an RLE Plot
#'
#' This function generates a Relative Log Expression (RLE) plot for visualizing
#' the distribution of expression data after normalization or log transformation.
#'
#' @param count_matrix A numeric matrix or data frame containing the expression data.
#' Each column represents a sample, and each row represents a feature.
#' @param id A vector of sample identifiers corresponding to the columns of `count_matrix`.
#' @param feature A vector of feature labels corresponding to the columns of `count_matrix`.
#' @param logged A logical value indicating whether `count_matrix` is already log-transformed.
#' Defaults to `FALSE`.
#'
#' @details
#' The function performs the following steps:
#' - Ensures input data is numeric and matches dimensions.
#' - Optionally log-transforms the data if `logged = FALSE`.
#' - Computes the RLE by subtracting the row medians from each value.
#' - Creates a boxplot using ggplot2 to visualize the distribution of RLE values.
#'
#' @return A ggplot object representing the RLE plot.
#'
#' @examples
#' # Example Data
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' count_matrix<-as.matrix(mac@assays$RNA$counts)
#' colnames(count_matrix)<-mac$Well_ID
#' rle_plot(count_matrix = count_matrix, id = mac$Well_ID, feature = mac$Row, logged=FALSE)
#'
#'
#' @importFrom Biobase rowMedians
#' @importFrom grDevices boxplot.stats
#' @importFrom grDevices colorRampPalette
#' @import ggplot2 tidyseurat
#' @importFrom RColorBrewer brewer.pal
#' @export
#'

rle_plot <- function(count_matrix, id, feature, logged = FALSE) {
  # Check if count_matrix is a numeric matrix or data frame
  if (!is.matrix(count_matrix) && !is.data.frame(count_matrix)) {
    stop("Error: 'count_matrix' must be a matrix or data frame.")
  }

  # Ensure count_matrix is numeric
  if (!all(sapply(count_matrix, is.numeric))) {
    stop("Error: 'count_matrix' must contain only numeric values.")
  }

  # Convert to matrix if data frame
  count_matrix <- as.matrix(count_matrix)

  # Validate id
  if (length(id) != ncol(count_matrix)) {
    stop("Error: Length of 'id' must match the number of columns in 'count_matrix'.")
  }

  # Convert id to a factor
  id <- as.factor(id)

  # Validate feature
  if (length(feature) != ncol(count_matrix)) {
    stop("Error: Length of 'feature' must match the number of columns in 'count_matrix'.")
  }

  # Convert feature to a factor
  feature <- as.factor(feature)

  # Add 1 to count_matrix and apply log2 transformation
  # Check if count_matrix is logged
  if (!logged) {
    count_matrix <- log2(count_matrix + 1)
  }

  # Compute RLE
  rle <- count_matrix - Biobase::rowMedians(count_matrix)

  # Sort RLE based on feature
  sort_index <- sort.list(feature)
  rle <- rle[, sort_index]

  # Create a data frame for RLE boxplot stats
  rledf <- t(apply(rle, 2, function(x) {
    grDevices::boxplot.stats(x)$stats
  }))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")

  # Add feature and sample information
  rledf$feature <- feature[sort_index]
  rledf$sample <- colnames(rle)

  # Reorder samples for plotting
  rledf$sample <- factor(rledf$sample, levels = unique(rledf$sample))

  # Plot the RLE
  tryCatch({
    p <- ggplot(rledf, aes(x = .data$sample, fill = .data$feature)) +
      geom_boxplot(
        aes(
            ymin = .data$ymin,
            lower = .data$lower,
            middle = .data$middle,
            upper = .data$upper,
            ymax = .data$ymax),
        stat = "identity"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_x_discrete(limits = rledf$sample) +
      scale_color_manual(values = colorRampPalette(
                                                   brewer.pal(12, "Paired"))
      (nlevels(feature))) +
      geom_hline(yintercept = 0,
                 linetype = "dotted",
                 col = "red",
                 linewidth = 1) +
      ylab("log2_expression_deviation")

    return(p)
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
