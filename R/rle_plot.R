#' Create an RLE Plot
#'
#' This function generates a Relative Log Expression (RLE) plot for visualizing
#' the distribution of expression data after normalization or log transformation.
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns "Well_ID", "Row", "Column".
#' @param barcodes A vector of sample barcodes corresponding to Cells(seurat_object).
#' @param labels A vector of labels of the same length as 'barcodes" to group the barcodes.
#' @param label_column A metadata column name to group the barcodes.
#' @param log A logical value indicating whether data should be log-transformed.
#' Defaults to `TRUE`.
#'
#' @details
#' The function performs the following steps:
#' - Ensures integrity of input data
#' - Log-transforms the data
#' - Computes the RLE by subtracting the row medians from each value.
#' - Creates a boxplot using ggplot2 to visualize the distribution of RLE values.
#'
#' @return A ggplot object representing the RLE plot.
#'
#' @examples
#' # Example Data
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' rle_plot(data = mac, barcodes = Seurat::Cells(mac), label_column = "Row", log=TRUE)
#'
#'
#' @importFrom Biobase rowMedians
#' @importFrom grDevices boxplot.stats
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr pull select
#' @import Seurat ggplot2 tidyseurat
#' @importFrom RColorBrewer brewer.pal
#' @export
#'

rle_plot <- function(data, barcodes=NULL, labels=NULL, label_column=NULL, log = TRUE) {

  # Check if mac is of the right format
  if (!inherits(data,"Seurat")) {
    stop("Error: 'data' must be a Seurat or TidySeurat object.")
  }

  # If there are no barcodes, use the whole set
  if (is.null(barcodes)) {
    barcodes <- Seurat::Cells(data)
  }

  if (!is.null(label_column) && inherits(label_column,"character") && length(label_column) == 1) {
    labels <- data %>%
      select({{label_column}}) %>%
      pull()
  } else if(is.null(labels)) {
    stop("The format of `label_column` should be a single character value.")
  }

  if (is.null(labels)) {
    stop("Either `labels` or `label_column` must be provided.")
  }

  # Ensure alignment
  if (length(labels) != ncol(data)) {
    stop("Labels must have the same length as the number of columns in the dataset")
  }

  # Ensure log is by default TRUE
  if (inherits(log,"function")) {
    log = TRUE
  }

  #fetch the count_matrix
  count_matrix<-as.matrix(data@assays$RNA$counts)

  # Validate id
  if (length(barcodes) != ncol(count_matrix)) {
    stop("Error: Length of 'id' must match the number of columns in 'count_matrix'.")
  }

  # Convert id to a factor
  barcodes <- as.factor(barcodes)

  # Validate feature
  if (length(labels) != ncol(count_matrix)) {
    stop("Error: Length of 'labels' must match the number of columns in 'count_matrix'.")
  }

  # Convert labels to a factor
  labels <- as.factor(labels)

  # Add 1 to count_matrix and apply log transformation
  # Check if count_matrix is logged
  if (log) {
    count_matrix <- log1p(count_matrix)
  }

  # Compute RLE
  rle <- count_matrix - Biobase::rowMedians(count_matrix)

  # Sort RLE based on feature
  sort_index <- sort.list(labels)
  rle <- rle[, sort_index]

  # Create a data frame for RLE boxplot stats
  rledf <- t(apply(rle, 2, function(x) {
    grDevices::boxplot.stats(x)$stats
  }))
  rledf <- as.data.frame(rledf)
  colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")

  # Add feature and sample information
  rledf$feature <- labels[sort_index]
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
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4)) +
      scale_x_discrete(limits = rledf$sample) +
      scale_color_manual(values = colorRampPalette(brewer.pal(12, "Paired"))(nlevels(rledf$feature))) +
      geom_hline(yintercept = 0,
                 linetype = "dotted",
                 col = "red",
                 linewidth = 1) +
      ylab("log_expression_deviation")

    return(p)
  }, error = function(e) {
    stop("Error in plotting: ", e$message)
  })
}
