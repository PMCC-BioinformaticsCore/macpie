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
#' \dontrun{
#' # Example Data
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' rle_plot(data = mac, label_column = "Row", log=TRUE)
#'
#' @importFrom Biobase rowMedians
#' @importFrom grDevices boxplot.stats
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr pull select
#' @import Seurat ggplot2 tidyseurat
#' @importFrom RColorBrewer brewer.pal
#' @importFrom colorspace darken
#' @export

# Main function
rle_plot <- function(data, barcodes = NULL, label_column = NULL, labels = NULL, log = TRUE) {

  # Helper function to validate input data
  validate_inputs <- function(data, barcodes, label_column, labels, log) {
    if (!inherits(data, "Seurat")) {
      stop("Error: 'data' must be a Seurat or TidySeurat object.")
    }
    if (is.null(barcodes)) {
      barcodes <- Seurat::Cells(data)
    }
    if (!is.null(label_column) && inherits(label_column, "character") && length(label_column) == 1) {
      labels <- data %>%
        select({{label_column}}) %>%
        pull()
    } else if (is.null(labels)) {
      stop("The format of `label_column` should be a single character value.")
    }
    if (is.null(labels)) {
      stop("Either `labels` or `label_column` must be provided.")
    }
    if (length(labels) != ncol(data)) {
      stop("Labels must have the same length as the number of columns in the dataset.")
    }
    list(barcodes = as.factor(barcodes), labels = as.factor(labels), log = ifelse(inherits(log, "function"), TRUE, log))
  }

  # Helper function to fetch and log-transform count matrix
  fetch_count_matrix <- function(data, log) {
    count_matrix <- as.matrix(data@assays$RNA$counts)
    if (log) {
      count_matrix <- log1p(count_matrix)
    }
    return(count_matrix)
  }

  # Helper function to compute RLE
  compute_rle_df <- function(count_matrix, labels) {
    # Compute RLE
    rle <- count_matrix - Biobase::rowMedians(count_matrix)
    sort_index <- sort.list(labels)
    rle <- rle[, sort_index]

    # Create RLE data frame
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

    return(rledf)
  }

  # Helper function to create the plot
  create_rle_plot <- function(rledf) {
    tryCatch({
      # Generate a color palette for fill
      fill_palette <- colorRampPalette(brewer.pal(12, "Paired"))(nlevels(rledf$feature))

      # Darken the palette for outlines
      outline_palette <- darken(fill_palette, amount = 0.1)

      ggplot(rledf, aes(x = .data$sample, fill = .data$feature, color = .data$feature)) +
        geom_boxplot(
          aes(
            ymin = .data$ymin,
            lower = .data$lower,
            middle = .data$middle,
            upper = .data$upper,
            ymax = .data$ymax
          ),
          stat = "identity"
        ) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4)) +
        scale_x_discrete(limits = rledf$sample) +
        scale_fill_manual(values = fill_palette) +
        scale_color_manual(values = outline_palette) +
        geom_hline(yintercept = 0, linetype = "dotted", col = "red", linewidth = 1) +
        ylab("Log expression deviation")
    }, error = function(e) {
      stop("Error in plotting: ", e$message)
    })
  }

  # Validate inputs
  validated <- validate_inputs(data, barcodes, label_column, labels, log)
  barcodes <- validated$barcodes
  labels <- validated$labels
  log <- validated$log

  # Fetch and transform count matrix
  count_matrix <- fetch_count_matrix(data, log)

  # Ensure alignment
  if (length(barcodes) != ncol(count_matrix)) {
    stop("Error: Length of 'barcodes' must match the number of columns in 'count_matrix'.")
  }
  if (length(labels) != ncol(count_matrix)) {
    stop("Error: Length of 'labels' must match the number of columns in 'count_matrix'.")
  }

  # Create RLE data frame
  rledf <- compute_rle_df(count_matrix, labels)

  # Create and return plot
  create_rle_plot(rledf)
}
