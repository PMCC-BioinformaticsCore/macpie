#' Create an RLE Plot
#'
#' This function generates a Relative Log Expression (RLE) plot for visualizing
#' the distribution of expression data after normalization or log transformation.
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param barcodes A vector of sample barcodes corresponding to
#'   Cells(seurat_object).
#' @param labels A vector of labels of the same length as 'barcodes" to group
#'   the barcodes.
#' @param label_column A metadata column name to group the barcodes.
#' @param log A logical value indicating whether data should be log-transformed.
#'   Defaults to `TRUE`.
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples.
#' @param normalisation One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "limma_trend", "limma_voom_zinb",
#'   "edgeR_zinb". If empty, defaults to raw.
#' @param spikes List of genes to use as spike controls in RUVg
#' @param num_cores Number of cores for edgeR and zinb calculations

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
#' data("mini_mac")
#' p <- plot_rle(mini_mac, label_column = "Row")
#'
#' @importFrom Biobase rowMedians
#' @importFrom dplyr pull select arrange
#' @import Seurat ggplot2 tidyseurat
#' @export

# Main function
plot_rle <- function(data, barcodes = NULL, 
                     label_column = NULL,
                     labels = NULL, 
                     log = TRUE,
                     batch = NULL, 
                     normalisation = NULL,
                     spikes = NULL,
                     num_cores = NULL) {
  
  req_pkgs <- c("colorspace", "grDevices")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "plot_rle(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }

  # Helper function to validate input data
  validate_inputs <- function(data, barcodes, label_column, labels, log, batch, 
                              normalisation, num_cores) {
    if (!inherits(data, "Seurat")) {
      stop("'data' must be a Seurat or TidySeurat object.")
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
    batch <- if (is.null(batch)) "1" else as.character(batch)
    normalisation <- if (is.null(normalisation)) "raw" else normalisation
    num_cores <- if (is.null(num_cores)) 1 else num_cores
    
    list(barcodes = as.factor(barcodes),
         labels = as.factor(labels),
         log = ifelse(inherits(log, "function"), TRUE, log),
         batch = batch,
         normalisation = normalisation,
         spikes = spikes,
         num_cores = num_cores)
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
    rle <- count_matrix - Biobase::rowMedians(as.matrix(count_matrix))
    sort_index <- sort.list(labels)
    rle <- rle[, sort_index]

    # Create RLE data frame
    rledf <- t(apply(rle, 2, function(x) {
      grDevices::boxplot.stats(x)$stats
    }))
    rledf <- as.data.frame(rledf)
    colnames(rledf) <- c("ymin", "lower", "middle", "upper", "ymax")

    # Map barcodes to well IDs using metadata
    barcode_to_well <- data@meta.data[colnames(rle), "Well_ID", drop = TRUE]  # Extract well IDs
    names(barcode_to_well) <- colnames(rle)
    
    # Add feature and sample information
    rledf$feature <- labels[sort_index]
    rledf$sample <- barcode_to_well[colnames(rle)]

    # Reorder samples for plotting
    rledf$sample <- factor(rledf$sample, levels = unique(rledf$sample))

    return(rledf)
  }

  # Compute average correlation coefficient
  compute_cv <- function(count_matrix) {

    # Convert to tibble for dplyr processing
    expression_df <- as.data.frame(count_matrix) %>% tibble::as_tibble()

    # Compute CV per sample (column-wise)
    cv_df <- expression_df %>%
      summarise(across(everything(),
                       ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE),
                       .names = "CV_{.col}"))

    avg_cv <- mean(as.numeric(cv_df), na.rm = TRUE)

    return(avg_cv)
  }


  # Helper function to create the plot
  create_rle_plot <- function(rledf, normalisation) {

    tryCatch({
      # Generate a color palette for fill
      fill_palette <- macpie_colours$discrete_40 #colorRampPalette(brewer.pal(12, "Paired"))(nlevels(rledf$feature))

      # Darken the palette for outlines
      outline_palette <- colorspace::darken(fill_palette, amount = 0.1)

      ggplot(rledf, aes(x = .data$sample, fill = .data$feature, color = .data$feature)) +
        geom_boxplot(aes(ymin = .data$ymin, lower = .data$lower, middle = .data$middle,
            upper = .data$upper, ymax = .data$ymax), stat = "identity") +
        geom_hline(yintercept = 0, linetype = "dotted", col = "red", linewidth = 1) +
        # Theme
        macpie_theme(show_x_title = TRUE, show_y_title = TRUE, legend_position_ = 'bottom', x_labels_angle = 45) +
        # Colours
        scale_x_discrete(limits = rledf$sample) +
        scale_fill_manual(values = fill_palette) +
        scale_color_manual(values = outline_palette) +
        # Titles
        geom_hline(yintercept = 0, linetype = "dotted", col = "red", linewidth = 1) +
        ggtitle(paste0("Normalisation method: ",
                       normalisation,
                       "\naverage Coef of Variation:",
                       sprintf("%.2f", compute_cv(count_matrix)))) +
        ylab("Log expression deviation")

    }, error = function(e) {
      stop("Problems in plotting ", e$message)
    })
  }

  # Validate inputs
  validated <- validate_inputs(data, barcodes, label_column, labels, log, batch, normalisation, num_cores)
  barcodes <- validated$barcodes
  labels <- validated$labels
  log <- validated$log
  batch <- validated$batch
  normalisation <- validated$normalisation
  num_cores <- validated$num_cores

  # Fetch and transform count matrix
  count_matrix <- switch(
    normalisation,
    raw = fetch_count_matrix(data, log),
    logNorm = compute_normalised_counts(data, method = "logNorm", batch = batch, num_cores = num_cores),
    cpm = compute_normalised_counts(data, method = "cpm", batch = batch, num_cores = num_cores),
    clr = compute_normalised_counts(data, method = "clr", batch = batch, num_cores = num_cores),
    SCT = compute_normalised_counts(data, method = "SCT", batch = batch, num_cores = num_cores),
    DESeq2 = compute_normalised_counts(data, method = "DESeq2", batch = batch, num_cores = num_cores),
    edgeR = compute_normalised_counts(data, method = "edgeR", batch = batch, num_cores = num_cores),
    limma_voom = compute_normalised_counts(data, method = "limma_voom", batch = batch, num_cores = num_cores),
    limma_trend = compute_normalised_counts(data, method = "limma_trend", batch = batch, num_cores = num_cores),
    RUVg = compute_normalised_counts(data, method = "RUVg", batch = batch, spikes = spikes, num_cores = num_cores),
    RUVs = compute_normalised_counts(data, method = "RUVs", batch = batch, num_cores = num_cores),
    RUVr = compute_normalised_counts(data, method = "RUVr", batch = batch, num_cores = num_cores),
    limma_voom_zinb = compute_normalised_counts(data, method = "limma_voom_zinb", batch = batch, num_cores = num_cores),
    edgeR_zinb = compute_normalised_counts(data, method = "edgeR_zinb", batch = batch, num_cores = num_cores),
    stop("Unsupported normalization method.")
  )

  # Ensure alignment
  if (length(barcodes) != ncol(count_matrix)) {
    stop("Length of 'barcodes' must match the number of columns in 'count_matrix'.")
  }
  if (length(labels) != ncol(count_matrix)) {
    stop("Length of 'labels' must match the number of columns in 'count_matrix'.")
  }

  # Create RLE data frame
  rledf <- compute_rle_df(count_matrix, labels)

  # Create and return plot
  out <- create_rle_plot(rledf, normalisation)
  out$SE <- compute_cv(count_matrix)
  return(out)
}
