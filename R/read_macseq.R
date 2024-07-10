
#' Load in data from Mac-seq experiments
#'
#' This functions is a derivative from Read10X function from the Seurat package.
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#' @author Mark Li, \email{justdont}
#'
#' @param data_dir Directory containing the matrix.mtx, genes.tsv
#' (or features.tsv), and barcodes.tsv files provided by nextflow
#' pipeline using Starsolo.
#' Functionality for reading several data directories was removed
#' to allow reading in the metadata label matching to the Macseq
#' plate with minimum confusion.
#' Next version will consider integration of other modalities.
#'
#' @param label_file File containing the matching metadata in csv format.
#' @param gene_column Specify which column of genes.tsv or features.tsv
#' to use for gene names; default is 2
#' @param cell_column Specify which column of barcodes.tsv
#' to use for cell names; default is 1
#' @param unique_features Make feature names unique (default TRUE)
#' @param strip_suffix Remove trailing "-1" if present in all cell barcodes.
#'
#'
#' @return If features.csv indicates the data has multiple data types, a list
#'   containing a sparse matrix of the data from each type will be returned.
#'   Otherwise a sparse matrix containing the expression data will be returned.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' #
#' data_dir <- 'path/to/data/directory'
#' list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
#' expression_matrix <- Readmacseq(data_dir = data_dir)
#' macpie_object = CreateMacPie(counts = expression_matrix)
#' }
#'
read_macseq <- function(
  data_dir,
  label_file,
  gene_column = 2,
  cell_column = 1,
  unique_features = TRUE,
  strip_suffix = FALSE
) {

  run <- data_dir
  metadata <- label_file

  #check for existance of input files
  if (!dir.exists(paths = run)) {
    stop("Macseq data directory provided does not exist")
  }
  barcode_loc <- file.path(run, "barcodes.tsv.gz")
  features_loc <- file.path(run, "features.tsv.gz")
  matrix_loc <- file.path(run, "matrix.mtx.gz")

  if (!file.exists(label_file)) {
    stop("Metadata file provided does not exist")
  }
  label_file <- list.files(metadata, pattern = "\\.csv$", full.names = TRUE)

  if (!any(file.exists(c(barcode_loc, features_loc, matrix_loc)))) {
    stop(paste(
               "Missing one or more of count matrix files. Expecting",
               basename(path = barcode_loc),
               basename(path = features_loc),
               basename(path = matrix_loc))
    )
  }

  if (!file.exists(label_file)) {
    stop("Metadata file missing. Expecting ", label_file)
  }

  #read in the matrix
  data <- readMM(file = matrix_loc)
  cell_barcodes <- as.data.frame(data.table::fread(barcode_loc, header = FALSE))

  #define cellular barcodes and strip suffix if necessary
  if (ncol(x = cell_barcodes) > 1) {
    cell_names <- cell_barcodes[, cell_column]
  } else {
    cell_names <- readLines(con = barcode_loc)
  }
  if (all(grepl(pattern = "\\-1$", x = cell_names)) && strip_suffix) {
    cell_names <- as.vector(x = as.character(x = sapply(
      X = cell_names,
      FUN = Seurat::ExtractField,
      field = 1,
      delim = "-"
    )))
  }

  #assign column names as cell barcodes
  colnames(x = data) <- cell_names

  feature_names <- as.data.frame(
    data.table::fread(features_loc, header = FALSE)
  )

  if (any(is.na(x = feature_names[, gene_column]))) {
    warning(
      "Some features names are NA. Replacing NA names with
      ID from the opposite column requested.",
      call. = FALSE,
      immediate. = TRUE
    )
    na_features <- which(is.na(feature_names[, gene_column]))
    replacement_column <- ifelse(gene_column == 2, yes = 1, no = 2)
    feature_names[na_features, gene_column] <-
      feature_names[na_features, replacement_column]
  }

  if (unique_features) {
    fcols <- ncol(x = feature_names)
    if (fcols < gene_column) {
      stop(
        paste0(
               "gene_column was set to ",
               gene_column,
               " but feature.tsv.gz (or genes.tsv) only has ",
               fcols,
               " columns. Try setting the gene_column argument to a value <= to ",
               fcols, ".")
      )
    }
    rownames(x = data) <- make.unique(names = feature_names[, gene_column])
  }

  list_of_data <- data
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}
