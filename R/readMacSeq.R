
#' Load in data from Mac-seq experiments
#'
#' This functions is a derivative from Read10X() function from the Seurat package.
#' Enables easy loading of sparse data matrices provided by 10X genomics.
#' @author Mark Li, \email{justdont}
#'
#' @param data.dir Directory containing the matrix.mtx, genes.tsv (or features.tsv), and barcodes.tsv
#' files provided by nextflow pipeline using Starsolo.
#' I removed the functionality for reading in several data directories as the original Read10X()
#' This is to allow me read in the metadata label matching to the Macseq plate with minimum confusion.
#' Next version will consider open opssibilities to integrate other omics data
#' @param label.dir Directory containing the matching metadata file in csv format.
#' @param gene.column Specify which column of genes.tsv or features.tsv to use for gene names; default is 2
#' @param cell.column Specify which column of barcodes.tsv to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
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
#' expression_matrix <- Readmacseq(data.dir = data_dir)
#' macpie_object = CreateMacPie(counts = expression_matrix)
#' }
#'
readMacSeq <- function(
    data.dir,
    label.dir,
    gene.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
) {

  full.data <- list()
  has_dt <- requireNamespace("data.table", quietly = TRUE) && requireNamespace("R.utils", quietly = TRUE)
  run <- data.dir
  metadata <- label.dir

  if (!dir.exists(paths = run)) {
    stop("Macseq data directory provided does not exist")
  }
  barcode.loc <- file.path(run, 'barcodes.tsv')
  features.loc <- file.path(run, 'features.tsv.gz')
  matrix.loc <- file.path(run, 'matrix.mtx')

  if (!dir.exists(paths = metadata)) {
    stop("Metadata directory provided does not exist")
  }
  label.loc <- list.files(metadata, pattern = "\\.csv$", full.names = TRUE)

  if (!file.exists(barcode.loc)) {
    stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
  }
  if (!file.exists(features.loc) ) {
    stop("Features file missing. Expecting ", basename(path = features.loc))
  }
  if (!file.exists(matrix.loc)) {
    stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
  }
  if (!file.exists(label.loc)) {
    stop("Metadata file missing. Expecting ", basename(path = label.loc))
  }
  warning("Make sure you only have 1 metadata file in the folder.")

  data <- readMM(file = matrix.loc)
  if (has_dt) {
    cell.barcodes <- as.data.frame(data.table::fread(barcode.loc, header = FALSE))
  } else {
    cell.barcodes <- read.table(file = barcode.loc, header = FALSE, sep = '\t', row.names = NULL)
  }

  if (ncol(x = cell.barcodes) > 1) {
    cell.names <- cell.barcodes[, cell.column]
  } else {
    cell.names <- readLines(con = barcode.loc)
  }
  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }
  if (is.null(x = names(x = data.dir))) {
    if (length(x = data.dir) < 2) {
      colnames(x = data) <- cell.names
    } else {
      colnames(x = data) <- paste0(i, "_", cell.names)
    }
  } else {
    colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
  }

  if (has_dt) {
    feature.names <- as.data.frame(data.table::fread(ifelse(test = pre_ver_3, yes = gene.loc, no = features.loc), header = FALSE))
  } else {
    feature.names <- read.delim(
      file = features.loc,
      header = FALSE,
      stringsAsFactors = FALSE
    )
  }

  if (any(is.na(x = feature.names[, gene.column]))) {
    warning(
      'Some features names are NA. Replacing NA names with ID from the opposite column requested',
      call. = FALSE,
      immediate. = TRUE
    )
    na.features <- which(x = is.na(x = feature.names[, gene.column]))
    replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
    feature.names[na.features, gene.column] <- feature.names[na.features, replacement.column]
  }
  if (unique.features) {
    fcols = ncol(x = feature.names)
    if (fcols < gene.column) {
      stop(paste0("gene.column was set to ", gene.column,
                  " but feature.tsv.gz (or genes.tsv) only has ", fcols, " columns.",
                  " Try setting the gene.column argument to a value <= to ", fcols, "."))
    }
    rownames(x = data) <- make.unique(names = feature.names[, gene.column])
  }
  # In cell ranger 3.0, a third column specifying the type of data was added
  # and we will return each type of data as a separate matrix
  # I decided to leave this piece in for now
  # as I did see "Gene Expression" column in the features.tsv file from our StarSolo output
  if (ncol(x = feature.names) > 2) {
    data_types <- factor(x = feature.names$V3)
    lvls <- levels(x = data_types)
    if (length(x = lvls) > 1 && length(x = full.data) == 0) {
      message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
    }
    expr_name <- "Gene Expression"
    if (expr_name %in% lvls) { # Return Gene Expression first
      lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
    }
    data <- lapply(
      X = lvls,
      FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      }
    )
    names(x = data) <- lvls
  } else{
    data <- list(data)

  }

  full.data[[length(x = full.data) + 1]] <- data

  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    ## Fix for Issue #913
    # I need to read into this Issue as
    # currently my R complains and errors out on as.sparse
    #list_of_data[[j]] <- as.sparse(x = list_of_data[[j]])
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise
  # a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

