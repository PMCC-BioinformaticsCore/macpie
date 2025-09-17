#'  Create a function to convert singlecellexperiment object to a Seurat object
#'  
#'  
#'  Convert an SCE to a Seurat object, copying a counts assay and (optionally)
#'  a log-normalized assay. Feature names are sanitized to match Seurat
#'  conventions (underscores → dashes; uniqueness enforced). Columns of
#'  `colData(sce)` become Seurat `meta.data`. The column `cell_id_col` is used
#'  for Seurat cell names.
#'  @param sce A SingleCellExperiment object.
#'  @param counts The name of the assay in the SingleCellExperiment object that contains the raw counts. 
#'  Default is "counts".
#'  @param log_counts The name of the assay in the SingleCellExperiment object that contains the 
#'  log-normalized counts.
#'  @param assay The name of the assay to create in the Seurat object. Default is "RNA".
#'  @param cell_id_col The name of the column in the colData of the SingleCellExperiment object 
#'  that contains the cell IDs.
#'  @param project_name The name of the project to assign to the Seurat object. Default is "project_name".
#'  @return A Seurat object.
#'  @importFrom SummarizedExperiment assay assayNames colData
#'  @importFrom Seurat CreateSeuratObject CreateAssayObject GetAssayData
#'  @export
#'  
#'  @examples
#' sim_counts <- matrix(rpois(10 * 20, lambda = 10), ncol = 20, nrow = 10)
#' colnames(sim_counts) <- paste0("Barcode", 1:20)
#' rownames(sim_counts) <- paste0("gene", 1:10)
#'
#'  metadata <- data.frame(
#'    Barcode = colnames(sim_counts),
#'    Condition = rep(c("A", "B"), each = 10)
#'  )
#'  sce <- SingleCellExperiment::SingleCellExperiment(
#'    assays = list(counts = sim_counts),
#'    colData = metadata
#'  )
#'  sce <- scuttle::logNormCounts(sce)
#'  sim_seurat <- sce_to_seurat(sce,
#'                              project_name = "SimulatedData")
#'  print(sim_seurat)
#'  stopifnot(inherits(sim_seurat,"Seurat"))
#'  

#' @keywords internal
#' @noRd
sanitize_features <- function(m) {
  if (is.null(m)) return(NULL)
  fn <- rownames(m)
  fn <- gsub("_", "-", fn)    # mimic Seurat underscore rule
  fn <- make.unique(fn)       # ensure uniqueness
  rownames(m) <- fn
  m
}

sce_to_seurat <- function(sce, 
                          counts = "counts",
                          log_counts = "logcounts",
                          assay = "RNA",
                          cell_id_col = "Barcode",
                          project_name = "project_name"){
  # convert SingleCellExperiment to Seurat object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("Input sce must be a SingleCellExperiment object.")
  }
  stopifnot(counts %in% SummarizedExperiment::assayNames(sce))
  cts <- SummarizedExperiment::assay(sce, counts)                   
  log_counts <- if (log_counts %in% SummarizedExperiment::assayNames(sce)) {
    SummarizedExperiment::assay(sce, log_counts)
  } else NULL
  md  <- as.data.frame(SummarizedExperiment::colData(sce))
  
  # sanitize feature names
  cts <- sanitize_features(cts)
  log_counts <- sanitize_features(log_counts)
  
  # check if cell_id_col exists in metadata
  if (!cell_id_col %in% colnames(md)) {
    stop(paste("Cell ID column", cell_id_col, "not found in metadata."))
  }
  # assign colnames to cts with cell_id_col
  colnames(cts) <- md[[cell_id_col]]
  
  # make sure cts - count matrix is a dgCMatrix
  if (!inherits(cts, "dgCMatrix")) {
    cts <- methods::as(cts, "dgCMatrix")
  }
  # create Seurat object
  seurat_obj <- CreateSeuratObject(counts = cts,
                                   meta.data = md,
                                   project = project_name)

  # add logcounts as 'logcounts' (Seurat v4) or as a layer (v5)
  if (!is.null(log_counts)) {
    if (utils::packageVersion("Seurat") < "5.0.0"
        ) {
      seurat_obj[[assay]] <- SeuratObject::SetAssayData(seurat_obj[[assay]], slot = "data", new.data = log_counts)
    } else {
      seurat_obj[[assay]] <- SeuratObject::SetAssayData(seurat_obj[[assay]], layer = "data", new.data = log_counts)
    }
  }
  return(seurat_obj)
}
