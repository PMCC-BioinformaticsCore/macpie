test_that("sce_to_seurat works", {
  # check packages and seurat version
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scuttle")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  is_v5 <- utils::packageVersion("Seurat") >= "5.0.0"
  
  # set up test data
  set.seed(123)
  counts <- matrix(rpois(12*8, lambda = 5), nrow = 12, ncol = 8)
  rownames(counts) <- c(
    "gene_a", "gene_a",      # duplicate
    "gene_b", "gene_b",      # duplicate
    "gene_c_1", "gene_c_1",  # duplicate + underscore
    "MT-ND1", "RpS3", "RPL5", "gene_d", "gene_e", "gene_f"
  )
  barcodes <- paste0("BC", seq_len(8))
  colnames(counts) <- barcodes
  metadata <- data.frame(
    Barcode = barcodes,
    Condition = rep(c("A", "B"), each = 4),
    rownames = barcodes
  )
  # create SCE object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = metadata
  )
  # normalize (adds 'logcounts' assay)
  sce <- scuttle::logNormCounts(sce)
  
  # convert to seurat
  se <- sce_to_seurat(
    sce,
    counts = "counts",
    log_counts = "logcounts",
    assay = "RNA",
    cell_id_col = "Barcode",
    project_name = "TestProj"
  )
  
  # basic structure
  expect_s4_class(se, "Seurat")
  expect_true("RNA" %in% names(se@assays))
  
  # helper accessors that work across Seurat v4/v5
  get_counts <- function(obj, assay = "RNA") {
    if (is_v5) {
      SeuratObject::LayerData(obj[[assay]], "counts")
    } else {
      SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
    }
  }
  get_data_layer <- function(obj, assay = "RNA") {
    if (is_v5) {
      if ("data" %in% SeuratObject::Layers(obj[[assay]])) {
        SeuratObject::LayerData(obj[[assay]], "data")
      } else {
        NULL
      }
    } else {
      SeuratObject::GetAssayData(obj, assay = assay, slot = "data")
    }
  }
  
  cnt <- get_counts(se, "RNA")
  # counts should be dgCMatrix (coerced by sce_to_seurat)
  expect_true(inherits(cnt, "dgCMatrix"))
  
  # dimensions: same as input
  expect_identical(dim(cnt), dim(counts))
  
  # cell names must come from metadata Barcode
  expect_identical(colnames(cnt), metadata$Barcode)
  
  # feature name sanitation (underscores replaced; uniqueness enforced)
  expect_false(any(grepl("_", rownames(cnt), fixed = TRUE)))
  expect_false(any(duplicated(rownames(cnt))))
  
  # logcounts were added as "data" (v4: slot; v5: layer)
  dat <- get_data_layer(se, "RNA")
  if (is_v5) {
    expect_true("data" %in% SeuratObject::Layers(se[["RNA"]]))
  } else {
    # In Seurat v4, "data" slot should exist and match dims
    expect_true(is.matrix(dat) || inherits(dat, "dgCMatrix"))
  }
  
  # Compare names & shape with the sanitized input logcounts
  # (sce_to_seurat sanitizes both counts and log_counts the same way)
  logc_in <- SummarizedExperiment::assay(sce, "logcounts")
  # reproduce the same sanitation used internally
  sanitize_features_local <- function(m) {
    if (is.null(m)) return(NULL)
    fn <- rownames(m)
    fn <- gsub("_", "-", fn)
    fn <- make.unique(fn)
    rownames(m) <- fn
    m
  }
  logc_sanitized <- sanitize_features_local(logc_in)
  
  # Ensure dims & names align
  if (!is.null(dat)) {
    dat_mat <- if (inherits(dat, "dgCMatrix")) as.matrix(dat) else dat
    expect_identical(rownames(dat_mat), rownames(cnt))
    expect_identical(colnames(dat_mat), colnames(cnt))
    
    # Numeric equality with the sanitized logcounts
    ref_mat <- as.matrix(logc_sanitized[rownames(dat_mat), colnames(dat_mat), drop = FALSE])
    expect_equal(dat_mat, ref_mat, tolerance = 1e-8, check.attributes = FALSE)
  }
})

test_that("sce_to_seurat handles missing logcounts gracefully", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  is_v5 <- utils::packageVersion("Seurat") >= "5.0.0"
  
  counts <- matrix(rpois(20, 5), nrow = 5)
  rownames(counts) <- paste0("g_", 1:5)
  colnames(counts) <- paste0("BC", 1:4)
  md <- data.frame(Barcode = colnames(counts), row.names = colnames(counts))
  
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = md
  )
  
  se <- sce_to_seurat(sce, log_counts = "DOES_NOT_EXIST", project_name = "NoData")
  
  # counts exist
  if (is_v5) {
    cnt <- SeuratObject::LayerData(se[["RNA"]], "counts")
    expect_true(inherits(cnt, "dgCMatrix"))
    # no "data" layer should have been added
    expect_false("data" %in% SeuratObject::Layers(se[["RNA"]]))
  } else {
    cnt <- SeuratObject::GetAssayData(se, assay = "RNA", slot = "counts")
    expect_true(inherits(cnt, "dgCMatrix"))
    # v4: data slot should be empty (0x0) or absent-like
    dat <- SeuratObject::GetAssayData(se, assay = "RNA", slot = "data")
    expect_true(nrow(dat) == 0 || ncol(dat) == 0)
  }
})

test_that("sce_to_seurat errors clearly on bad inputs", {
  skip_if_not_installed("SingleCellExperiment")
  
  counts <- matrix(rpois(20, 5), nrow = 5)
  rownames(counts) <- paste0("g", 1:5)
  colnames(counts) <- paste0("BC", 1:4)
  md <- data.frame(Barcode = colnames(counts), row.names = colnames(counts))
  
  # not an SCE
  expect_error(
    sce_to_seurat(counts),
    "SingleCellExperiment",
    fixed = FALSE
  )
  
  # missing cell_id_col
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = data.frame(NotBarcode = colnames(counts), row.names = colnames(counts))
  )
  expect_error(
    sce_to_seurat(sce, cell_id_col = "Barcode"),
    "Cell ID column.*not found",
    ignore.case = TRUE
  )
  
  # counts assay name missing
  sce2 <- SingleCellExperiment::SingleCellExperiment(
    assays = list(other = counts),
    colData = md
  )
  expect_error(
    sce_to_seurat(sce2, counts = "counts"),
    "counts",
    fixed = FALSE
  )
})


