#' Retrieve normalised counts of MAC-seq data
#'
#' This function retrieves  counts from a number of methods that are available
#' for normalisation, with the default of limma-voomLmFit.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param method One of "raw", "logNorm", "cpm", "clr", "SCT", "DEseq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq methods, check RUVSeq tutorial
#' @importFrom limma voom
#' @import DESeq2
#' @importFrom stats model.matrix

#' @returns Data frame of normalised counts
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' fetch_normalised_counts(mac)
fetch_normalised_counts <- function(data = NULL, method = NULL, batch = NULL, k = NULL) {
  # Helper function to validate input data
  validate_inputs <- function(data, method, batch, k) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma_voom" else method
    if (!method %in% c("raw", "logNorm",
                       "cpm", "clr", "SCT",
                       "DEseq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom")) {
      stop("Your normalization method is not available.")
    }
    batch <- if (is.null(batch)) "1" else as.character(batch)
    k <- if (is.null(k)) 2 else k
    list(data = data, batch = batch, k = k, method = method)
  }

  # Sub-functions for each normalization method
  normalize_raw <- function(data) {
    data@assays$RNA$counts
  }

  normalize_lognorm <- function(data) {
    lognorm <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
    lognorm@assays$RNA$data
  }

  normalize_cpm <- function(data) {
    cpm <- NormalizeData(data, normalization.method = "RC", scale.factor = 1e6)
    cpm@assays$RNA$data
  }

  normalize_clr <- function(data) {
    clr <- NormalizeData(data, normalization.method = "CLR", margin = 2)
    clr@assays$RNA$data
  }

  normalize_sct <- function(data, batch) {
    sct <- SCTransform(data, do.scale = TRUE,
                       return.only.var.genes = FALSE,
                       vars.to.regress = c(batch), verbose = FALSE)
    sct@assays$SCT$data
  }

  normalize_deseq2 <- function(data, batch) {
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
    design <- if (length(batch) == 1) model.matrix(~data$combined_id) else model.matrix(~data$combined_id + batch)
    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts, colData = coldata, design = design)
    dds <- estimateSizeFactors(dds)
    counts(dds, normalized = TRUE)
  }

  normalize_edger <- function(data, batch) {
    dge <- DGEList(counts = data@assays$RNA$counts, samples = data$combined_id, group = data$combined_id)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- if (length(batch) == 1) model.matrix(~data$combined_id) else model.matrix(~data$combined_id + batch)
    dge <- estimateDisp(dge, design)
    cpm(dge, log = FALSE)
  }

  normalize_limma_voom <- function(data, batch) {
    dge <- DGEList(counts = data@assays$RNA$counts, samples = data$combined_id, group = data$combined_id)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- if (length(batch) == 1) model.matrix(~data$combined_id) else model.matrix(~data$combined_id + batch)
    dge <- voom(dge, design)
    dge$E
  }

  # Main function logic
  validated <- validate_inputs(data, method, batch, k)
  data <- validated$data
  batch <- validated$batch
  k <- validated$k
  method <- validated$method
  data <- data %>%
    mutate(combined_id = apply(pick(starts_with("Treatment_") | starts_with("Concentration_")),
                               1, function(row) paste(row, collapse = "_"))) %>%
    mutate(combined_id = gsub(" ", "", .data$combined_id))

  # Select the appropriate normalization method
  norm_data <- switch(
    method,
    raw = normalize_raw(data),
    logNorm = normalize_lognorm(data),
    cpm = normalize_cpm(data),
    clr = normalize_clr(data),
    SCT = normalize_sct(data, batch),
    DEseq2 = normalize_deseq2(data, batch),
    edgeR = normalize_edger(data, batch),
    limma_voom = normalize_limma_voom(data, batch),
    stop("Unsupported normalization method.")
  )

  return(norm_data)
}
