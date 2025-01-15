#' Retrieve normalised counts of MAC-seq data
#'
#' This function retrieves  counts from a number of methods that are available
#' for normalisation, with the default of limma-voomLmFit.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param method One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq methods, check RUVSeq tutorial
#' @param spikes List of genes to use as spike controls
#' @importFrom limma voom
#' @import DESeq2
#' @import RUVSeq
#' @importFrom EDASeq newSeqExpressionSet normCounts
#' @importFrom Biobase pData
#' @importFrom stats model.matrix residuals
#'
#' @returns Data frame of normalised counts
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' fetch_normalised_counts(mac)
fetch_normalised_counts <- function(data = NULL,
                                    method = NULL,
                                    batch = NULL,
                                    k = NULL,
                                    spikes = NULL) {
  # Helper function to validate input data
  validate_inputs <- function(data, method, batch, k) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma_voom" else method
    if (!method %in% c("raw", "logNorm",
                       "cpm", "clr", "SCT",
                       "DESeq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom")) {
      stop("Your normalization method is not available.")
    }
    batch <- if (is.null(batch)) "1" else as.character(batch)
    k <- if (is.null(k)) 2 else k
    list(data = data, batch = batch, k = k, method = method)
  }

  #create an unique identifier based on combined annotation
  batch <- if (is.null(batch)) "1" else as.character(batch)
  data <- data %>%
    mutate(combined_id = apply(pick(starts_with("Treatment_") | starts_with("Concentration_")),
                               1, function(row) paste(row, collapse = "_"))) %>%
    mutate(combined_id = gsub(" ", "", .data$combined_id))

  #define pheno data and model matrix
  #replicates but only one data type
  if (length(batch) == 1 && length(unique(data$combined_id)) == 1) {
    model_matrix <- matrix(1, ncol = 1, nrow = ncol(data))
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(colnames(data)))
  } else if (length(batch) == 1 && (length(unique(data$combined_id)) > 1)) {
    model_matrix <- model.matrix(~data$combined_id)
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
  } else if (length(batch) %% length(colnames(data)) == 0 && length(unique(data$combined_id)) > 1) {
    model_matrix <- model.matrix(~data$combined_id + batch)
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
  } else {
    stop("Insufficient number of factors for definition of model matrix.")
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
                       vars.to.regress = if (length(batch) == 1) NULL else batch, verbose = FALSE)
    sct@assays$SCT$data
  }

  normalize_deseq2 <- function(data, batch) {
    coldata <- coldata
    design <- model_matrix
    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts, colData = coldata, design = design)
    dds <- estimateSizeFactors(dds)
    counts(dds, normalized = TRUE)
  }

  normalize_edger <- function(data, batch) {
    if (ncol(data) > 100) {
      print("Warning: EdgeR with over 100 samples takes a long time. Consider reducing the number of samples or genes.")
    }
    dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMM")
    design <- model_matrix
    dge <- estimateDisp(dge, design)
    cpm(dge, log = FALSE)
  }

  normalize_limma_voom <- function(data, batch) {
    dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- model_matrix
    dge <- voom(dge, design)
    dge$E
  }

  normalize_ruvg <- function(data, batch, spikes, k) {
    counts <- data@assays$RNA$counts
    if (length(spikes) == 0) {
      stop("List of control genes not provided for RUVg.")
    }
    if (!all(spikes %in% row.names(data@assays$RNA$counts))) {
      stop("Some or all of your control genes are not present in the dataset.")
    }
    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(counts = as.matrix(counts),
                               phenoData = data.frame(condition = coldata$condition,
                                                      row.names = colnames(counts)))
    set <- RUVg(set, cIdx = spikes, k = k)
    normCounts(set)
  }

  normalize_ruvs <- function(data, batch, k) {
    counts <- data@assays$RNA$counts
    genes <- rownames(counts)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(as.matrix(counts),
                               phenoData = data.frame(condition = coldata$condition,
                                                      row.names = colnames(counts)))
    differences <- model_matrix
    set <- RUVs(set, cIdx = genes, k = k, scIdx = differences)

    normCounts(set)
  }

  normalize_ruvr <- function(data, batch, k) {
    if (ncol(data) > 100) {
      print("Warning: EdgeR with over 100 samples takes very long time. Consider reducing the number of samples.")
    }
    counts <- as.matrix(data@assays$RNA$counts)
    genes <- rownames(counts)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(counts,
                               phenoData = data.frame(condition = coldata$condition,
                                                      row.names = colnames(counts)))
    design <- model_matrix
    y <- DGEList(counts = counts(set), group = pData(set)$condition)
    y <- calcNormFactors(y, method = "TMMwsp")
    y <- estimateGLMCommonDisp(y, design)
    y <- estimateGLMTagwiseDisp(y, design)
    fit <- glmFit(y, design)
    res <- residuals(fit, type = "deviance")
    set <- RUVr(set, genes, k = k, res)
    normCounts(set)
  }

  # Main function logic
  validated <- validate_inputs(data, method, batch, k)
  data <- validated$data
  batch <- validated$batch
  k <- validated$k
  method <- validated$method

  # Select the appropriate normalization method
  norm_data <- switch(
    method,
    raw = normalize_raw(data),
    logNorm = normalize_lognorm(data),
    cpm = normalize_cpm(data),
    clr = normalize_clr(data),
    SCT = normalize_sct(data, batch),
    DESeq2 = normalize_deseq2(data, batch),
    edgeR = normalize_edger(data, batch),
    limma_voom = normalize_limma_voom(data, batch),
    RUVg = normalize_ruvg(data, batch, spikes, k),
    RUVs = normalize_ruvs(data, batch, k),
    RUVr = normalize_ruvr(data, batch, k),
    stop("Unsupported normalization method.")
  )
  return(norm_data)
}
