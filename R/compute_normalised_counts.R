#' Retrieve normalised counts of MAC-seq data
#'
#' This function retrieves  counts from a number of methods that are available
#' for normalisation, with the default of limma-voomLmFit.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param method One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "limma_trend", "zinb"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq and zinb methods
#' @param spikes List of genes to use as spike controls
#' @param max_counts Maximum count for a gene across all samples
#' @param num_cores Number of cores
#' @importFrom limma voom
#' @importFrom Seurat as.SingleCellExperiment
#' @import DESeq2
#' @import RUVSeq
#' @importFrom Biobase pData
#' @importFrom stats model.matrix residuals
#' @importFrom limma lmFit
#' @importFrom parallel makeCluster stopCluster
#'
#' @returns Data frame of normalised counts
#' @export
#'
#' @examples
#' data(mini_mac)
#' compute_normalised_counts(mini_mac)
compute_normalised_counts <- function(data = NULL,
                                    method = NULL,
                                    batch = NULL,
                                    k = NULL,
                                    spikes = NULL,
                                    max_counts = NULL,
                                    num_cores = NULL) {
  req_pkgs <- c("SingleCellExperiment", "EDASeq", "BiocParallel", "doParallel","zinbwave")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "compute_normalised_counts(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }
  # Helper function to validate input data
  validate_inputs <- function(data, method, batch, k, max_counts, num_cores) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma_voom" else method
    if (!method %in% c("raw", "logNorm",
                       "cpm", "clr", "SCT",
                       "DESeq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom", "limma_trend", "zinb")) {
      stop("Your normalization method is not available.")
    }
    batch <- if (is.null(batch)) "1" else as.character(batch)
    k <- if (is.null(k)) 2 else k
    max_counts <- if (is.null(max_counts)) 100 else as.numeric(max_counts)
    num_cores <- if (is.null(num_cores)) 1 else num_cores
    list(data = data, batch = batch, k = k, method = method, num_cores=num_cores)
  }

  #create an unique identifier based on combined annotation
  batch <- if (is.null(batch)) "1" else as.character(batch)
  data <- data %>%
    mutate(combined_id = apply(pick(starts_with("Treatment_") | starts_with("Concentration_")),
                               1, function(row) paste(row, collapse = "_"))) %>%
    mutate(combined_id = gsub(" ", "", .data$combined_id))

  #define pheno data and model matrix
  #replicates but only one data type
  model_matrix <- NULL
  if (length(batch) == 1 && length(unique(data$combined_id)) == 1) {
    model_matrix <- matrix(1, ncol = 1, nrow = ncol(data))
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(colnames(data)))
  } else if (length(batch) == 1 && (length(unique(data$combined_id)) > 1)) {
    model_matrix <- model.matrix(~data$combined_id)
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
  } else if (length(batch) %% length(colnames(data)) == 0 && length(unique(data$combined_id)) > 1) {
    model_matrix <- model.matrix(~data$combined_id + batch)
    coldata <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
  } else if (length(batch) %% length(colnames(data)) == 0 && length(unique(data$combined_id == 1))) {
    model_matrix <- model.matrix(~batch)
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
    return(lognorm@assays$RNA$data)
  }

  normalize_cpm <- function(data) {
    cpm <- NormalizeData(data, normalization.method = "RC", scale.factor = 1e6)
    return(cpm@assays$RNA$data)
  }

  normalize_clr <- function(data) {
    clr <- NormalizeData(data, normalization.method = "CLR", margin = 2)
    return(clr@assays$RNA$data)
  }

  normalize_sct <- function(data, batch) {
    sct <- SCTransform(data, do.scale = TRUE,
                       return.only.var.genes = FALSE,
                       vars.to.regress = if (length(batch) == 1) NULL else batch, verbose = FALSE)
    return(sct@assays$SCT$data)
  }

  normalize_deseq2 <- function(data, batch) {
    coldata <- coldata
    design <- model_matrix
    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts, colData = coldata, design = design)
    dds <- estimateSizeFactors(dds)
    return(counts(dds, normalized = TRUE))
  }

  normalize_edger <- function(data, batch) {
    if (ncol(data) > 100) {
      message("EdgeR with over 100 samples takes a long time. Consider reducing the number of samples or genes.")
    }
    cl <- makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    p <- BiocParallel::DoparParam()
    dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMM")
    design <- model_matrix
    dge <- estimateDisp(dge, design, BPPARAM = p)
    return(cpm(dge, log = FALSE))
  }

  normalize_limma_voom <- function(data, batch) {
    dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- model_matrix
    dge <- voom(dge, design)
    return(dge$E)
  }

  normalize_limma_trend <- function(data, batch) {
    dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- model_matrix
    logCPM <- cpm(dge, log=TRUE, prior.count=3)
    fit <- lmFit(logCPM, design)
    fit <- eBayes(fit, trend=TRUE)
    dge <- fit
    return(logCPM)
  }
  
  normalize_ruvg <- function(data, batch, spikes, k) {
    counts <- data@assays$RNA$counts
    if (length(spikes) == 0) {
      warning("List of control genes not provided for RUVg, using default human housekeeping genes.")
      spikes <- c(
        "ACTB",
        "GAPDH",
        "RPLP0",
        "B2M",
        "HPRT1",
        "PGK1",
        "TBP",
        "UBC",
        "YWHAZ",
        "PPIA",
        "RPL19",
        "EEF1A1",
        "RPS18",
        "TFRC"
      )
    }
    if (!all(spikes %in% row.names(data@assays$RNA$counts))) {
      warning("Some or all of your control genes are not present in the dataset.")
      spikes <- intersect(spikes, rownames(counts(set)))
    }
    #k defines number of sources of variation, two have been chosen for row and column
    set <- EDASeq::newSeqExpressionSet(counts = as.matrix(counts),
                               phenoData = data.frame(condition = coldata$condition,
                                                      row.names = colnames(counts)))
    set <- RUVg(set, cIdx = spikes, k = k)
    EDASeq::normCounts(set)
  }

  normalize_ruvs <- function(data, batch, k) {
    counts <- data@assays$RNA$counts
    genes <- rownames(counts)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- EDASeq::newSeqExpressionSet(as.matrix(counts),
                               phenoData = data.frame(condition = coldata$condition,
                                                      row.names = colnames(counts)))
    differences <- model_matrix
    differences <- lapply(seq_len(ncol(model_matrix)), function(j) {
      which(model_matrix[, j] == 1)
    })
    max_len <- max(lengths(differences))
    scIdx <- matrix(NA, nrow = max_len, ncol = length(differences))
    for (i in seq_along(differences)) {
      scIdx[seq_along(differences[[i]]), i] <- differences[[i]]
    }
    scIdx<-t(scIdx)
    colnames(scIdx) <- colnames(model_matrix)
    set <- RUVs(set, cIdx = genes, k = k, scIdx = scIdx)

    EDASeq::normCounts(set)
  }

  normalize_ruvr <- function(data, batch, k) {
    if (ncol(data) > 100) {
      message("EdgeR with over 100 samples takes very long time. Consider reducing the number of samples.")
    }
    counts <- as.matrix(data@assays$RNA$counts)
    genes <- rownames(counts)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- EDASeq::newSeqExpressionSet(counts,
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
    EDASeq::normCounts(set)
  }

  normalize_zinb <- function(data, batch) {

    message("Please allow extra time for zinb mode.")
    if (ncol(data) > 50) {
      message("zinb with over 50 samples takes a long time. Consider reducing the number of samples or genes.")
    }
    data_sce <- as.SingleCellExperiment(data)
    filtered_sce <- subset(data_sce, rowSums(as.data.frame(counts(data_sce))) > 0)
    cl <- makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    p <- BiocParallel::DoparParam()
    system.time(zinb <- zinbwave::zinbwave(filtered_sce, K = k,
                                 epsilon=12000,
                                 BPPARAM = p,
                                 observationalWeights = TRUE))
    counts <- zinb@assays@data$counts
    weights <- zinb@assays@data$weights 
    
    # approximate denoised counts (downweighting dropouts)
    #dge <- DGEList(counts = data@assays$RNA$counts, samples = coldata$condition, group = coldata$condition)
    #dge <- calcNormFactors(dge, methods = "TMM")
    #design <- model_matrix
    #dge <- estimateDisp(dge, design, BPPARAM = p)
    #dge$weights <- assay(zinb, "weights")
    #fit <- glmQLFit(dge, design, BPPARAM = p)
    #norm_counts <- fitted(fit)
    
    dge <- DGEList(counts = counts(filtered_sce), samples = coldata$condition, group = coldata$condition)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- model_matrix
    dge <- voom(dge, design, weights = weights)
    normalised_values <- dge$E
    
    #normalised_values <- zinb@assays@data$normalizedValues
    stopCluster(cl)
    doParallel::registerDoParallel()
    return(normalised_values)
  }

  # Main function logic
  validated <- validate_inputs(data, method, batch, k, max_counts, num_cores)
  data <- validated$data
  batch <- validated$batch
  k <- validated$k
  method <- validated$method
  num_cores <- validated$num_cores
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
    limma_trend = normalize_limma_trend(data, batch),
    RUVg = normalize_ruvg(data, batch, spikes, k),
    RUVs = normalize_ruvs(data, batch, k),
    RUVr = normalize_ruvr(data, batch, k),
    zinb = normalize_zinb(data, batch),
    stop("Unsupported normalization method.")
  )
  return(norm_data)
}
