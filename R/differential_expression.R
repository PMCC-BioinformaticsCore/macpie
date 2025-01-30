#' Retrieve normalised counts of MAC-seq data
#'
#' This function retrieves  counts from a number of methods that are available
#' for normalisation, with the default of limma-voomLmFit.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param treatment_samples Value in the column "combined_id" representing replicates of treatment samples in the data
#' @param control_samples Value in the column "combined_id"  representing replicates of control samples in the data
#' @param method One of "Seurat_wilcox", "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq methods, check RUVSeq tutorial
#' @param spikes List of genes to use as spike controls
#' @importFrom limma makeContrasts eBayes contrasts.fit topTable
#' @importFrom tibble rownames_to_column
#' @import DESeq2
#' @import RUVSeq
#' @importFrom stats model.matrix
#' @importFrom dplyr rename select
#'
#' @returns Data frame of DE counts
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' treatment_samples="Staurosporine_0.1"
#' control_samples<-"DMSO_0"
#' top_table <- differential_expression(mac, treatment_samples, control_samples, method = "limma_voom")

differential_expression <- function(data = NULL,
                                    treatment_samples = NULL,
                                    control_samples = NULL,
                                    method = NULL,
                                    batch = 1,
                                    k = 2,
                                    spikes = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, method, treatment_samples, control_samples) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma_voom" else method
    if (!method %in% c("Seurat_wilcox", "DESeq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom")) {
      stop("Your normalization method is not available.")
    }
    if (is.null(treatment_samples) || is.null(control_samples)) {
      stop("Missing the vectors of treatment and control samples.")
    }
    if (!"combined_id" %in% colnames(data@meta.data)) {
      data <- data %>%
        mutate(combined_id = apply(select(starts_with("Treatment_") | starts_with("Concentration_")),
                                   1, paste, collapse = "_")) %>%
        mutate(combined_id = gsub(" ", "", .data$combined_id))
    }
    if (length(treatment_samples) == 1 && length(control_samples) == 1) {
      treatment_samples_list <- grepl(treatment_samples, data$combined_id)
      control_samples_list <- grepl(control_samples, data$combined_id)
      if (any(sum(treatment_samples_list) == 0, sum(control_samples_list) == 0)) {
        stop("The combined id of your samples (format: 'treatment'_'concentration') is not valid.")
      }
    }
  }

  # Helper: Prepare data and pheno_data
  prepare_data <- function(data, treatment_samples, control_samples, batch) {
    data <- data[, grepl(paste0(treatment_samples, "|", control_samples), data$combined_id)]
    if (length(unique(data$combined_id)) < 2) {
      stop("Insufficient factors for differential expression analysis.")
    }
    pheno_data <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
    return(list(data = data, pheno_data = pheno_data))
  }

  de_limma_voom <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)
    dge <- DGEList(counts = data@assays$RNA$counts,
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- estimateDisp(dge, model_matrix)
    dge <- calcNormFactors(dge, method = "TMMwsp")
    fit <- voomLmFit(dge, model_matrix)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    tmp <- contrasts.fit(fit, contrasts)
    tmp <- eBayes(tmp, robust = TRUE)
    top_table <- topTable(tmp, number = Inf, sort.by = "P") %>%
      select("logFC", "t", "P.Value", "adj.P.Val") %>%
      rename("log2FC" = "logFC", "metric" = "t", "p_value" = "P.Value", "p_value_adj" = "adj.P.Val") %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }

  de_edger <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)
    dge <- DGEList(counts = data@assays$RNA$counts,
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- calcNormFactors(dge, method = "TMMwsp")
    dge <- estimateDisp(dge, model_matrix)
    fit <- glmQLFit(dge, model_matrix)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    qlf <- glmQLFTest(fit, contrast = contrasts)
    top_table <- topTags(qlf, n = nrow(data)) %>%
      as.data.frame() %>%
      select("logFC", "F", "PValue", "FDR") %>%
      rename("log2FC" = "logFC", "metric" = "F", "p_value" = "PValue", "p_value_adj" = "FDR") %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }

  #DEseq produces NA adjusted p-values if
  de_deseq2 <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts,
                                  colData = pheno_data,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", treatment_samples, control_samples))
    top_table <- as.data.frame(res) %>%
      select("log2FoldChange", "stat", "pvalue", "padj") %>%
      rename("log2FC" = "log2FoldChange", "metric" = "stat", "p_value" = "pvalue", "p_value_adj" = "padj") %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }

  de_seurat <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    data <- NormalizeData(data)
    Idents(data) <- combined_id
    data$batch <- batch
    top_table <- FindMarkers(data,
                             ident.1 = treatment_samples,
                             ident.2 = control_samples,
                             latent.vars = "batch",
                             test.use = "wilcox") %>%
      select("avg_log2FC", "p_val", "p_val_adj") %>%
      rename("log2FC" = "avg_log2FC", "p_val" = "p_val", "p_value_adj" = "p_val_adj") %>%
      mutate(metric = 1) %>%
      rownames_to_column("gene")

    return(as.data.frame(top_table))
  }

  de_ruvg <- function(data, pheno_data, treatment_samples, control_samples, batch, spikes, k) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)
    if (length(spikes) == 0) {
      warning("List of control genes not provided for RUVg, using default.")
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
      stop("Some or all of your control genes are not present in the dataset.")
    }
    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(counts = as.matrix(data@assays$RNA$counts),
                               phenoData = pheno_data)
    set <- RUVg(set, cIdx = spikes, k = k)

    dge <- DGEList(counts = normCounts(set),
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- calcNormFactors(dge, method = "upperquartile")
    dge <- estimateGLMCommonDisp(dge, model_matrix)
    dge <- estimateGLMTagwiseDisp(dge, model_matrix)
    fit <- glmFit(dge, model_matrix)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    lrt <- glmLRT(fit, contrast = contrasts)
    top_table <- topTags(lrt, n = nrow(data)) %>%
      as.data.frame() %>%
      select("logFC", "PValue", "FDR") %>%
      rename("log2FC" = "logFC", "p_value" = "PValue", "p_value_adj" = "FDR") %>%
      mutate(metric = 1) %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }

  de_ruvs <- function(data, pheno_data, treatment_samples, control_samples, batch, k) {
    combined_id <- data$combined_id
    genes <- row.names(data@assays$RNA$counts)
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(counts = as.matrix(data@assays$RNA$counts),
                               phenoData = pheno_data)
    differences <- makeGroups(combined_id)
    set <- RUVs(set, cIdx = genes, k = k, scIdx = differences)

    dge <- DGEList(counts = normCounts(set),
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- calcNormFactors(dge, method = "upperquartile")
    dge <- estimateGLMCommonDisp(dge, model_matrix)
    dge <- estimateGLMTagwiseDisp(dge, model_matrix)
    fit <- glmFit(dge, model_matrix)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    lrt <- glmLRT(fit, contrast = contrasts)
    top_table <- topTags(lrt, n = nrow(data)) %>%
      as.data.frame() %>%
      select("logFC", "PValue", "FDR") %>%
      rename("log2FC" = "logFC", "p_value" = "PValue", "p_value_adj" = "FDR") %>%
      mutate(metric = 1) %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }

  de_ruvr <- function(data, pheno_data, treatment_samples, control_samples, batch, k) {
    if (ncol(data) > 100) {
      print("Warning: EdgeR with over 100 samples takes very long time. Consider reducing the number of samples.")
    }
    combined_id <- data$combined_id
    genes <- row.names(data@assays$RNA$counts)
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)

    #k defines number of sources of variation, two have been chosen for row and column
    set <- newSeqExpressionSet(counts = as.matrix(data@assays$RNA$counts),
                               phenoData = pheno_data)
    dge <- DGEList(counts = data@assays$RNA$counts,
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- calcNormFactors(dge, method = "TMMwsp")
    dge <- estimateGLMCommonDisp(dge, model_matrix)
    dge <- estimateGLMTagwiseDisp(dge, model_matrix)
    fit <- glmFit(dge, model_matrix)
    res <- residuals(fit, type = "deviance")
    set <- RUVr(set, genes, k = k, res)

    dge <- DGEList(counts = normCounts(set),
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- calcNormFactors(dge, method = "upperquartile")
    dge <- estimateGLMCommonDisp(dge, model_matrix)
    dge <- estimateGLMTagwiseDisp(dge, model_matrix)
    fit <- glmFit(dge, model_matrix)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    lrt <- glmLRT(fit, contrast = contrasts)
    top_table <- topTags(lrt, n = nrow(data)) %>%
      as.data.frame() %>%
      select("logFC", "PValue", "FDR") %>%
      rename("log2FC" = "logFC", "p_value" = "PValue", "p_value_adj" = "FDR") %>%
      mutate(metric = 1) %>%
      rownames_to_column("gene")
    return(as.data.frame(top_table))
  }


  # Main function
  validate_inputs(data, method, treatment_samples, control_samples)
  prepared <- prepare_data(data, treatment_samples, control_samples, batch)
  data <- prepared$data
  pheno_data <- prepared$pheno_data

  # Select the appropriate normalization method
  de_data <- switch(
    method,
    limma_voom = de_limma_voom(data, pheno_data, treatment_samples, control_samples),
    edgeR = de_edger(data, pheno_data, treatment_samples, control_samples),
    DESeq2 = de_deseq2(data, pheno_data, treatment_samples, control_samples),
    Seurat_wilcox = de_seurat(data, pheno_data, treatment_samples, control_samples),
    RUVg = de_ruvg(data, pheno_data, treatment_samples, control_samples, batch, spikes, k),
    RUVs = de_ruvs(data, pheno_data, treatment_samples, control_samples, batch, k),
    RUVr = de_ruvr(data, pheno_data, treatment_samples, control_samples, batch, k),
    stop("Unsupported DE method.")
  )
  return(de_data)
}
