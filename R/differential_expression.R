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
#' fetch_normalised_counts(mac)

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
    design <- model_matrix
    fit <- voomLmFit(dge, design)
    myargs <- list(paste0("combined_id",
                          treatment_samples, "-",
                          paste0("combined_id", control_samples)),
                   levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    tmp <- contrasts.fit(fit, contrasts)
    tmp <- eBayes(tmp)
    top_table <- topTable(tmp, number = Inf) %>%
      select(logFC, P.Value, adj.P.Val) %>%
      rename("log2FC" = logFC, "p_value" = P.Value, "p_value_adj" = adj.P.Val)
    return(as.data.frame(top_table))
  }

  de_edger <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)
    group <- pheno_data$condition
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
    top_tags <- topTags(qlf, n = nrow(data)) %>%
      as.data.frame() %>%
      select(logFC, PValue, FDR) %>%
      rename("log2FC" = logFC, "p_value" = PValue, "adj.p.value" = FDR)
    return(as.data.frame(top_tags))
  }

  #DEseq produces NA adjusted p-values if
  de_deseq2 <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)
    dds <- DESeqDataSetFromMatrix(countData = data@assays$RNA$counts,
                                  colData = pheno_data,
                                  design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", treatment_samples, control_samples))
    #res.ape <- lfcShrink(dds, coef = "condition_Staurosporine_0.1_vs_DMSO_0", type = "apeglm")
    res <- as.data.frame(res) %>%
      as.data.frame() %>%
      select(log2FoldChange, pvalue, padj) %>%
      rename("log2FC" = log2FoldChange, "p_value" = pvalue, "adj.p.value" = padj)
    return(as.data.frame(top_tags))
  }

  de_seurat <- function(data, pheno_data, treatment_samples, control_samples) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
      model.matrix(~0 + combined_id + batch)

    data <- NormalizeData(data)
    Idents(data) <- combined_id
    data$batch = batch
    de.markers <- FindMarkers(data,
                              ident.1 = treatment_samples,
                              ident.2 = control_samples,
                              latent.vars = "batch",
                              test.use = "DESeq2") %>%
      select(avg_log2FC, p_val, p_val_adj) %>%
      rename("log2FC" = avg_log2FC, "p_val" = p_val, "adj.p.value" = p_val_adj)

    return(as.data.frame(de.markers))
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
    Seurat_wilcox = de_edger(data, pheno_data, treatment_samples, control_samples),
    stop("Unsupported DE method.")
  )
  return(de_data)
}
