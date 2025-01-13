#' Retrieve normalised counts of MAC-seq data
#'
#' This function retrieves  counts from a number of methods that are available
#' for normalisation, with the default of limma-voomLmFit.
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param treatment_samples A logical vector representing replicates of treatment samples in the data
#' @param control_samples A logical vector representing replicates of control samples in the data
#' @param method One of "Seurat", "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq methods, check RUVSeq tutorial
#' @param spikes List of genes to use as spike controls
#' @importFrom limma voom makeContrasts
#' @import DESeq2
#' @import RUVSeq
#' @importFrom EDASeq newSeqExpressionSet normCounts
#' @importFrom Biobase pData
#' @importFrom stats model.matrix residuals
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
                                    batch = NULL,
                                    k = NULL,
                                    spikes = NULL) {
  # Helper function to validate input data
  validate_inputs <- function(data, treatment_samples, control_samples, method, k) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    method <- if (is.null(method)) "limma_voom" else method
    if (!method %in% c("SCT", "DESeq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom")) {
      stop("Your normalization method is not available.")
    }
    if(is.null(treatment_samples) || is.null(control_samples)){
      stop("Missing the vectors of treatment and control samples.")
    }

    if(length(treatment_samples) == 1 && length(control_samples) == 1){
      treatment_samples_list <- grepl(treatment_samples, data$combined_id)
      control_samples_list <- grepl(control_samples, data$combined_id)
      if(any(sum(treatment_samples_list) == 0, sum(control_samples_list) == 0)){
        stop("The combined id of your samples (format: Treatment_Concentration) is not valid.")
      }
    }

    k <- if (is.null(k)) 2 else k
    list(data = data,
         treatment_samples_list = treatment_samples,
         control_samples_list = control_samples,
         k = k, method = method, spikes = spikes)
  }

  de_limma_voom <- function(data, model_matrix, pheno_data, treatment_samples, control_samples) {
    dge <- DGEList(counts = data@assays$RNA$counts,
                   samples = pheno_data$condition,
                   group = pheno_data$condition)
    dge <- estimateDisp(dge, model_matrix)
    dge <- calcNormFactors(dge, methods = "TMMwsp")
    design <- model_matrix
    fit <- glmQLFit(dge, design)
    myargs = list(paste0("combined_id", treatment_samples, "-", paste0("combined_id", control_samples)), levels = model_matrix)
    contrasts <- do.call(makeContrasts, myargs)
    qlf <- glmQLFTest(fit, contrast=contrasts)
    topTags <- topTags(qlf,n=length(qlf$df.total))
    return(as.data.frame(topTags))
  }

  # Main function

  #create an unique identifier based on combined annotation
  batch <- if (is.null(batch)) "1" else as.character(batch)
  if(!any(colnames(data@meta.data) %in% "combined_id")){
    data <- data %>%
    mutate(combined_id = apply(pick(starts_with("Treatment_") | starts_with("Concentration_")),
                               1, function(row) paste(row, collapse = "_"))) %>%
    mutate(combined_id = gsub(" ", "", .data$combined_id))
  }

  #validate inputs
  validated <- validate_inputs(data, treatment_samples, control_samples, method, k)
  data <- validated$data
  treatment_samples_list <- validated$treatment_samples_list
  control_samples_list <- validated$control_samples_list
  method <- validated$method
  k <- validated$k

  #define pheno data and model matrix
  data<-data[,grepl(paste0(treatment_samples,"|",control_samples), data$combined_id)]

  #replicates but only one data type
  if (length(unique(data$combined_id)) == 1) {
    stop("Only one factor given for the DE analysis.")
  } else if (length(unique(data$combined_id)) > 1) {
    combined_id <- data$combined_id
    model_matrix <- if (length(batch) == 1) model.matrix(~0+combined_id) else model.matrix(~0+combined_id)
    pheno_data <- data.frame(batch = as.factor(batch), condition = as.factor(data$combined_id))
  } else {
    stop("Insufficient number of factors for definition of model matrix.")
  }

    # Select the appropriate normalization method
  de_data <- switch(
    method,
    limma_voom = de_limma_voom(data = data,
                               model_matrix = model_matrix,
                               pheno_data = pheno_data,
                               treatment_samples = treatment_samples,
                               control_samples = control_samples),
    stop("Unsupported DE method.")
  )
  return(de_data)
}
