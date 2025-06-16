#' Model Gene Dose-Response Curves Across Treatments
#'
#' This function fits dose-response models for a set of genes across different treatments
#' using the `drc` package. It returns EC50 values per gene per treatment.
#'
#' @param data A Seurat or TidySeurat object containing expression data and metadata.
#' @param genes A character vector of gene names to model. If NULL, all significant DE genes across comparisons are used.
#' @param normalisation A character string indicating the normalization method. One of: "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "zinb". Default is "limma_voom".
#' @param control_value A string indicating the control condition in "Treatment_1". Default is "DMSO".
#' @param batch Batch variable to use for normalization if applicable. Default is 1.
#' @param k Number of unwanted factors for RUV normalization. Default is 2.
#' @param num_cores Number of CPU cores to use in parallel model fitting. Default is 1.
#'
#' @import drc
#' @importFrom dplyr filter mutate select pull distinct bind_rows
#' @importFrom parallel mclapply
#' @importFrom utils capture.output
#' @return A data frame of EC50 values per gene and treatment.
#'
#' @examples
#' \dontrun{
#' rds_file <- system.file("extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(rds_file)
#' res <- compute_multiple_dose_response(
#'   data = mac,
#'   genes = c("PTPRA", "MYC"),
#'   normalisation = "limma_voom",
#'   treatment_value = "Camptothecin",
#'   num_cores = 2
#' )
#' head(res)
#' }
#' @export
compute_multiple_dose_response <- function(data,
                                           genes = NULL,
                                           normalisation = "limma_voom",
                                           control_value = "DMSO",
                                           batch = 1,
                                           k = 2,
                                           num_cores = 1) {
  
  if (!inherits(data, "Seurat")) {
    stop("Input data must be a Seurat object")
  }
  
  supported_methods <- c("raw", "logNorm", "cpm", "clr", "SCT",
                         "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr",
                         "limma_voom", "zinb")
  if (!normalisation %in% supported_methods) {
    stop("Unsupported normalization method.")
  }
  
  assay_data <- switch(
    normalisation,
    raw = as.matrix(data@assays$RNA$counts),
    logNorm = compute_normalised_counts(data, method = "logNorm", batch = batch),
    cpm = compute_normalised_counts(data, method = "cpm", batch = batch),
    clr = compute_normalised_counts(data, method = "clr", batch = batch),
    SCT = compute_normalised_counts(data, method = "SCT", batch = batch),
    DESeq2 = compute_normalised_counts(data, method = "DESeq2", batch = batch),
    edgeR = compute_normalised_counts(data, method = "edgeR", batch = batch),
    limma_voom = compute_normalised_counts(data, method = "limma_voom", batch = batch),
    RUVg = compute_normalised_counts(data, method = "RUVg", batch = batch, spikes = spikes),
    RUVs = compute_normalised_counts(data, method = "RUVs", batch = batch),
    RUVr = compute_normalised_counts(data, method = "RUVr", batch = batch),
    zinb = compute_normalised_counts(data, method = "zinb", batch = batch)
  )
  
  treatments <- data %>%
    select("Treatment_1") %>%
    filter(!grepl("control_value", "Treatment_1")) %>%
    pull() %>%
    unique()
  
  if (is.null(genes)) {
    all_de <- bind_rows(data@tools$diff_exprs, .id = "comparison")
    genes <- all_de %>%
      filter(p_value_adj < 0.001) %>%
      distinct(gene) %>%
      pull(gene)
  }
  
  results <- list()
  for (treatment_value in treatments) {
    meta <- data@meta.data %>%
      mutate(barcode = rownames(.)) %>%
      filter(Treatment_1 %in% c(treatment_value, control_value))
    
    if (nrow(meta) < 3) next
    
    results_genes <- mclapply(genes, function(gene) {
      if (!gene %in% rownames(assay_data)) return(NULL)
      expr <- as.numeric(assay_data[gene, meta$barcode])
      names(expr) <- meta$combined_id
      if (all(expr == 0)) return(NULL)
      
      meta$concentration <- ifelse(meta$Treatment_1 == control_value, 0,
                                   as.numeric(as.character(meta$Concentration_1)))
      df <- data.frame(
        expression = expr,
        concentration = meta$concentration,
        replicate = meta$combined_id
      )
      
      model <- tryCatch({
        suppressWarnings(suppressMessages(
          drm(expression ~ concentration, data = df, fct = LL.4())
        ))
      }, error = function(e) NULL)
      
      if (!is.null(model)) {
        ed <- tryCatch({
          suppressWarnings(suppressMessages(
            invisible(capture.output(ED(model, 50, interval = "delta")))
          ))
        }, error = function(e) NA)
        
        tokens <- strsplit(ed[5], "\\s+")[[1]]
        tokens <- tokens[tokens != ""]
        second_value <- tokens[2]
        ec50 <- if (length(second_value) == 1) as.numeric(second_value) else NA
      } else {
        ec50 <- NA
      }
      
      return(ec50)
    }, mc.cores = num_cores)
    
    names(results_genes) <- genes
    results[[treatment_value]] <- unlist(results_genes)
  }
  
  df <- do.call("cbind", results)
  return(df)
}
