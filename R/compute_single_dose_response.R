#' Model Gene Dose-Response Curve Using drc
#'
#' @param data A Seurat or TidySeurat object containing expression and metadata.
#' @param gene A gene name (must match a row name in the object).
#' @param normalisation One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "zinb". If empty, defaults to cpm
#' @param treatment_value A character string matching one value in metadata column "Treatment_1".
#' @import drc
#' @import tidyverse
#' @return A list with drc model, predicted values, and ggplot curve
#' @examples
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' res <- compute_single_dose_response(data = mac, gene = "SOX12", normalisation = "cpm", treatment_value = "Camptothecin")
#' res$plot
#' }
#' 
#' @export
compute_single_dose_response <- function(data, gene, normalisation = "cpm", treatment_value, batch=1) {
  
  # Sanity checks
  if (!gene %in% rownames(data)) stop("Gene not found in expression matrix.")
  if (!"Treatment_1" %in% colnames(data@meta.data)) stop("Metadata must include 'Treatment_1' column.")
  # Subset metadata
  meta <- data@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    filter(Treatment_1 == treatment_value | Treatment_1 == "DMSO")
  
  if (nrow(meta) < 3) stop("Not enough cells in this treatment group.")
  
  # Get normalised expression values
  assay_data <- switch(
    normalisation,
    raw = fetch_count_matrix(data, log),
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
    zinb = compute_normalised_counts(data, method = "zinb", batch = batch),
    stop("Unsupported normalization method.")
  )
  
  expr <- as.numeric(assay_data[gene, meta$barcode])
  names(expr) <- meta$combined_id
  if (all(expr == 0)) stop("Gene not expressed in selected treatment.")
  
  # Assign 0 concentration for DMSO samples
  meta$concentration <- ifelse(meta$Treatment_1 == "DMSO", 0, as.numeric(meta$Concentration_1))
  
  # Data frame for modeling
  df <- data.frame(
    expression = expr,
    concentration = meta$concentration,
    replicate = meta$combined_id
  )
  
  # Fit 4-parameter logistic curve using concentration (including 0)
  model <- tryCatch({
    drm(expression ~ concentration, data = df, fct = LL.4())
  }, error = function(e) {
    warning("Sigmoidal fit failed: ", e$message)
    return(NULL)
  })  
  
  if (!is.null(model)) {
    newdata <- data.frame(concentration = seq(min(df$concentration), max(df$concentration), length.out = 100))
    newdata$predicted <- predict(model, newdata)
    
    p <- ggplot(df, aes(x = concentration, y = expression)) +
      geom_point(size = 2) +
      geom_line(data = newdata, aes(x = concentration, y = predicted), color = "blue", size = 1) +
      scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
      theme_minimal() +
      labs(
        title = paste("Sigmoidal fit for", gene, "with", treatment_value),
        x = "Concentration (log scale)",
        y = "Expression"
      )
  } else {
    p <- ggplot(df, aes(x = concentration, y = expression)) +
      geom_point(size = 2, color = "red") +
      theme_minimal() +
      labs(title = "Model failed to fit", x = "Concentration", y = "Expression")
  }
  
  p
  return(list(model = model, predictions = newdata, plot = p))
}
