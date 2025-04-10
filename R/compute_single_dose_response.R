#' Model Gene Dose-Response Curve Using drc
#'
#' @param data A Seurat or TidySeurat object containing expression and metadata.
#' @param gene A gene name (must match a row name in the object).
#' @param normalisation One of "raw", "logNorm", "cpm", "clr", "SCT", "DESeq2",
#'   "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom", "zinb". If empty, defaults to cpm
#' @param treatment_value A character string matching one value in metadata column "Treatment_1".
#' @import drc
#' @import tidyverse
#' @importFrom scales pseudo_log_trans
#' @return A list with drc model, predicted values, and ggplot curve
#' @examples
#' rds_file<-system.file("/extdata/PMMSq033/PMMSq033.rds", package = "macpie")
#' mac<-readRDS(rds_file)
#' res <- compute_single_dose_response(data = mac,
#' gene = "SOX12",
#' normalisation = "limma_voom",
#' treatment_value = "Staurosporine")
#' res$plot
#' res <- compute_single_dose_response(data = mac,
#' pathway = "Myc Targets V1",
#' treatment_value = "Camptothecin")
#' res$plot
#' @export
compute_single_dose_response <- function(data,
                                         gene = NULL,
                                         pathway = NULL,
                                         normalisation = "limma_voom",
                                         treatment_value,
                                         control_value = "DMSO",
                                         batch = 1,
                                         k = 2) {
  # Helper function to fetch and log-transform count matrix
  fetch_count_matrix <- function(data, log) {
    count_matrix <- as.matrix(data@assays$RNA$counts)
    if (log) {
      count_matrix <- log1p(count_matrix)
    }
    return(count_matrix)
  }
  # Helper function to validate input data
  validate_inputs <- function(data, gene, pathway, normalisation, treatment_value, control_value, batch, k) {
    if (!inherits(data, "Seurat")) {
      stop("Error: 'data' must be a Seurat or TidySeurat object.")
    }
    if (!is.null(gene)) {
      if (!gene %in% row.names(data@assays$RNA$counts)) {
        stop("Error: Your gene is not present in the dataset.")
      }
    }
    if (!is.null(pathway)) {
      if (!pathway %in% data@tools$pathway_enrichment$Term) {
        stop("Error: Your pathway was not present in the list 
             of enriched pathways. Check mac@tools$pathway_enrichment.")
      }
    }
    if (!is.null(gene) && !is.null(pathway)) {
      stop("Error: Please select only gene OR pathway, not both.")
    }
    if (!treatment_value %in% data$Treatment_1) {
      stop("Error: Your treatment_value was not present in data$Treatment_1.")
    }
    if (!control_value %in% data$Treatment_1) {
      stop("Error: Your control_value was not present in data$Treatment_1.")
    }
    normalisation <- if (is.null(normalisation)) "limma_voom" else normalisation
    if (!normalisation %in% c("raw", "logNorm",
                              "cpm", "clr", "SCT",
                              "DESeq2", "edgeR",
                              "RUVg", "RUVs", "RUVr",
                              "limma_voom", "zinb")) {
      stop("Your normalization method is not available.")
    }
    batch <- if (is.null(batch)) "1" else as.character(batch)
    k <- if (is.null(k)) 2 else k
  }
  validate_inputs(data, gene, pathway, normalisation, treatment_value, control_value, batch, k)
  # Subset metadata
  data <- data %>%
    filter(.data$Treatment_1 == treatment_value | .data$Treatment_1 == control_value)
  meta <- data@meta.data %>%
    mutate(barcode = rownames(.)) %>%
    filter(.data$Treatment_1 == treatment_value | .data$Treatment_1 == control_value)
  if (nrow(meta) < 3) stop("Not enough cells in this treatment group.")
  if (!is.null(gene)) {
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
    # Assign 0 concentration for control_value samples
    meta$concentration <- ifelse(meta$Treatment_1 == control_value, 0, as.numeric(as.character(meta$Concentration_1)))
    # Data frame for modeling
    df <- data.frame(
      expression = expr,
      concentration = meta$concentration,
      replicate = meta$combined_id
    )
  } else if (!is.null(pathway)) {
    pathway_enrichment <- data@tools$pathway_enrichment %>%
      filter(.data$Term == .env$pathway) %>%
      select(.data$combined_id, .data$Combined.Score)
    meta <- data@meta.data %>%
      select(.data$Concentration_1, .data$combined_id, .data$Treatment_1) %>%
      unique() %>%
      left_join(., pathway_enrichment, join_by(combined_id)) %>%
      filter(grepl(.env$treatment_value, combined_id)) %>%
      unique()
    meta$Combined.Score[is.na(meta$Combined.Score)] <- 0
    expr <- expr_pathway$Combined.Score
    names(expr) <- expr_pathway$combined_id
    df <- data.frame(
      expression = meta$Combined.Score,
      concentration = as.numeric(as.character(meta$Concentration_1))  # Convert from factor to numeric
    ) %>%
      bind_rows(data.frame(expression = 0, concentration = 0))
  }
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
    ed <- ED(model, 50, interval = "delta")
    ec50 <- ed[1, "Estimate"]
    # Dynamically determine y-axis label
    if (is.null(gene) && !is.null(pathway)) {
      y_label <- "Enrichment score"
      marker <- pathway
    } else if (!is.null(gene)) {
      y_label <- paste0("Expression (", normalisation, ")")
      marker <- gene
    } else {
      y_label <- "Metric"
      marker <- column_name
    }
    p <- ggplot(df, aes(x = concentration, y = expression)) +
      geom_point(size = 2) +
      geom_line(data = newdata, aes(x = concentration, y = predicted), color = "blue", size = 1) +
      scale_x_continuous(
        trans = pseudo_log_trans(base = 10),
        breaks = sort(unique(df$concentration)),
        labels = function(x) {
          sapply(x, function(val) {
            if (is.na(val)) {
              NA
            } else if (val < 1) {
              val
            } else {
              as.character(round(val))
            }
          })
        }
      ) +
      geom_vline(xintercept = ec50, linetype = "dashed", color = macpie_colours$high) +
      annotate("text",
               x = ec50,
               y = max(df$expression, na.rm = TRUE) * 0.95,
               label = sprintf("EC50 = %.2f\n", ec50),
               color = "black",
               size = 4,
               hjust = 0.5,
               vjust = 1) +
      theme_minimal() +
      labs(
        title = paste("Sigmoidal fit for", marker, "with", treatment_value),
        x = "Concentration (log scale)",
        y = y_label
      )
  } else {
    p <- ggplot(df, aes(x = concentration, y = expression)) +
      geom_point(size = 2, color = macpie_colours$high) +
      theme_minimal() +
      labs(
        title = "Model failed to fit",
        x = "Concentration",
        y = y_label
      )
  }
  p
  return(list(model = model, predictions = newdata, plot = p))
}
