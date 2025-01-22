
#' Perform DE of multiple samples in a screen vs control
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column"
#' @param treatment_samples Value in the column "combined_id" representing replicates of treatment samples in the data
#' @param control_samples Value in the column "combined_id"  representing replicates of control samples in the data
#' @param method One of "Seurat_wilcox", "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"
#' @param batch Either empty, a single value, or a vector corresponding to the
#'   number of samples
#' @param k Parameter k for RUVSeq methods, check RUVSeq tutorial
#' @param spikes List of genes to use as spike controls
#' @param num_cores Number of cores
#' @importFrom parallel detectCores
#' @importFrom mcprogress pmclapply
#' @returns List of DE counts vs control
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)[50:150]
#' control_samples<-"DMSO_0"
#' treatment_samples <- mac$combined_id[!grepl("DMSO", mac$combined_id)]
#' de_list <- multi_DE(mac, treatment_samples, control_samples, num_cores = 2, method = "edgeR")

multi_DE <- function(data = NULL,
                     treatment_samples = NULL,
                     control_samples = NULL,
                     method = "edgeR",
                     num_cores = 2,
                     batch = 1,
                     k = 2,
                     spikes = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, treatment_samples, control_samples, method, num_cores) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (!"combined_id" %in% colnames(data@meta.data)) {
      data <- data %>%
        mutate(combined_id = apply(select(starts_with("Treatment_") | starts_with("Concentration_")),
                                   1, paste, collapse = "_")) %>%
        mutate(combined_id = gsub(" ", "", .data$combined_id))
    }
    method <- if (is.null(method)) "edgeR" else method
    if (!method %in% c("Seurat_wilcox", "DESeq2", "edgeR",
                       "RUVg", "RUVs", "RUVr",
                       "limma_voom")) {
      stop("Your normalization method is not available.")
    }
    if (is.null(control_samples)) {
      stop("Missing control samples.")
    }
    if (is.null(treatment_samples)) {
      warn("Missing treatment samples, using all that differ from control.")
      treatment_samples <- data %>%
        select(.data$combined_id) %>%
        filter(!grepl(control_samples, combined_id)) %>%
        pull() %>%
        unique()
    }
    if (length(treatment_samples) == 1 && length(control_samples) == 1) {
      treatment_samples_list <- grepl(treatment_samples, data$combined_id)
      control_samples_list <- grepl(control_samples, data$combined_id)
      if (any(sum(treatment_samples_list) == 0, sum(control_samples_list) == 0)) {
        stop("The combined id of your samples (format: 'treatment'_'concentration') is not valid.")
      }
    }
    num_cores <- if (is.null(num_cores)) (detectCores() - 1) else num_cores
    return(list(data = data, treatment_samples = treatment_samples, num_cores = num_cores))
  }

  validated <- validate_inputs(data, treatment_samples, control_samples, method, num_cores)
  data <- validated$data
  treatment_samples <- validated$treatment_samples
  num_cores <- validated$num_cores

  de_list <- pmclapply(treatment_samples, function(x) {
    result <- differential_expression(data, x, control_samples, method, batch, k, spikes)
    result$combined_id <- x
    return(result)
  }, mc.cores = num_cores)
  return(de_list)
}
