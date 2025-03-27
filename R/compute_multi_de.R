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
#' @import glmGamPoi
#' @import doParallel
#' @import foreach
#' @import progressr
#' @returns List of DE counts vs control
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)[50:150]
#' control_samples <- "DMSO_0"
#' treatment_samples <- mac$combined_id[!grepl("DMSO", mac$combined_id)]
#' mac <- compute_multi_de(mac, treatment_samples, control_samples, num_cores = 1, method = "edgeR")

#' Perform DE of multiple samples in a screen vs control (Windows-compatible)
#' @export
compute_multi_de <- function(data = NULL,
                             treatment_samples = NULL,
                             control_samples = NULL,
                             method = "edgeR",
                             num_cores = 2,
                             batch = 1,
                             k = 2,
                             spikes = NULL) {
  
  # Validate input
  validate_inputs <- function(data, treatment_samples, control_samples, method, num_cores) {
    if (!inherits(data, "Seurat")) stop("data must be a Seurat or TidySeurat object.")
    
    if (!"combined_id" %in% colnames(data@meta.data)) {
      data@meta.data$combined_id <- apply(
        dplyr::select(data@meta.data, dplyr::starts_with("Treatment_") | dplyr::starts_with("Concentration_")),
        1, paste, collapse = "_"
      ) |> gsub(" ", "", x = _)
    }
    
    if (is.null(control_samples)) stop("Missing control samples.")
    
    if (is.null(treatment_samples)) {
      message("Missing treatment samples — inferring from data.")
      treatment_samples <- data@meta.data %>%
        dplyr::filter(!grepl(control_samples, .data$combined_id)) %>%
        dplyr::pull(.data$combined_id) %>%
        unique()
    }
    
    num_cores <- if (is.null(num_cores)) max(1, parallel::detectCores() - 1) else num_cores
    method <- match.arg(method, choices = c("Seurat_wilcox", "DESeq2", "edgeR", "RUVg", "RUVs", "RUVr", "limma_voom"))
    
    list(data = data, treatment_samples = treatment_samples, num_cores = num_cores)
  }
  
  validated <- validate_inputs(data, treatment_samples, control_samples, method, num_cores)
  data <- validated$data
  treatment_samples <- validated$treatment_samples
  num_cores <- validated$num_cores
  
  # Setup parallel backend 
  future::plan(future::multisession, workers = num_cores)
  
  progressr::handlers(global = TRUE)
  progressr::handlers("progress")
  options(progressr.clear = FALSE)
  
  de_list <- suppressWarnings(
    progressr::with_progress({
      p <- progressr::progressor(steps = length(treatment_samples))
      
      furrr::future_map(
        treatment_samples,
        function(x) {
          result <- tryCatch({
            out <- compute_single_de(data, x, control_samples, method, batch, k, spikes)
            out$combined_id <- x
            out
          }, error = function(e) {
            message("Error in sample ", x, ": ", e$message)
            tibble(log2FC = NA_real_, metric = NA_real_, p_value = NA_real_, p_value_adj = NA_real_, combined_id = x)
          })
          
          p()  # progress tick
          result
        },
        .options = furrr::furrr_options(seed = TRUE)
      )
    })
  )
  
  names(de_list) <- treatment_samples
  data@tools[["diff_exprs"]] <- de_list
  return(data)
}

