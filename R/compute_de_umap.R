
#' Prepare DE-based umap
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @importFrom dplyr bind_rows select
#' @importFrom tidyr pivot_wider
#' @importFrom umap umap
#' @returns A tidyseurat object with umap_de data frame in slot tools
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' mac <- compute_de_umap(mac)

compute_de_umap <- function(data = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(de_list) {
    if (!inherits(de_list, "list") && length(de_list) > 0) {
      stop("Error: argument 'data' must contain a list of DE comparisons in the slot tool.")
    }
  }

  de_list <- data@tools$diff_exprs
  validate_inputs(de_list)

  df <- bind_rows(de_list)
  df_wide <- df %>%
    select("gene", "combined_id", "metric") %>%
    pivot_wider(names_from = "combined_id", values_from = "metric")

  set.seed(1)
  df_umap <- umap(t(df_wide[, -1]))

  #plot controls on UMAP
  df_umap_data <- df_umap$layout %>%
    as.data.frame() %>%
    rownames_to_column("combined_id") %>%
    rename("UMAP_1" = "V1", "UMAP_2" = "V2")

  data@tools[["umap_de"]] <- df_umap_data
  return(data)
}
