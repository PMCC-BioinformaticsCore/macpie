#' Create a distance heatmap
#'
#' Plot a heatmap to show poisson distance on count matrix of pre-processed data.
#' The poisson distance matrix is calculated on subset of wells that specify by two
#' other parameters. By specifying group_by to a column in the metadata and a specific
#' treatment of interest, poisson distance is then calculated for corresponding wells.
#' @param data A Seurat object
#' @param group_by A metadata column name to group data
#' @param treatment to specify one treatment group in the group_by parameter
#'
#' @return it returns a pheatmap object
#' @examples
#' data(mini_mac)
#' p <- plot_distance(mini_mac, group_by = "combined_id", treatment = "DMSO_0")
#' @export

plot_distance <- function(data = NULL, group_by = NULL, treatment = NULL) {
  req_pkgs <- c("PoiClaClu", "pheatmap")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "plot_distance(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }
  validate_inputs <- function(data, group_by, treatment) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    all_treatments <- data %>% pull(group_by) %>% unique() %>% as.character()
    if (!treatment %in% all_treatments) {
      stop("Your treatment names are not present in the data or metadata.")
    }
    column_names <- data %>%
      head() %>%
      colnames()
    if (!all(c(group_by) %in% column_names)) {
      stop("Your column names are not present in the data or metadata.")
    }
    list(data = data, group_by = group_by, treatment = treatment)
  }

  validated <- validate_inputs(data, group_by, treatment)
  group_by <- validated$group_by
  treatment <- validated$treatment
  mac_val <- validated$data

  # Subset on group_by
  mac_val  <- mac_val  %>% filter(.data[[group_by]] == treatment)

  # Poisson distance
  poisd <- PoiClaClu::PoissonDistance(t(as.matrix(mac_val@assays$RNA$counts)))
  poisd_matrix <- as.matrix(poisd$dd)
  rownames(poisd_matrix) <- paste0(mac_val@meta.data[[group_by]], "_", mac_val$Well_ID)
  # Order the matrix
  well_order <- as.numeric(order(str_replace(rownames(poisd_matrix), paste0(treatment, "_"), "")))
  poisd_matrix <- poisd_matrix[well_order, ]
  poisd_matrix <- poisd_matrix[, well_order]
  colnames(poisd_matrix) <- rownames(poisd_matrix)
  poisd_matix_reordered <- poisd_matrix[rownames(poisd_matrix), ]

  p <- pheatmap::pheatmap(poisd_matix_reordered,
                clustering_distance_rows = poisd$dd,
                clustering_distance_cols = poisd$dd,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                color = macpie_colours$continuous_rev)
  return(p)
}
