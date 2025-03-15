#' Calculate clusters for umap based on DE analysis
#'
#' @param data tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param k Number of nearest neighbors in buildSNNGraph used to construct a
#'   shared nearest-neighbor (SNN) graph
#' @importFrom scran buildSNNGraph

#' @returns A tidyseurat object with cluster information in the metadata slot
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' mac <- compute_de_umap(mac)
#' mac <- find_clusters_de_umap(mac)
#' plot_de_umap(mac, group_by = "cluster", max_overlaps = 10)
find_clusters_de_umap <- function(data = NULL, k = 10) {

  validate_inputs <- function(df_umap_data, k) {
    if (!inherits(df_umap_data, "data.frame") && length(df_umap_data) > 0) {
      stop("Error: argument 'data' must contain a data frame of UMAP coordinates in the slot tool.")
    }
    if (!inherits(k, "numeric")) {
      stop("Error: argument 'k' must contain be numeric.")
    }
  }

  #fetch data from the tidyseurat object
  df_umap_data <- data@tools$umap_de

  validate_inputs(df_umap_data, k)

  ## Set seed
  set.seed(1)
  ## Set number of Nearest-Neighbours (NNs)
  ## Build the k-NN graph
  g <- buildSNNGraph(t(df_umap_data[, 2:3]), k = k)
  ## Run walktrap clustering
  g_walk <- igraph::cluster_walktrap(g)
  ## Get the cluster labels
  clus <- g_walk$membership

  # Add cluster information to the UMAP data frame
  df_umap_data$cluster <- as.character(clus)

  data@meta.data <- data@meta.data %>% select(-any_of(starts_with(c("cluster",
                                                                    "UMAP_1",
                                                                    "UMAP_2"))))
  data@meta.data <- data@meta.data %>%
    left_join(., df_umap_data, join_by("combined_id"))

  return(data)
}
