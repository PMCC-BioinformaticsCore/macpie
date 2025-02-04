
#' Perform enrichR-style analysis on a screen
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param genesets Named list of genes
#' @param p_value_cutoff Cutoff for adjusted p-value (column p_value_adj), default 0.01
#' @param n_distinct Minimum number of genes in a geneset, default 5
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @importFrom dplyr reframe
#'
#' @returns A tidyseurat object with pathway_enrichment list in slot tools
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' file_path <- system.file("extdata", "PMMSq033/pathways.Rds", package = "macpie")
#' genesets <- readRDS(file_path)
#' multi_enrich_pathways(mac, genesets = genesets)
multi_enrich_pathways <- function(data,
                                  genesets = enrichr_genesets,
                                  species = NULL,
                                  p_value_cutoff = 0.01,
                                  n_distinct = 5) {

  # Helper function to validate input data
  validate_inputs <- function(data, genesets, species, p_value_cutoff, n_distinct) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (!inherits(genesets, "list")) {
      stop("Error: argument 'genesets' must be a named list of genes.")
    }
    if (!is.null(species) && (!species %in% c("human", "mouse", "fly", "yeast", "worm", "fish"))) {
      stop("Your values for species should be in: fly, yeast, worm, fish, human, mouse. ")
    }
    if (!inherits(p_value_cutoff, "numeric")) {
      stop("Error: argument 'p_value_cutoff' must be numeric.")
    }
    if (!inherits(n_distinct, "numeric")) {
      stop("Error: argument 'n_distinct' must be numeric.")
    }
    if (is.null(species) || length(species) > 1) {
      stop("Missing species information.")
    }
    if (length(data@tools$diff_exprs) == 0) {
      stop("Missing information on DE genes. Run multi_DE first.")
    }
  }

  species <- unique(mac$Species)
  validate_inputs(data, genesets, species, p_value_cutoff, n_distinct)

  #extract de information from the object data
  de_df <- bind_rows(data@tools$diff_exprs)

  enriched_pathways <- de_df %>%
    filter(.data$p_value_adj < .env$p_value_cutoff) %>% #select DE genes based on FDR < 0.01
    group_by(.data$combined_id) %>%
    filter(n_distinct(.data$gene) > .env$n_distinct) %>% #filter out samples with less than 5 DE genes
    reframe(enrichment = hyper_enrich_bg(.data$gene, genesets = .env$genesets, background = .env$species)) %>%
    unnest(enrichment)
  mac@tools[["pathway_enrichment"]] <- enriched_pathways
  return(mac)
}
