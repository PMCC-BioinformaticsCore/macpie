utils::globalVariables(c("enrichment"))
#' Perform enrichR-style analysis on a screen
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param genesets Named list of genes
#' @param direction Direction of the differentially expressed genes, 
#'  one of "up", "down", "both" (default).
#' @param p_value_cutoff Cutoff for adjusted p-value (column p_value_adj), default 0.01
#' @param n_distinct Minimum number of genes in a geneset, default 5
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @importFrom dplyr reframe
#' @importFrom tidyr unnest
#' @importFrom rlang .env
#'
#' @returns A tidyseurat object with appended pathway_enrichment dataframe in slot tools
#' @export
#'
#' @examples
#' data(mini_mac)
#' data(genesets)
#' compute_multi_enrichr(mini_mac, genesets = genesets)
compute_multi_enrichr <- function(data,
                                  genesets = NULL,
                                  species = NULL,
                                  direction = "both",
                                  p_value_cutoff = 0.01,
                                  n_distinct = 10) {

  # Helper function to validate input data
  validate_inputs <- function(data, genesets, species, direction, p_value_cutoff, n_distinct) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (!inherits(genesets, "list")) {
      stop("Error: argument 'genesets' must be a named list of genes.")
    }
    if (!is.null(species) && (!species %in% c("human", "mouse", "fly", "yeast", "worm", "fish"))) {
      stop("Your values for species should be in: fly, yeast, worm, fish, human, mouse. ")
    }
    if (!is.null(direction) && !direction %in% c("up", "down", "both")) {
      stop("Value of the direction paramater should be either up, down or both.")
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

  species <- tolower(unique(data$Species))
  validate_inputs(data, genesets, species, direction, p_value_cutoff, n_distinct)

  #modify the analysis if the species is mouse
  if(species == "mouse"){
    genesets <- lapply(genesets, convert_human_to_mouse)
  }
  
  #extract de information from the object data
  de_df <- bind_rows(data@tools$diff_exprs)

  enriched_pathways <- de_df %>%
    filter(.data$p_value_adj < .env$p_value_cutoff) %>% #select DE genes based on FDR < 0.01
    filter(
      case_when(
        .env$direction == "up" ~ .data$log2FC > 0,
        .env$direction == "down" ~ .data$log2FC < 0,
        .env$direction == "both" ~ TRUE, # Keep all genes
        TRUE ~ FALSE # Default case, in case `direction` is misspecified
      )
    ) %>%
    group_by(.data$combined_id) %>%
    filter(n_distinct(.data$gene) > .env$n_distinct) %>% #filter out samples with less than 5 DE genes
    reframe(enrichment = compute_hyper_enrich_bg(.data$gene, 
                                                 genesets = .env$genesets, 
                                                 background = .env$species)) %>%
    unnest(enrichment)
  data@tools[["pathway_enrichment"]] <- enriched_pathways
  
  return(data)
}
