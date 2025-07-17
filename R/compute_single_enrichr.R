#' Pathway enrichment analysis
#'
#' pathway_enrichment performs hypergeometric test of enrichment of a set of genes
#' in user-provided genesets of enricher databases
#' @param genes Differentially expressed genes for the enrichment analysis
#' @param db Valid name of an enrichR database
#' @param genesets Named list of genes
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @importFrom dplyr filter distinct
#' @returns Data frame of pathway-enrichment statistics
#' @export
#'
#' @examples
#' data(mini_mac)
#' top_genes <- mini_mac@tools$diff_exprs[["Staurosporine_10"]]$gene
#' data(genesets)
#' compute_single_enrichr(genes = top_genes, genesets = genesets, species="human")

compute_single_enrichr <- function(genes = NULL, db = NULL, genesets = NULL, species = NULL) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    stop(
      "compute_single_enrichr(): the following package is required but not installed: enrichR",
      "\nPlease install via `install.packages()`.")
  }

  #first validate the inputs
  validate_inputs <- function(genes, db, genesets, species) {
    if (!inherits(genes, "character")) {
      stop("Error: argument 'genes' must be a character dataset.")
    }
    if (is.null(db) && is.null(genesets)) {
      stop("Please provide either the name of the enrichr database or your own list of genes.")
    }
    if (!is.null(db) && !is.null(genesets)) {
      stop("Please provide only the database name or the list of genes, not both.")
    }
    if (!is.null(species) && (!species %in% c("human", "mouse", "fly", "yeast", "worm", "fish"))) {
      stop("Your values for species should be in: fly, yeast, worm, fish, human, mouse. ")
    }
    if (is.null(species)) stop("Missing species information.")
  }

  validate_inputs(genes, db, genesets, species)

  if (length(db) > 0) {
    cat("Downloading gene sets from the Enrichr server.\n")
    genesets <- download_geneset(species, db)
  }

  #perform enrichment analysis
  gene_enrichment_scores <- compute_hyper_enrich_bg(genes, genesets)
  return(gene_enrichment_scores)
}
