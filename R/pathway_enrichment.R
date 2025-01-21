#' Pathway enrichment analysis
#'
#' pathway_enrichment performs hypergeometric test of enrichment of a set of genes
#' in user-provided genesets of enricher databases
#' @param genes Differentially expressed genes for the enrichment analysis
#' @param db Valid name of an enrichR database
#' @param genesets Named list of genes
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @importFrom httr2 request req_perform resp_body_string
#' @importFrom purrr map list_rbind
#' @importFrom dplyr filter distinct
#' @returns Data frame of pathway-enrichment statistics
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' treatment_samples="Staurosporine_0.1"
#' control_samples<-"DMSO_0"
#' top_table <- differential_expression(mac, treatment_samples, control_samples, method = "limma_voom")
#' top_genes <- top_genes<-top_table %>%
#'   filter(p_value_adj<0.01) %>%
#'     select(gene) %>%
#'     pull()
#' pathway_enrichment(genes = top_genes, db = "MSigDB_Hallmark_2020", species = "human") %>%
#'    head()
pathway_enrichment <- function(genes = NULL, db = NULL, genesets = NULL, species = NULL) {

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

  #download enrichr gene sets
  geneset_download<-function(url, geneset){
    results_table <- request(paste0(url,"geneSetLibrary?mode=text&libraryName=",geneset)) %>%
      req_perform() %>%           #perform request
      resp_body_string() %>%      #extract body of the response
      strsplit(split = "\n") %>%  #split by newline
      .[[1]] %>%                #take the first element
      map(~ {
        x <- strsplit(.x, "\t")[[1]] #split each pathway into components
        data.frame( # Convert to a tibble per list element for easier manipulation
          pathway_name = x[1],
          description = x[2],
          gene = x[-(1:2)]
        )
      }) %>%
      list_rbind() %>% #join into a data frame to
      filter(gene != "") %>% #eliminate empty genes per pathway
      distinct(pathway_name, gene)

    results_list <- split(results_table$gene, results_table$pathway_name)

    return(results_list)
  }

  validate_inputs(genes, db, genesets, species)

  if(length(db) > 0){
    cat("Downloading gene sets from the Enrichr server.\n")
    genesets <- geneset_download(url=paste0("https://maayanlab.cloud/", species, "Enrichr/"), db)
  }

  #perform enrichment analysis
  gene_enrichment_scores <- hyper_enrich_bg(genes, genesets)
  return(gene_enrichment_scores)
}
