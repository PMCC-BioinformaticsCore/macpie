#' Pathway enrichment analysis
#'
#' @param genes Differentially expressed genes for the enrichR analysis
#' @param db Valid name of an enrichR database
#' @param genesets List of genes
#' @param species One of "Fly", "Yeast", "Worm", "Fish", "human" or "mouse"
#'
#' @returns
#' @export
#'
#' @examples
#'
#'
#' add functionality about the species
enrichr_pathways <- function(genes = NULL, db = NULL, genesets = NULL, species = NULL) {

  #first validate the inputs
  validate_inputs <- function(genes, db, genesets) {
    if (!inherits(genes, "character")) {
      stop("Error: argument 'genes' must be a character dataset.")
    }
    if (is.null(db) && is.null(genesets)) {
      stop("Please provide either the name of the enrichr database or your own list of genes.")
    }
    if (length(db) > 0 && length(genesets) > 0) {
      stop("Please provide only the database name or the list of genes, not both.")
    }
    if (length(species) > 0 && (!species %in% c("Fly", "Yeast", "Worm", "Fish", "human", "mouse"))) {
      stop("Your value for species is not in the accepted list for enrichr. ")
    }
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
        tibble( # Convert to a tibble per list element for easier manipulation
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

  validate_inputs(genes, db, genesets)

  if(length(db) > 0){
    cat("Downloading gene sets from the Enrichr server.")
    genesets <- geneset_download(url=paste0("https://maayanlab.cloud/", species, "Enrichr/"), db)
  }

  gene_enrichment_scores <- hyper_enrich_bg(genes, genesets)

}
