# Download enrichr gene sets
#' Title
#'
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @param db Valid name of an enrichR database
#'
#' @importFrom httr2 request req_perform resp_body_string
#' @importFrom purrr map list_rbind
#' @returns List of genes per geneset
#' @export
#'
#' @examples
#' genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
#' head(genesets[["Adipogenesis"]])

download_geneset <- function(species = "human", db = "MSigDB_Hallmark_2020") {

  #first validate the inputs
  validate_inputs <- function(species, db) {
    if (is.null(db)) {
      stop("Please provide the name of the enrichr database.")
    }
    if (!is.null(species) && (!species %in% c("human", "mouse", "fly", "yeast", "worm", "fish"))) {
      stop("Your values for species should be in: fly, yeast, worm, fish, human, mouse. ")
    }
    if (is.null(species)) stop("Missing species information.")
  }

  validate_inputs(species, db)

  url <- paste0("https://maayanlab.cloud/", species, "Enrichr/")

  results_table <- request(paste0(url, "geneSetLibrary?mode=text&libraryName=", db)) %>%
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
    filter(.data$gene != "") %>% #eliminate empty genes per pathway
    distinct(.data$pathway_name, .data$gene)

  results_list <- split(results_table$gene, results_table$pathway_name)

  return(results_list)
}
