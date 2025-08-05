# Download enrichr gene sets
#' Title
#'
#' @param species One of "human", "mouse", "fly", "yeast", "worm" or "fish"
#' @param db Valid name of an enrichR database
#'
#' @returns List of genes per geneset
#' @export
#'
#' @examples
#' genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
#' head(genesets[["Adipogenesis"]])

download_geneset <- function(species = "human", db = "MSigDB_Hallmark_2020") {
  req_pkgs <- c("httr2", "purrr")
  missing <- req_pkgs[!vapply(req_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "download_geneset(): the following packages are required but not installed: ",
      paste(missing, collapse = ", "),
      "\nPlease install via `install.packages()`."
    )
  }

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

  if(!species %in% c("human", "mouse")){
    url <- paste0("https://maayanlab.cloud/", species, "Enrichr/")
  } else {
    url <- paste0("https://maayanlab.cloud/Enrichr/")
  }
  tryCatch({
    results_table <- httr2::request(paste0(url, "geneSetLibrary?mode=text&libraryName=", db)) %>%
      httr2::req_perform() %>%           #perform request
      httr2::resp_body_string() %>%      #extract body of the response
      strsplit(split = "\n") %>%  #split by newline
      .[[1]] %>%                #take the first element
      purrr::map(~ {
        x <- strsplit(.x, "\t")[[1]] #split each pathway into components
        data.frame( # Convert to a tibble per list element for easier manipulation
          pathway_name = x[1],
          description = x[2],
          gene = x[-seq_len(2)]
        )
      }) %>%
      purrr::list_rbind() %>% #join into a data frame to
      filter(.data$gene != "") %>% #eliminate empty genes per pathway
      distinct(.data$pathway_name, .data$gene)
  }, error = function(e) {
    stop("Problems connecting to enrichR database, check connection ", e$message)
  })

  results_list <- split(results_table$gene, results_table$pathway_name)

  return(results_list)
}
