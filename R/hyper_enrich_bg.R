#' Calculate enrichment of DE genes in a gene set
#' Internal function
#' @param deg vector of differentially expressed genes
#' @param genesets list of genes per pathway from enrichr
#' @param background defaults to the number of genes in the geneset, otherwise a
#'   number of genes (integer) or species
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @returns enricment stats
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' treatment_samples="Staurosporine_0.1"
#' control_samples<-"DMSO_0"
#' top_table <- differential_expression(mac, treatment_samples, control_samples,
#' method = "limma_voom")
#' top_genes <- top_genes<-top_table %>%
#'   filter(.data$p_value_adj<0.01) %>%
#'     select(gene) %>%
#'     pull()
#' genesets <- geneset_download(url=paste0("https://maayanlab.cloud/", NULL,
#' "Enrichr/"), "MSigDB_Hallmark_2020")
#' hyper_enrich_bg(top_genes, genesets) %>% head()

hyper_enrich_bg <-function(deg = NULL, #vector of DEGs
                          genesets = NULL, #pathway from enrichr
                          background = "human" #define the number of genes in your experiment (interger), or "human" to use human genome or
                          #"geneset" to use number genes in the gene set collection
) {
  #checking input
  if (!is(deg, "vector")) {
    stop("DEGs expected to be a vector of gene symbols\n")
  }
  if (!is(genesets, "list")) {
    stop("Genesets expected to be a list of non-null names\n")
  }

  if (typeof(background) != "integer") {
    #set background/universe gene number
    background <- switch(background,
                         "human" = as.character(20000),
                         #as enrichr set
                         "geneset" = length(unique(unlist(genesets))),)
    background <- as.numeric(background)

  } else{
    background = background
  }

  deg <- unique(deg)
  genesets <- lapply(genesets, unique)
  deg_genesets <- deg[deg %in% unique(unlist(genesets))]
  n_hits <- sapply(genesets, function(x, y)
    length(intersect(x, y)), deg_genesets) #q
  n_hits_updatebg <- n_hits != 0 #updating background
  genesets <- genesets[n_hits_updatebg] #updating background
  n_hits <- n_hits[n_hits > 0] #exlude 0 overlapping terms


  n_genesets <- sapply(genesets, length)#m
  n_minus <- background - n_genesets#n
  n_deg <- length(deg)#k

  #hyper-geometric
  pvals <- stats::phyper(
    q = n_hits - 1,
    m = n_genesets,
    n = n_minus,
    k = n_deg,
    lower.tail = FALSE
  )

  #put into a dataframe
  res_data <- data.frame(
    term = names(genesets),
    overlap = paste0(n_hits, "/", n_genesets),
    p_value = pvals,
    p_value_adjust = p.adjust(pvals, "BH"),
    genes = sapply(genesets, function(x, y)
      paste(intersect(x, y), collapse = ";"), deg_genesets)
  )

  res_data <- arrange(res_data, p_value)
  return(results = res_data)
}
