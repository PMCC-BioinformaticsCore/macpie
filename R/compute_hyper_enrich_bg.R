#' Calculate enrichment of DE genes in a gene set
#' Internal function
#' @param deg vector of differentially expressed genes
#' @param genesets list of genes per pathway from enrichr
#' @param background defaults to the number of genes in the geneset, otherwise a
#'   number of genes (integer) or species
#' @importFrom stats p.adjust
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @returns enrichment stats
#' @export
#'
#' @examples
#' data(mini_mac)
#' treatment_samples="Staurosporine_0.1"
#' control_samples<-"DMSO_0"
#' top_table <- compute_single_de(mini_mac, treatment_samples, control_samples,
#' method = "limma_voom")
#' top_genes <- top_table$gene[top_table$p_value_adj<0.01]
#' data(genesets)
#' results <- compute_hyper_enrich_bg(top_genes, genesets)
#' head(results)

compute_hyper_enrich_bg <- function(deg = NULL, #vector of DEGs
  genesets = NULL, #pathway from enrichr
  background = "human"
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
                         "mouse" = as.character(22000),
                         #as enrichr set
                         "geneset" = length(unique(unlist(genesets))))
    background <- as.numeric(background)

  } else {
    background <- background
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

  # Calculate z-scores
  expected_hits <- n_genesets * (n_deg / background)
  std_deviation_hits <- sqrt(n_genesets * (n_deg / background) * (1 - n_deg / background))
  z_scores <- (n_hits - expected_hits) / std_deviation_hits

  # Calculate Combined Scores
  combined_scores <- -log(pvals) * z_scores

  #put into a dataframe
  res_data <- data.frame(
    Term = names(genesets),
    Overlap = paste0(n_hits, "/", n_genesets),
    P.value = pvals,
    Adjusted.P.value = p.adjust(pvals, "BH"),
    Genes = sapply(genesets, function(x, y) {
      paste(intersect(x, y), collapse = ";")
    }, deg_genesets),
    Combined.Score = combined_scores
  )

  res_data <- res_data %>%
    arrange(.data$P.value)
  return(results = res_data)
}
