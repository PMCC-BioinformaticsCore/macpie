
#' Find similarities between expression profiles with fgsea
#'
#' @param data List of expression profiles, output of the multi_DE function
#' @param target Value of a target expression profile in data (target) against
#'   which the others should be compared to
#' @param geneset List of genes whose enrichment in expression profiles should
#'   be evaluated
#' @param n_genes_profile Number of genes for a target profile, default 200
#' @param direction Direction of the expression of DE genes in the target
#'   profile: on of "up", "down" and "both", default both.
#' @param num_cores Number of cores
#' @importFrom fgsea fgsea
#' @importFrom stats setNames
#'
#' @returns Data with fgsea enrichment score
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/de_screen.Rds", package = "macpie")
#' de_list <- readRDS(file_path)
#' fgsea_results <- screen_profile(de_list, target = "Staurosporine_10",
#' n_genes_profile = 100, direction = "up", num_cores = 2)

screen_profile <- function(data = NULL,
                           target = NULL,
                           geneset = NULL,
                           n_genes_profile = NULL,
                           direction = NULL,
                           num_cores = NULL) {

  # Helper function to validate input data
  validate_inputs <- function(data, target, geneset, n_genes_profile, direction, num_cores) {
    if (!inherits(data, "list")) {
      stop("Error: argument 'data' must be a list of DE comparisons vs DMSO.")
    }
    target <- if (is.null(target)) NULL else target
    geneset <- if (is.null(geneset)) NULL else geneset
    n_genes_profile <- if (is.null(n_genes_profile)) 200 else n_genes_profile
    num_cores <- if (is.null(num_cores)) (detectCores() - 1) else num_cores

    direction <- if (is.null(direction)) "both" else direction
    if (is.null(target) && is.null(geneset)) {
      stop("Both the target and geneset parameters are empty.")
    }
    if (!is.null(target) && !is.null(geneset)) {
      stop("Both the target and geneset parameters are present, please select only one")
    }
  }

  # Validate inputs
  validate_inputs(data, target, geneset, n_genes_profile, direction, num_cores)

  #### prepare your gene list
  if (!is.null(target)) {
    geneset <- list(
      data[[target]] %>%
        filter(!grepl("mt-", .data$gene, ignore.case = TRUE)) %>%
        filter(!grepl("Rp(s|l)", .data$gene, ignore.case = TRUE)) %>%
        filter(
          direction == "both" | (direction == "up" & .data$log2FC > 0) | (direction == "down" & .data$log2FC < 0)
        ) %>%
        arrange(if (direction == "both") desc(abs(.data$log2FC)) else desc(.data$log2FC)) %>%
        slice(1:n_genes_profile) %>%
        pull(.data$gene)
    ) %>%
      setNames(target)
  }

  #perform the enrichment analysis
  fgsea_list <- list()

  fgsea_list <- pmclapply(data, function(x) {
    ordered_genes <- x %>%
      filter(!grepl("mt-", .data$gene, ignore.case = TRUE)) %>%
      filter(!grepl("Rp(s|l)", .data$gene, ignore.case = TRUE))
    ranks <- ordered_genes$log2FC
    names(ranks) <- ordered_genes$gene
    result <- fgsea(geneset, ranks, minSize = 15, maxSize = 500)
    result$target <- unique(x$target)
    return(result)
  }, mc.cores = num_cores)

  fgsea_df <- do.call("rbind", fgsea_list) %>%
    mutate(NES = -.data$NES) %>%
    filter(!is.na(.data$NES))

  return(fgsea_df)

}
