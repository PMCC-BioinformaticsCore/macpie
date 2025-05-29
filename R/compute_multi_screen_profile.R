utils::globalVariables(c("NES"))
#' Find similarities between expression profiles with fgsea
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param target Value of a target expression profile in data (target) against
#'   which the others should be compared to
#' @param geneset List of genes whose enrichment in expression profiles should
#'   be evaluated
#' @param n_genes_profile Number of genes for a target profile, default 200
#' @param direction Direction of the expression of DE genes in the target
#'   profile: one of "up", "down" and "both", default both.
#' @param num_cores Number of cores
#' @importFrom fgsea fgsea
#' @importFrom stats setNames
#' @importFrom tidyr drop_na
#'
#' @returns A tidyseurat object with screen_profile data frame in slot tools
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac@tools$diff_exprs<-mini_mac@tools$diff_exprs[1:2]
#' mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10",
#' n_genes_profile = 10, direction = "up", num_cores = 1)

compute_multi_screen_profile <- function(data = NULL,
                                         target = NULL,
                                         geneset = NULL,
                                         n_genes_profile = 200,
                                         direction = "both",
                                         num_cores = parallel::detectCores() - 1) {

  # Helper function to validate input data
  validate_inputs <- function(de_list, target, geneset, n_genes_profile, direction, num_cores) {
    if (!inherits(de_list, "list")) {
      stop("Error: argument 'data' must be a list of DE comparisons.")
    }
    target <- if (is.null(target)) NULL else target
    geneset <- if (is.null(geneset)) NULL else geneset
    if (is.null(target) && is.null(geneset)) {
      stop("Both the target and geneset parameters are empty.")
    }
    if (!is.null(target) && !is.null(geneset)) {
      stop("Both the target and geneset parameters are present, please select only one")
    }
  }

  de_list <- data@tools$diff_exprs

  # Validate inputs
  validate_inputs(de_list, target, geneset, n_genes_profile, direction, num_cores)

  #### prepare your gene list
  if (!is.null(target)) {
    geneset <- list(
      de_list[[target]] %>%
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

  #define when it's no support for mclapply
  fgsea_fun <- function(x) {
    ordered_genes <- x %>%
      filter(!grepl("mt-", .data$gene, ignore.case = TRUE)) %>%
      filter(!grepl("Rp(s|l)", .data$gene, ignore.case = TRUE))
    ranks <- ordered_genes$log2FC
    names(ranks) <- ordered_genes$gene
    result <- fgsea(geneset, ranks, minSize = 15, maxSize = 500)
    result$target <- unique(x$combined_id)
    return(result)
  }
  
  #if it's on Windows, use lapply
  if (.Platform$OS.type== "windows" || num_cores <= 1) {
    fgsea_list <- lapply(de_list, fgsea_fun)
  } else {
      fgsea_list <- pmclapply(de_list, fgsea_fun, mc.cores= num_cores)
    }
  

  #reorder the NES direction (why though?)
  fgsea_df <- bind_rows(fgsea_list) %>%
    mutate(NES = -.data$NES) %>%
    drop_na(NES)

  data@tools[["screen_profile"]] <- fgsea_df

  return(data)
}
