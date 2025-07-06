utils::globalVariables(c("NES"))

#' Find similarities between expression profiles with fgsea. 
#' Mitochondrial/ribosomal genes are filtered from the analysis
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param target Value of a target expression profile in data (target) against
#'   which the others should be compared to.
#' @param geneset List of genes whose enrichment in expression profiles should
#'   be evaluated.
#' @param n_genes_profile Number of genes to take for a target profile from a 
#'   ranked list, default 200.
#' @param direction Direction of the expression of DE genes in the target
#'   profile: one of "up", "down" and "both", default "both".
#' @param minSize Minimum size of gene sets to consider (default 15).
#' @param maxSize Maximum size of gene sets to consider (default 500).
#' @param num_cores Number of cores to use for parallel processing.
#' 
#' @importFrom fgsea fgsea
#' @importFrom stats setNames
#' @importFrom tidyr drop_na
#' @returns A tidyseurat object with screen_profile data frame in slot tools.
#' @export
#'
#' @examples
#' data(mini_mac)
#' mini_mac@tools$diff_exprs <- mini_mac@tools$diff_exprs[1:2]
#' mini_mac <- compute_multi_screen_profile(mini_mac, target = "Staurosporine_10",
#' n_genes_profile = 10, direction = "up", num_cores = 1)

compute_multi_screen_profile <- function(data = NULL,
                                         target = NULL,
                                         geneset = NULL,
                                         n_genes_profile = 200,
                                         direction = "both",
                                         minSize = 15,
                                         maxSize = 500,
                                         num_cores = 1) {
  
  # Helper: validate input
  validate_inputs <- function(de_list, target, geneset, n_genes_profile, direction, num_cores) {
    if (!inherits(de_list, "list")) {
      stop("Error: argument 'data' must be a list of DE comparisons.")
    }
    if (is.null(target) && is.null(geneset)) {
      stop("Both the target and geneset parameters are empty.")
    }
    if (!is.null(target) && !is.null(geneset)) {
      stop("Both the target and geneset parameters are present. Please select only one.")
    }
  }
  
  de_list <- data@tools$diff_exprs
  validate_inputs(de_list, target, geneset, n_genes_profile, direction, num_cores)
  
  # If target is defined, construct gene set from it
  if (!is.null(target)) {
    geneset <- list(
      de_list[[target]] %>%
        filter(!grepl("mt-", .data$gene, ignore.case = TRUE)) %>%
        filter(!grepl("Rp(s|l)", .data$gene, ignore.case = TRUE)) %>%
        filter(
          direction == "both" |
            (direction == "up" & .data$log2FC > 0) |
            (direction == "down" & .data$log2FC < 0)
        ) %>%
        arrange(if (direction == "both") desc(abs(.data$log2FC)) else desc(.data$log2FC)) %>%
        slice(1:n_genes_profile) %>%
        pull(.data$gene)
    ) %>% setNames(target)
  } else {
    if(inherits(geneset, "character")){
      geneset <- list(gene_set = geneset)
    }
  }
  
  # fgsea wrapper
  fgsea_fun <- function(x) {
    ordered_genes <- x %>%
      filter(!grepl("^mt-|^MT-", .data$gene, ignore.case = TRUE)) %>%
      filter(!grepl("^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA", .data$gene, ignore.case = TRUE))
    ranks <- ordered_genes$log2FC + rnorm(nrow(ordered_genes), sd = 1e-6)
    names(ranks) <- ordered_genes$gene
    result <- fgsea(geneset, ranks, minSize = minSize, maxSize = maxSize)
    result$target_id <- unique(x$combined_id)
    return(result)
  }
  
  # Run in parallel or not
  fgsea_list <- if (.Platform$OS.type == "windows" || num_cores <= 1) {
    lapply(de_list, fgsea_fun)
  } else {
    pmclapply(de_list, fgsea_fun, mc.cores = num_cores)
  }
  
  # Merge and normalize NES if target was used
  fgsea_df <- bind_rows(fgsea_list) %>%
    mutate(NES = -.data$NES) %>%
    drop_na(NES)
  
  if (!is.null(target)) {
    target_nes <- fgsea_df %>%
      as.data.frame() %>%
      filter(target_id == target) %>%
      pull(NES)
    
    if (length(target_nes) != 1 || is.na(target_nes)) {
      stop("Could not extract a unique NES value for the target.")
    }
    
    fgsea_df <- fgsea_df %>%
      mutate(percent_target_activity = 100 * NES / target_nes)
  }
  
  data@tools[["screen_profile"]] <- fgsea_df
  return(data)
}
