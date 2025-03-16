utils::globalVariables(c(".data"))
#' Generate heatmap of DE genes from multiple treatments
#' 
#' This is the function to generate a heatmap of DE genes from running compute_multi_DE 
#' that shared by more than one treatment group. 
#'
#' @param data A tidyseurat object merged with metadata. Must contain columns
#'   "Well_ID", "Row", "Column".
#' @param p_value_cutoff Cutoff for adjusted p-value (column p_value_adj), default 0.01
#' @param direction Direction to select up or down regulated genes indicating by log2FC > 0 or < 0  
#' @param n_genes Top n genes (ordered by Log2FC) to be extracted from each treatment comparison
#' @param control The control group to be included in the final heatmap, usually DMSO_0
#' 
#' @import pheatmap
#' @import dplyr 
#'
#' @returns a heatmap by pheatmap 
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "PMMSq033/PMMSq033.rds", package = "macpie")
#' mac <- readRDS(file_path)
#' plot_multi_de(mac, p_value_cutoff = 0.01, direction="up", n_genes = 10, control = "DMSO_0")


plot_multi_de <- function(data,
                                  p_value_cutoff = 0.01,
                                  direction = "up", n_genes =10, control="DMSO_0" ) {
  
  # Helper function to validate input data
  validate_inputs <- function(data, p_value_cutoff, direction, n_genes, control) {
    if (!inherits(data, "Seurat")) {
      stop("Error: argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (!inherits(p_value_cutoff, "numeric")) {
      stop("Error: argument 'p_value_cutoff' must be numeric.")
    }
    if (!inherits(direction, "character")) {
      stop("Error: argument 'direction' must be characters either up or down.")
    }
    if (!inherits(n_genes, "numeric")) {
      stop("Error: argument 'n_genes' must be numeric.")
    }
    if (!inherits(control, "character")) {
      stop("Error: argument 'control' must be control group in same format as combined_id.")
    }
    
    if (length(data@tools$diff_exprs) == 0) {
      stop("Missing information on DE genes. Run multi_DE first.")
    }
  }
  

  validate_inputs(data, p_value_cutoff, direction, n_genes, control)
  
  
  #extract de information from the object data
  all_de <- mac@tools$diff_exprs
  de_df <- bind_rows(all_de)
  
  #filter DE on adjusted p value and log fold-change direction
  filtered_de_df <- de_df %>%
    filter(.data$p_value_adj < .env$p_value_cutoff) %>%
    filter(if(direction=="up") .data$log2FC > 0 else .data$log2FC <0)
  
  
  top_genes_per_combined_id <- filtered_de_df %>%
    group_by(.data$combined_id) %>%
    slice_max(order_by = .data$log2FC, n = n_genes)
  
  common_genes <- top_genes_per_combined_id %>%
    group_by(.data$gene) %>%
    filter(n_distinct(.data$combined_id) > 1)
  
  features <- common_genes$gene
  
  common_genes_treatments <- common_genes %>% dplyr::distinct(combined_id) %>%pull()%>%unique()
  #add dmso
  common_genes_treatments <- c(common_genes_treatments, control)
    
  #calculate cpm on full samples
  lcpm <- cpm(mac@assays$RNA$counts, log = T)
  combined_id_barcodes <- mac@meta.data$combined_id
  names(combined_id_barcodes)<-rownames(mac@meta.data)
  colnames(lcpm) <- combined_id_barcodes[colnames(lcpm)]
  
  sub_lcpm <- lcpm[features, colnames(lcpm) %in% common_genes_treatments]
    

  
  col_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(common_genes_treatments))
  p <- pheatmap(sub_lcpm, 
           cexRow=0.3,
           cexCol=0.4)
  
  return(p)
  
  

}


