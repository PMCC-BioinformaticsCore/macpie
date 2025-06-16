build_mofa <- function(data, 
                             combined_ids = NULL, 
                             metadata_col = "Treatment_1",
                             metric = "metric", 
                             control_treatment = "DMSO",
                             pval_thresh = 0.01,
                             pathway_pval_thresh = 0.05,
                             padj_col = "p_value_adj",
                             pathway_pval_col = "Adjusted.P.value",
                             score_col = "Combined.Score") {
  
  meta <- data@meta.data
  tools <- data@tools
  
  # If no combined IDs provided, default to highest concentration or DMSO treatments
  if (is.null(combined_ids)) {
    combined_ids <- data %>%
      filter("Concentration_1" == max(as.numeric(as.character(meta$Concentration_1)))) %>%
      pull(combined_id) %>% unique()
  }
  
  # Pathway Enrichment 
  pathway_mat <- tools$pathway_enrichment %>%
    filter(combined_id %in% combined_ids, !!sym(pathway_pval_col) < pathway_pval_thresh) %>%
    dplyr::select("combined_id", "Term", "Combined.Score") %>%
    pivot_wider(names_from = "Term", values_from = !!sym(score_col)) %>%
    left_join(meta %>% dplyr::select(combined_id, !!sym(metadata_col)) %>% distinct(), by = "combined_id") %>%
    relocate(!!sym(metadata_col)) %>%
    column_to_rownames(metadata_col) %>%
    dplyr::select(-combined_id) %>%
    t()
  
  # Gene Expression
  all_de <- bind_rows(tools$diff_exprs, .id = "comparison")
  signif_genes <- all_de %>%
    filter(!!sym(padj_col) < pval_thresh) %>%
    pull(gene) %>% unique()
  
  gene_mat <- all_de %>%
    filter(gene %in% signif_genes, combined_id %in% combined_ids) %>%
    dplyr::select(combined_id, gene, !!sym(metric)) %>%
    pivot_wider(names_from = gene, values_from = !!sym(metric)) %>%
    left_join(meta %>% dplyr::select(combined_id, !!sym(metadata_col)) %>% distinct(), by = "combined_id") %>%
    relocate(!!sym(metadata_col)) %>%
    column_to_rownames(metadata_col) %>%
    dplyr::select(-combined_id) %>%
    t()
  
  # Chemical descriptors
  desc_mat <- tools$chem_descriptors %>%
    dplyr::select(-"clean_compound_name") %>%
    column_to_rownames(metadata_col) %>%
    t()
  
  # Additional descriptors
  phenotype_mat <- data@meta.data %>%
    filter(combined_id %in% combined_ids) %>%
    dplyr::select(combined_id, all_of("additional_features")) %>%
    group_by(combined_id) %>%
    summarise(across(all_of("additional_features"), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
    left_join(
      data@meta.data %>% dplyr::select("combined_id", "Treatment_1") %>% distinct(),
      by = "combined_id"
    ) %>%
    relocate("Treatment_1") %>%
    dplyr::select(-combined_id) %>%
    column_to_rownames("Treatment_1") %>%
    t()
  
  # Match columns
  common_cols <- Reduce(intersect, list(
    colnames(gene_mat),
    colnames(pathway_mat),
    colnames(desc_mat),
    colnames(phenotype_mat)
  ))
  
  phenotype_views <- lapply(rownames(phenotype_mat), function(feature) {
    mat <- phenotype_mat[feature, common_cols, drop = FALSE]
    colnames(mat) <- common_cols
    mat
  })
  names(phenotype_views) <- rownames(phenotype_mat)  # e.g. "nCount_RNA", "CTG", etc.
  
  # Combine all views
  views <- c(
    phenotype_views,
    list(
      genes = gene_mat[, common_cols, drop = FALSE],
      pathways = pathway_mat[, common_cols, drop = FALSE],
      chem_descriptors = desc_mat[, common_cols, drop = FALSE]
    )
  )
  
  # Ensure all views have matching column names
  views <- lapply(views, function(view) {
    colnames(view) <- common_cols
    return(view)
  })
  
  # Create MOFA object
  mofa_obj <- MOFA2::create_mofa(views)
  return(mofa_obj)
}
