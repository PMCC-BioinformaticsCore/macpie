build_mofa <- function(data, 
                             combined_ids = NULL, 
                             treatment_ID = "Compound_ID",
                             metadata_columns = NULL,
                             de_metric = "metric", 
                             de_pval_thresh = 0.01,
                             de_padj_col = "p_value_adj",
                             pathway_pval_thresh = 0.05,
                             pathway_pval_col = "Adjusted.P.value",
                             pathway_score_col = "Combined.Score") {
  
  meta <- data@meta.data
  tools <- data@tools
  
  # If no combined IDs provided, default to highest concentration or DMSO treatments
  if (is.null(combined_ids)) {
    stop("Please provide a list of combined_ids, these should correspond to your DE comparisons.")
  }
  
  if(treatment_ID == "combined_id"){
    cat("Warning: your parameter treatment_ID should not point to a column that contains 
        concentration information.")
  }
  
  # Pathway Enrichment 
  pathway_mat <- tools$pathway_enrichment %>%
    filter(combined_id %in% combined_ids, !!sym(pathway_pval_col) < pathway_pval_thresh) %>%
    dplyr::select("combined_id", "Term", "Combined.Score") %>%
    pivot_wider(names_from = "Term", values_from = !!sym(pathway_score_col)) %>%
    left_join(meta %>% dplyr::select(combined_id, !!sym(treatment_ID)) %>% distinct(), by = "combined_id") %>%
    relocate(!!sym(treatment_ID)) %>%
    column_to_rownames(treatment_ID) %>%
    dplyr::select(-combined_id) %>%
    t()
  
  # Gene Expression
  all_de <- bind_rows(tools$diff_exprs, .id = "comparison")
  signif_genes <- all_de %>%
    filter(!!sym(de_padj_col) < de_pval_thresh) %>%
    pull(gene) %>% unique()
  
  gene_mat <- all_de %>%
    filter(gene %in% signif_genes, combined_id %in% combined_ids) %>%
    dplyr::select(combined_id, gene, !!sym(de_metric)) %>%
    pivot_wider(names_from = gene, values_from = !!sym(de_metric)) %>%
    left_join(meta %>% dplyr::select(combined_id, !!sym(treatment_ID)) %>% distinct(), by = "combined_id") %>%
    relocate(!!sym(treatment_ID)) %>%
    column_to_rownames(treatment_ID) %>%
    dplyr::select(-combined_id) %>%
    t()

  
  # Metadata-based descriptors
  has_phenotype <- all(!is.na(metadata_columns)) && length(metadata_columns) > 0
  
  if (has_phenotype) {
      phenotype_mat <- data@meta.data %>%
      filter(combined_id %in% combined_ids) %>%
      dplyr::select(combined_id, all_of(metadata_columns)) %>%
      group_by(combined_id) %>%
      summarise(across(all_of(metadata_columns), ~ median(.x, na.rm = TRUE)), .groups = "drop") %>%
      left_join(
        data@meta.data %>% dplyr::select(combined_id, !!sym(treatment_ID)) %>% distinct(),
        by = "combined_id"
      ) %>%
      relocate({{treatment_ID}}) %>%
      dplyr::select(-combined_id) %>%
      column_to_rownames({{treatment_ID}}) %>%
      t()
  }
  
  # Collect present matrices
  view_matrices <- list(
    genes = gene_mat,
    pathways = pathway_mat
  )
  
  if ("chem_descriptors" %in% names(tools)) {
    # Chemical descriptors
    desc_mat <- tools$chem_descriptors %>%
      dplyr::select(-"clean_compound_name") %>%
      column_to_rownames("Treatment_1") %>%
      t()
    view_matrices$chem_descriptors <- desc_mat
  }
  
  if (has_phenotype && exists("phenotype_mat") && !is.null(phenotype_mat)) {
    phenotype_views <- lapply(rownames(phenotype_mat), function(feature) {
      mat <- phenotype_mat[feature, , drop = FALSE]
      mat
    })
    names(phenotype_views) <- rownames(phenotype_mat)
    view_matrices <- c(phenotype_views, view_matrices)
  }
  
  #add any additional data from @toolls
  known_keys <- c("diff_exprs", "pathway_enrichment", "chem_descriptors")
  additional_views <- tools[setdiff(names(tools), known_keys)]
  
  for (nm in names(additional_views)) {
    tbl <- additional_views[[nm]]
    if (is.data.frame(tbl) && "Treatment" %in% colnames(tbl)) {
      df <- tbl %>%
        column_to_rownames("Treatment") %>%
        t()
      view_matrices[[nm]] <- df
    }
  }
  
  # Match columns
  common_cols <- Reduce(intersect, lapply(view_matrices, colnames))

  # Subset all views to common columns
  views <- lapply(view_matrices, function(mat) {
    mat[, common_cols, drop = FALSE]
  })
  
  # Ensure all views have matching column names
  views <- lapply(views, function(view) {
    colnames(view) <- common_cols
    return(view)
  })
  
  # Remove features (rows) with all NA or zero across all samples
  views <- lapply(views, function(view) {
    is_all_na_or_zero <- function(x) all(is.na(x) | x == 0)
    view[!apply(view, 1, is_all_na_or_zero), , drop = FALSE]
  })
  
  # Create MOFA object
  mofa_obj <- MOFA2::create_mofa(views)
  return(mofa_obj)
}
