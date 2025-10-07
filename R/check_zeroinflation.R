#' Quick group-aware zero-inflation check (NB baseline via edgeR)
#'
#' @description
#' Computes a **group-aware Zero-Inflation (ZI) index** for each gene using a
#' negative-binomial (NB) baseline fitted with **edgeR**. For each group
#' (e.g., drug condition), the function:
#' 1) estimates gene-wise tagwise dispersions with edgeR (using all selected groups),
#' 2) builds NB-expected zero probabilities from TMMwsp-scaled means, and
#' 3) returns per-gene ZI (observed zeros minus NB-expected zeros) and
#'    per-group summaries (e.g., % genes with ZI > 0.05).
#'
#' This is intended as a **fast screening diagnostic** to decide whether
#' standard NB GLM methods (edgeR/DESeq2) are adequate or whether a
#' zero-aware workflow (e.g., ZINB-WaVE) might be warranted.
#'
#' @important
#' This helper **relies on edgeR** to estimate dispersion. The current
#' implementation requires **≥2 groups** in the design so that edgeR can
#' stabilize gene-wise dispersions across groups. If you only have a single
#' group and still want a design-aware baseline for expected zeros, fit a
#' Gamma–Poisson/NB GLM yourself (e.g., with **glmGamPoi**) and compute the
#' expected zero probabilities from its fitted means and overdispersion.
#'
#' @param data Seurat object.
#' @param group_by Character, column in `data@meta.data` that defines groups
#'   (default: `"combined_id"`).
#' @param samples Character vector of group labels/patterns to include. If
#'   `NULL` or if none match, all groups in `group_by` are used.
#' @param batch Optional batch indicator; if length 1, an intercept-free design
#'   is used with group dummies.
#'
#' @returns A list with:
#' * `gene_metrics_by_group`: long data frame (group × gene) with `p0_obs`,
#'   `p0_nb`, `ZI`, and counts.
#' * `summary_by_group`: one row per group with medians and % ZI thresholds,
#'   plus observed/expected zero **counts** for the group.
#'
#' @note
#' - This is a **screening** tool; it is not a replacement for fitting a full
#'   GLM with your actual design. If strong covariates exist, a GLM baseline
#'   (e.g., `glmGamPoi::glm_gp`) will yield more faithful expected-zero rates.
#' - For single-group experiments, consider either adding a reference group or
#'   switching to a GLM-based baseline that does not require multiple groups.
#'
#' @examples
#' # check_zeroinflation(mini_mac, group_by = "combined_id",
#' #                           samples = c("DrugA_10", "DMSO_0"))



check_zeroinflation <- function(data = NULL,
                                     group_by = NULL,
                                     samples = NULL,
                                     batch = 1
){
  validate_inputs <- function(data, group_by, samples) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    group_by <- if (is.null(group_by)) "combined_id" else group_by
    
    # check samples in combined_id column
    meta_groups <- as.character(data@meta.data[[group_by]])
    matched_groups <- !is.null(samples) && any(grepl(samples, meta_groups))
    if (!matched_groups){
      # all samples included
      samples <- unique(data@meta.data[[group_by]])
      cat("All samples will be included in the combined_id column.")
    } 
    # need at least two groups for edgeR dispersion estimation
    if (length(samples) == 1) {
      stop("Two treatment groups are needed to calculate dispersion using edgeR.")
    } 
    return(list(data = data, group_by = group_by, samples = samples))
  }
  
  
  
  validated <- validate_inputs(data, group_by, samples)
  data <- validated$data
  group_by <- validated$group_by
  samples <- validated$samples
  
  # samples must be specified
  mac_data <- subset(data, subset = combined_id %in% samples)
  count_matrix <- GetAssayData(mac_data, assay = "RNA", layer = "counts")
  count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
  obs_zero <- Matrix::rowMeans(count_matrix == 0)
  
  # Negative binomial expected zeros
  # using edgeR for dispersion estimation
  dge <- edgeR::DGEList(counts = count_matrix)
  dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")
  
  #design matrix
  combined_id <- mac_data$combined_id
  #make up batch parameter
  model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
    model.matrix(~0 + combined_id + batch)
  # tagwise dispersion
  dge <- edgeR::estimateDisp(dge, design = model_matrix)  
  phi <- dge$tagwise.dispersion  # NB variance: mu + phi * mu^2  (phi >= 0)
  # Build per-sample NB mean mu_gj using TMMwsp-scaled library sizes
  # Effective library sizes
  eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
  
  per_group_gene_metrics <- lapply(samples, function(g){
    
    idx <- which(combined_id == g)
    n_wells <- length(idx)
    # sub count matrix for group g
    count_matrix_g <- count_matrix[, idx, drop=FALSE]
    # Observed zeros within group g
    p0_obs_g <- Matrix::rowMeans(count_matrix_g==0)
    
    # count zeros per gene within group g
    # sum later for summary
    obs_zero_num_g <- Matrix::rowSums(count_matrix_g==0)
    
    # Group-specific q_{g,g} using only wells in group g
    eff_lib_g <- eff_lib[idx]
    total_eff_lib_g <- sum(eff_lib_g)
    total_counts_per_gene_g <- Matrix::rowSums(count_matrix_g)
    q_g_g <- as.numeric(total_counts_per_gene_g) / total_eff_lib_g
    
    # NB-expected zeros within group g (average over wells in g)
    eps <- 1e-12
    phi_safe <- pmax(phi, eps)
    inv_phi <- 1 / phi_safe
    
    # Fast loop over wells in g, no GxJ materialization
    p0_nb_sum_g <- numeric(nrow(count_matrix))
    for (j in seq_along(idx)) {
      Lj <- eff_lib_g[j]
      mu_gj <- q_g_g * Lj
      p0_nb_sum_g <- p0_nb_sum_g + (1 + phi_safe * mu_gj)^(-inv_phi)
    }
    p0_nb_g <- p0_nb_sum_g / length(idx)
    
    # Poisson fallback where phi ~ 0
    poi_idx <- which(phi < 1e-8)
    if (length(poi_idx)) {
      mu_bar_g <- q_g_g * mean(eff_lib_g)
      p0_nb_g[poi_idx] <- exp(-mu_bar_g[poi_idx])
    }
    
    # ZI within group g
    zi_g <- p0_obs_g - p0_nb_g
    
    data.frame(
      group = g,
      gene  = rownames(count_matrix),
      mean_count_group = total_counts_per_gene_g / length(idx),
      dispersion = phi,
      p0_obs = p0_obs_g,
      obs_zeros_num = obs_zero_num_g,
      p0_nb  = p0_nb_g,
      expected_zeros_num = p0_nb_g*n_wells,
      ZI     = zi_g,
      stringsAsFactors = FALSE
    )
  })
  
  gene_metrics_by_group <- do.call(rbind, per_group_gene_metrics)
  
  # Per-group summaries (one row per group)
  summary_by_group <- do.call(rbind, lapply(split(gene_metrics_by_group, gene_metrics_by_group$group), function(df){
    list(
      group            = unique(df$group),
      n_genes          = nrow(df),
      n_wells          = sum(combined_id == unique(df$group)),
      median_p0_obs    = median(df$p0_obs),
      median_p0_nb     = median(df$p0_nb),
      median_ZI        = median(df$ZI),
      pct_ZI_gt_0.05   = mean(df$ZI > 0.05),
      pct_ZI_gt_0.10   = mean(df$ZI > 0.1),
      observed_zeros_num       = sum(df$obs_zeros_num),
      expected_zeros_num  = sum(df$expected_zeros_num) 
    )
  })) |>
    as.data.frame()
  
  # Return gene-level metrics and sample group-level summaries
  list(
    gene_metrics_by_group = gene_metrics_by_group %>% head(10),     
    summary_by_group      = summary_by_group
  )
  
}