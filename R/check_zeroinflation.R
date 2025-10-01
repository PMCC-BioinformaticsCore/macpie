#' @tite Check zero-inflated counts 
#'
#' This offers a fast, method-agonistic way to check whether the data is zero-inflated 
#' relative to a Negative Binomial distribution. 
#' 
#' It computes observed zero fraction and NB-expected zero fraction per gene,
#' then reports a ZI index = p0_obs - p0_NB. Uses edgeR for dispersion.
#' 
#' Other methods to check zero-inflation include: GLM-fitted (glmGamPoi) when
#' your data has strong effects (plate/row/column/treatment) and you want 
#' the most accurate estimation.
#' 


check_zeroinflation <- function(data = NULL,
                                treatment_samples = NULL,
                                control_samples = NULL,
                                batch = 1
                                ){
  
  validate_inputs <- function(data, treatment_samples, control_samples) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    if (is.null(treatment_samples) || is.null(control_samples)) {
      stop("Missing the vectors of treatment and control samples.")
    }
    if (!"combined_id" %in% colnames(data@meta.data)) {
      data <- data %>%
        mutate(combined_id = apply(select(starts_with("Treatment_") | starts_with("Concentration_")),
                                   1, paste, collapse = "_")) %>%
        mutate(combined_id = gsub(" ", "", .data$combined_id))
    }
    if (length(treatment_samples) == 1 && length(control_samples) == 1) {
      treatment_samples_list <- grepl(treatment_samples, data$combined_id)
      control_samples_list <- grepl(control_samples, data$combined_id)
      if (any(sum(treatment_samples_list) == 0, sum(control_samples_list) == 0)) {
        stop("Your treatment and control samples are not in your combined_id column.")
      }
    }
  }
  
  
  count_matrix <- GetAssayData(data, assay = "RNA", layer = "counts")
  count_matrix <- Matrix::Matrix(count_matrix, sparse = TRUE)
  obs_zero <- Matrix::rowMeans(count_matrix == 0)
  
  # Negative binomial expected zeros
  # using edgeR for dispersion estimation
  dge <- edgeR::DGEList(counts = count_matrix)
  dge <- edgeR::calcNormFactors(dge, method = "TMMwsp")
  
  #design matrix
  combined_id <- data$combined_id
  #model matrix with the batch parameter
  model_matrix <- if (length(batch) == 1) model.matrix(~0 + combined_id) else
    model.matrix(~0 + combined_id + batch)
  # tagwise dispersion
  dge <- edgeR::estimateDisp(dge, design = model_matrix)  
  phi <- dge$tagwise.dispersion  # NB variance: mu + phi * mu^2  (phi >= 0)
  
  
  # ---- 4) Build per-sample NB mean mu_gj using TMM-scaled library sizes
  # Effective library sizes
  eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
  # Per-gene proportion q_g = total counts for gene / total eff_lib across samples
  total_eff_lib <- sum(eff_lib)
  total_counts_per_gene <- Matrix::rowSums(count_matrix)
  q_g <- as.numeric(total_counts_per_gene) / total_eff_lib
  # expected mean cont for each gene in each sample
  mu <- q_g %o% eff_lib      # (genes x samples) without materializing too big objects in typical 384-well
  
  
  # ---- 5) NB-expected zero per gene, averaged across samples
  # p0_NB_gj = (1 + phi_g * mu_gj)^(-1/phi_g); then average over j
  # guard against phi==0 (Poisson limit): use exp(-mu) when phi ~ 0
  eps <- 1e-12
  phi_safe <- pmax(phi, eps)
  # vectorized: compute per gene using rowMeans on transformed mu
  inv_phi <- 1 / phi_safe
  # (1 + phi * mu)^(-1/phi)  -> do per gene with broadcasting
  # Use sweep for efficiency
  one_plus <- sweep(mu, 1, phi_safe, `*`)
  one_plus <- 1 + one_plus
  p0_nb_mat <- sweep(one_plus, 1, inv_phi, `^`)  # (1 + phi*mu)^(-1/phi)
  p0_nb <- rowMeans(1 / p0_nb_mat)               # since above computed (1+phi*mu)^(+1/phi); invert
  # NOTE: If memory is a concern, chunk over columns; okay for 384 wells.
  
  # Poisson fallback where phi was ~0 (replace those rows with exp(-rowMeans(mu)))
  poi_idx <- which(phi < 1e-8)
  if (length(poi_idx)) {
    mu_bar <- rowMeans(mu[poi_idx, , drop = FALSE])
    p0_nb[poi_idx] <- exp(-mu_bar)
  }
  
  # ---- 6) ZI index per gene
  zi <- obs_zero - p0_nb
  
  # ---- 7) Summaries
  libsize <- Matrix::colSums(count_matrix)
  summary <- list(
    n_genes          = nrow(count_matrix),
    n_samples        = ncol(count_matrix),
    median_libsize   = stats::median(libsize),
    median_p0_obs    = stats::median(obs_zero),
    median_p0_nb     = stats::median(p0_nb),
    median_ZI        = stats::median(zi),
    pct_ZI_gt_0.05   = mean(zi > 0.05),
    pct_ZI_gt_0.10   = mean(zi > 0.10),
    mean_observed_zeros       = mean(obs_zero*ncol(count_matrix)),
    mean_expected_zeros       = mean(p0_nb*ncol(count_matrix))
  )
  
  thresholds_used <- list(
    zi_margin_small = 0.05,
    zi_margin_mod   = 0.10
  )
  
  gene_metrics <- data.frame(
    gene = rownames(count_matrix),
    mean_count = total_counts_per_gene / ncol(count_matrix),
    dispersion = phi,
    p0_obs = p0_obs,
    p0_nb  = p0_nb,
    ZI     = zi,
    stringsAsFactors = FALSE
  )
  
  list(
    gene_metrics   = gene_metrics,
    summary        = summary,
    thresholds_used = thresholds_used
  )
  
}