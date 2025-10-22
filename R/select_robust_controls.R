#' Select high-quality control replicates via TMMwsp log-CPM correlation
#'
#' @description
#' For a given control group (e.g., DMSO) on a specific plate/batch, this function
#' ranks samples by their average correlation (Fisher z-averaged) to all *other*
#' samples using edgeR's TMMwsp-normalized log2-CPM. It returns the ranking and (optionally)
#' plots per-sample expression distributions and sample-sample correlation heatmaps.
#'
#' @param data A tidyseurat object containing an RNA assay with a **counts** layer.
#' @param samples the control/treatment label to keep in column samples
#'   (e.g., "CB_43_EP73_0"). Only cells/samples with this label are considered.
#' @param orig_ident Character scalar: the plate/batch identifier to keep
#'   (e.g., "VH02012942"). Only cells/samples from this batch are considered.
#' @param cpm_filter Numeric scalar; CPM threshold used for gene filtering prior to
#'   normalization (default 1).
#' @param min_samps Integer; a gene must be expressed (CPM > cpm_filter) in at least
#'   this many samples to be retained (default 16).
#' @param corr_method Correlation type used for ranking; one of
#'   c("spearman","pearson") (default "spearman").
#' @param top_n Integer; the number of top-ranked samples to report in topN.
#'   Ties at the cutoff are kept (default 5).
#' @param make_plots Logical; if TRUE, print a log2-CPM boxplot and Pearson/Spearman
#'   correlation heatmaps (default TRUE).
#'
#' @details
#' Workflow:
#' 1) Subset to the specified samples and orig_ident (plate/batch).
#' 2) Build an edgeR::DGEList, filter lowly expressed genes using CPM and min_samps.
#' 3) Normalize with TMMwsp and compute log2-CPM.
#' 4) Rank samples by mean Fisher z transformed correlation to all other samples
#'    (according to corr_method).
#' 5) Return the ranking, correlation matrices, the normalized matrix, and (optionally)
#'    plots for QC.
#'
#' Column names of the counts matrix are rewritten to "<orig.ident>_<Well_ID>"
#' for easier visual inspection in plots.
#'
#' @return A list with elements:
#' * subset_obj: The Seurat object subset used for analysis.
#' * dge: The filtered edgeR::DGEList
#' * log_cpm_tmm: Matrix of TMMwsp log2-CPM.
#' * boxplot_df: Long-format data frame used for the boxplot (gene, sample, log_cpm).
#' * cor_pearson: Sample-sample Pearson correlation matrix.
#' * cor_spearman: Sample-sample Spearman correlation matrix.
#' * ranking_method: The correlation method used for ranking.
#' * scores_mean_to_others: Named numeric vector of mean Fisher-z back-transformed
#'   correlations (higher = better), sorted decreasing.
#' * topN: Named numeric vector of the top-ranked samples (ties at the cutoff kept).

#' 
#' 
#' @examples
#' data(mini_mac)
#' res <- select_robust_controls(mini_mac,samples = "DMSO_0", orig_ident = "PMMSq033_mini")
#'
#' 
#' @importFrom rlang .data
#' @importFrom edgeR DGEList calcNormFactors cpm
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @export

select_robust_controls <- function(
    data,
    samples,                 # e.g. "CB_43_EP73_0"
    orig_ident,                  # e.g. "VH02012942"
    cpm_filter    = 1,           # CPM threshold for gene filtering
    min_samps     = 16,          # number of samples a gene must be expressed in
    corr_method   = c("spearman","pearson"),
    top_n         = 5,
    make_plots    = TRUE
){
  validate_inputs <- function(data, samples, orig_ident) {
    if (!inherits(data, "Seurat")) {
      stop("argument 'data' must be a Seurat or TidySeurat object.")
    }
    
    # check samples and orig_ident columns
    if (colnames(data@meta.data)%in% c("combined_id","orig.ident") %>% sum() < 2) {
      stop("The 'data' object must contain 'combined_id' and 'orig.ident' columns in its metadata.")
    }
    # check samples in samples column
    if (is.null(samples)){
      stop("Please provide a value for 'samples'.")
    } else if (!all(samples %in% unique(data$combined_id))) {
      stop("Some values in 'samples' are not found in the 'combined_id' column of 'data'.")
    }
    # check orig.ident in the orig.ident column
    if (is.null(orig_ident)){
      stop("Please provide a value for 'orig_ident'.")
    } else if (!orig_ident %in% unique(data$orig.ident)) {
      stop("The value of 'orig_ident' is not found in the 'orig.ident' column of 'data'.")
    }
    return(list(data = data, samples = samples, orig_ident = orig_ident))
  }
  validated <- validate_inputs(data = data, samples = samples, orig_ident = orig_ident)
  data <- validated$data
  group_by <- validated$orig_ident
  samples <- validated$samples
  
  corr_method <- match.arg(corr_method)
  sel_cells <- colnames(data)[data$combined_id == samples &
                                data$orig.ident  == orig_ident]
  if (length(sel_cells) == 0L) {
    stop("No cells/samples match the specified 'combined_id' and 'orig_ident'.")
  }
  subgroup <- subset(data, cells = sel_cells)
  # Counts and human-friendly column names
  counts_d <- Seurat::GetAssayData(subgroup, assay = "RNA", layer = "counts")
  well_colnames <- paste0(subgroup$orig.ident, "_", subgroup$Well_ID)
  names(well_colnames) <- rownames(subgroup@meta.data)
  colnames(counts_d) <- well_colnames[colnames(counts_d)]
  # edgeR container + gene filtering
  y <- edgeR::DGEList(counts_d, group = subgroup$orig.ident)
  keep <- rowSums(edgeR::cpm(y) > cpm_filter) >= min_samps
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMMwsp")
  log_cpm_tmm <- edgeR::cpm(y, log = TRUE, normalized.lib.sizes = TRUE)
  # Long data for boxplot 
  df_long <- as.data.frame(log_cpm_tmm) |>
    tibble::rownames_to_column(var = "gene") |>
    tidyr::pivot_longer(-gene, names_to = "sample", values_to = "log_cpm")
  if (make_plots) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' not available; skipping boxplot.")
    } else {
      print(
        ggplot2::ggplot(df_long, ggplot2::aes(x = .data$sample, y = .data$log_cpm)) +
          ggplot2::geom_boxplot(outlier.size = 0.5) +
          ggplot2::labs(x = "Sample", y = "log2 CPM",
                        title = "Boxplot of log2-CPM (TMMwsp)") +
          ggplot2::theme_classic() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      )
    }
  }
  # Correlation matrices
  cors_pear  <- stats::cor(log_cpm_tmm, method = "pearson")
  cors_spear <- stats::cor(log_cpm_tmm, method = "spearman")
  if (make_plots) {
    if (!requireNamespace("pheatmap", quietly = TRUE)) {
      warning("Package 'pheatmap' not available; skipping heatmaps.")
    } else {
      pheatmap::pheatmap(cors_pear,  main = "Pearson correlation")
      pheatmap::pheatmap(cors_spear, main = "Spearman correlation")
    }
  }
  # Ranking by mean Fisher-z correlation to all *other* samples
  R <- stats::cor(log_cpm_tmm, method = corr_method, use = "pairwise.complete.obs")
  diag(R) <- NA_real_
  # Clip to (-1,1), Fisher z-transform, average, back-transform
  Z <- atanh(pmin(pmax(R, -0.999999), 0.999999))
  score_z <- rowMeans(Z, na.rm = TRUE)
  score_r <- tanh(score_z)
  # Top-N names and scores (keep ties at the cutoff)
  ord   <- order(score_r, decreasing = TRUE, na.last = NA)
  srt   <- score_r[ord]
  if (length(srt) == 0L) {
    stop("No samples available after filtering; adjust 'cpm_filter'/'min_samps'.")
  }
  k      <- min(top_n, length(srt))
  cutoff <- srt[k]
  keepN  <- srt >= cutoff
  topN   <- srt[keepN]
  # Return everything useful
  list(
    subset_obj              = subgroup,
    dge                     = y,
    log_cpm_tmm             = log_cpm_tmm,
    boxplot_df              = df_long,
    cor_pearson             = cors_pear,
    cor_spearman            = cors_spear,
    ranking_method          = corr_method,
    scores_mean_to_others   = sort(score_r, decreasing = TRUE),
    topN                    = topN
  )
}
