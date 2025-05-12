# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation
# Set seed
set.seed(123456)

# Load libraries ===================================================
library(devtools)
library(tidyverse)
library(Seurat)
library(gtools)
library(tidyseurat)
library(DESeq2)
library(edgeR)
library(RUVSeq)
library(variancePartition)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(lintr)
library(enrichR)
library(jsonlite)
library(parallel)
library(mcprogress)
library(httr2)
library(clusterProfiler)
library(ggsci)
library(pheatmap)
library(umap)
library(doParallel)
library(pbapply)
library(zinbwave)
#library(SingleCellExperiment)
library(gdtools)
library(ggiraph)
library(drc)
library(webchem)
library(randomForest)
library(MOFA2)

# Define longer length for description files
custom_linters <- lintr::linters_with_defaults(
  line_length_linter = lintr::line_length_linter(120) # Set max line length to 120
)

# Load all previous libraries
load_all()

# Select which sample dataset to use
project_name <- "PMMSq033"

# Directory with data
data_dir <- "inst/extdata"

# Metadata ===================================================
# Load metadata
project_metadata <- system.file("extdata/PMMSq033/PMMSq033_metadata_drugnames.csv", package = "macpie")
metadata <- read_metadata(project_metadata)

#load in the metadata on cell viability and confluence
cell_viability <- read.csv("inst/extdata/PMMSq033/PMMSq033_CTG_cellcount.csv") %>%
  select("Well_ID","CTG","Confluence")

metadata <- metadata %>%
  left_join(., cell_viability, join_by(Well_ID))

# Validate metadata
validate_metadata(metadata)

# Heatmap visualisation of metadata
plot_metadata_heatmap(metadata)

# Only create an id if there are multiple plates
#TODO: test multiple plates
if(length(project_metadata) > 1){
  metadata <- metadata %>%
    mutate(id = gsub("Plate", "", Plate_ID)) %>%
    mutate(Barcode = paste0(id, "_", Barcode))
}

################## reads to Seurat object ##################
# TO-DO: check for multiple folders of data that should be imported at the same time
project_rawdata <- file.path(data_dir, project_name, "raw_matrix")
raw_counts_total <- Read10X(data.dir = project_rawdata)
keep <- rowSums(cpm(raw_counts_total) >= 10) >= 2
raw_counts <- raw_counts_total[keep, ]

# Create tidySeurat object
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 5,
                          min.features = 5)

# Calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|
                                                ^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#TO-DO: verify that the mouse and human regexp for ribo proteins is correct (grepl("^RP",row.names(mac@assays$RNA$counts)))

# Add metadata to the sequencing data and select only a specific project
mac <- mac %>%
  inner_join(metadata,by = c(".cell" = "Barcode")) %>%
  mutate("Barcode" = .cell) %>%
  filter(Project == "Current")

# Create an ID that uniquely identifies samples based on the
# Combination of treatment and treatment concentration
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id)) %>%
  mutate(combined_id = make.names(combined_id))


# QC ===================================================
# QC plot plate layout (all metadata columns can be used):
plot_plate_layout(mac,"nCount_RNA","Sample_type")

# Example of Seurat function being used
VlnPlot(mac,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"),
        group.by = "Sample_type",
        ncol = 4)

# Example of MDS function, using limma
p <- plot_mds(mac)
girafe(ggobj = p) %>%
  girafe_options(
  opts_tooltip(css = "font-family: 'Arial', sans-serif; font-size: 14px; background: white; border-radius: 5px; padding: 5px;")
)

# Compare normalisation methods using the RLE function
#To-do: verify that different plates would be plotted side-by-side
#To-do: add all Xin's plots
#To-do: QC summary

# RLE plot
mac_dmso <- mac %>%
  filter(Treatment_1 == "DMSO")
plot_rle(mac_dmso, label_column = "Row", normalisation = "limma_voom")

# Differential expression ===================================================

## DE Single -----------------

treatment_samples <- "Staurosporine_10"
control_samples <- "DMSO_0"

# Perform differential expression
top_table <- compute_single_de(mac, treatment_samples, control_samples, method = "limma_voom")
plot_volcano(top_table)

# Perform pathway enrichment
top_genes <- top_table %>%
  filter(p_value_adj<0.01) %>%
  select(gene) %>%
  pull()

# Basic pathway enrichment (MSigDB_Hallmark_2020, LINCS_L1000_CRISPR_KO_Consensus_Sigs)
enrichR_pathway <- "MSigDB_Hallmark_2020"
enrichr_results <- enrichr(top_genes, enrichR_pathway)
plotEnrich(enrichr_results[[1]], title = enrichR_pathway)

# Process locally
enriched <- compute_single_enrichr(top_genes, "MSigDB_Hallmark_2020", species = "human")
plotEnrich(enriched)

## Screen-level analysis --------------------------

### DGE analysis ------------------
#to-do: shorten
treatments <- mac %>%
  filter(Concentration_1 == "10") %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()

# Load genesets from enrichr for a specific species or define your own
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")

# Update the mac object with differential expression
mac <- compute_multi_de(mac, treatments, control_samples = "DMSO_0", method = "limma_voom")
mac <- compute_multi_enrichr(mac, genesets = enrichr_genesets)
mac <- compute_multi_screen_profile(mac, target = "Staurosporine_10")


# Get all the differential expression information in a tabular format
de_genes_per_comparison <- bind_rows(mac@tools$diff_exprs)
enriched_pathways_per_comparison <- mac@tools$pathway_enrichment
screen_profile_per_comparison <- mac@tools$screen_profile

### Plot DGE results ------------------
# Plot the results for pathways across all the comparisons
enriched_pathways_mat <- enriched_pathways_per_comparison %>%
  select(combined_id, Term, Combined.Score) %>%
  pivot_wider(
    names_from = combined_id,
    values_from = Combined.Score
  ) %>%
  column_to_rownames(var = "Term") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, log1p(.)))) %>%  # Replace NA with 0 across all columns
  as.matrix()
pheatmap(enriched_pathways_mat)

# Screen for similarity of profiles -----------
# to-do add coloring by data type
screen_profile_per_comparison %>%
  mutate(logPadj = c(-log10(padj))) %>%
  arrange(desc(NES)) %>%
  mutate(target = factor(target, levels = unique(target))) %>%
  ggplot(.,aes(target,NES))+
  #geom_point(aes(size = logPadj)) +
  geom_point() +
  facet_wrap(~pathway,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


######## Cheminformatics
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
#to-do: shorten
treatments <- mac %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
# Update the mac object with differential expression
mac <- compute_multi_de(mac, treatments, control_samples = "DMSO_0", method = "limma_voom", num_cores = 1)
mac <- compute_multi_enrichr(mac, genesets = enrichr_genesets)

# Add smiles (warning this can take a while)
mac <- compute_smiles(mac)

# Calculate descriptors
mac <- compute_chem_descriptors(mac)

# Join with target variable (e.g. pathway score)
model_df <- mac@tools$pathway_enrichment %>%
  filter(Term == "Estrogen Response Early") %>%
  left_join(., mac@meta.data, join_by(combined_id)) %>%
  filter(Concentration_1 == 10) %>%
  select(Treatment_1, Combined.Score) %>%
  unique() %>%
  left_join(., mac@tools$chem_descriptors, join_by(Treatment_1)) %>%
  drop_na()

# Train random forest
rf_model <- randomForest(Combined.Score ~ ., data = model_df, importance = TRUE, na.action = na.omit)

# Get importance scores
rf_importance <- importance(rf_model, type = 1)  # %IncMSE = predictive power
rf_ranked <- sort(rf_importance[, 1], decreasing = TRUE)

# Top 20 important descriptors
head(rf_ranked, 20)
varImpPlot(rf_model, n.var = 20, main = "Top 20 Random Forest Features")


# Rank features by lowest p-value
lm_ranked <- sort(lm_pvals, na.last = NA)
head(lm_ranked, 20)  # Top 20 most significant predictors


# Convert fingerprints to binary matrix
#fp_matrix <- do.call(rbind, lapply(fp_list, function(x) {
#  fp_vec <- rep(0, 1024)
#  if (!is.null(x)) fp_vec[as.numeric(x)] <- 1
#  fp_vec
#}))
#
#colnames(fp_matrix) <- paste0("FP_", seq_len(ncol(fp_matrix)))
#fingerprint_df <- as.data.frame(fp_matrix)
#fingerprint_df$clean_compound_name <- names(fp_list)
#
## Combine all features
#combined_features <- descriptor_df %>%
#  left_join(.,fingerprint_df, join_by(clean_compound_name))


############ MOFA

#load external data
file_path <- system.file("extdata", "PMMSq033/PMMSq033_CTG_cellcount.csv", package = "macpie")

cell_viability <- read.csv(file_path) %>%
  mutate(
    cell_viability = scale(cell_viability)[, 1],
    cell_confluence = scale(cell_confluence)[, 1]
  ) %>%
  dplyr::select(Well_ID, cell_viability, cell_confluence)

#add to metadata
mac@meta.data <- mac@meta.data %>%
  rownames_to_column("Barcode") %>%
  left_join(cell_viability, by = "Well_ID") %>%
  column_to_rownames("Barcode")



compounds_10um <- mac %>%
  filter(Concentration_1 == 10) %>%
  dplyr::select(combined_id) %>%
  pull()
         
pathway_mat <- mac@tools$pathway_enrichment %>%
  filter(combined_id %in% compounds_10um) %>%
  filter(Adjusted.P.value < 0.05) %>%
  dplyr::select(combined_id, Term, Combined.Score) %>%
  pivot_wider(names_from = Term, values_from = Combined.Score) %>%
  column_to_rownames("combined_id") %>%
  as.data.frame() %>%
  rownames_to_column("combined_id") %>%
  left_join(
    mac@meta.data %>% dplyr::select(combined_id, Treatment_1) %>% distinct(),
    by = "combined_id"
  ) %>%
  dplyr::select(-combined_id) %>%
  relocate(Treatment_1) %>%
  column_to_rownames("Treatment_1") %>%
  t()

pathway_mat[is.na(pathway_mat)]<-0

# 1. Combine all DE results
all_de <- bind_rows(mac@tools$diff_exprs, .id = "comparison")

# 2. Identify genes that are significant in at least one comparison
signif_genes <- all_de %>%
  filter(p_value_adj < 0.01) %>%
  #filter(grepl("ENSG", gene)) %>%
  pull(gene) %>%
  unique() 

# 3. Filter and reshape the matrix for only significant genes and relevant compounds
gene_mat <- all_de %>%
  filter(gene %in% signif_genes, combined_id %in% compounds_10um) %>%
  dplyr::select(combined_id, gene, metric) %>%  # Change 'logFC' to your desired metric
  pivot_wider(names_from = gene, values_from = metric) %>%
  column_to_rownames("combined_id") %>%
  as.data.frame() %>%
  rownames_to_column("combined_id") %>%
  left_join(
    mac@meta.data %>% dplyr::select(combined_id, Treatment_1) %>% distinct(),
    by = "combined_id"
  ) %>%
  dplyr::select(-combined_id) %>%
  relocate(Treatment_1) %>%
  column_to_rownames("Treatment_1") %>%
  t()

# 2. View 2: chemical descriptors
desc_mat <- mac@tools$chem_descriptors %>%
  dplyr::select(-clean_compound_name) %>%
  column_to_rownames("Treatment_1") %>%
  t() %>%
  scale()

# 3. Optional View 3: chemical fingerprints
#fp_mat <- mac@tools$chem_fingerprints %>%
#  column_to_rownames("Treatment_1") %>%
#  t()

reads_counts <- mac %>%
  filter(Concentration_1 == 10 | Treatment_1 == "DMSO") %>%
  dplyr::select(Treatment_1, nCount_RNA) %>%
  group_by(Treatment_1) %>%
  summarise(read_count = log10(median(nCount_RNA))) %>%
  ungroup() %>%
  deframe()

cell_viability <- mac %>%
  filter(Concentration_1 == 10 | Treatment_1 == "DMSO") %>%
  dplyr::select(Treatment_1, cell_viability) %>%
  group_by(Treatment_1) %>%
  summarise(cell_viability = median(cell_viability)) %>%
  ungroup() %>%
  deframe() 

cell_count <- mac %>%
  filter(Concentration_1 == 10 | Treatment_1 == "DMSO") %>%
  dplyr::select(Treatment_1, cell_confluence) %>%
  group_by(Treatment_1) %>%
  summarise(cell_viability = median(cell_confluence)) %>%
  ungroup() %>%
  deframe() 

reads_fc <- reads_counts - reads_counts["DMSO"]
viab_fc <- cell_viability 
cell_count <- cell_count 


# 4. Match sample names across views
common_cols <- Reduce(intersect, 
                      list(names(reads_counts), 
                           names(reads_fc),
                           names(viab_fc),
                           names(cell_count),
                           colnames(gene_mat), colnames(pathway_mat), colnames(desc_mat)))  # add fp_mat if used

views <- list(
  #reads = t(as.matrix(reads_fc[common_cols])),
  #viability = t(as.matrix(viab_fc[common_cols])),
  cell_count = t(as.matrix(cell_count[common_cols])),
  genes = gene_mat[, common_cols],
  pathways = pathway_mat[, common_cols],
  descriptors = desc_mat[, common_cols]
  # , fingerprints = fp_mat[, common_cols]
)

# Set proper colnames on all matrices
views <- lapply(views, function(view) {
  colnames(view) <- common_cols
  return(view)
})

mofa_obj <- create_mofa(views)
model_opts <- get_default_model_options(mofa_obj)
train_opts <- get_default_training_options(mofa_obj)
train_opts$seed <- 1  # for reproducibility

mofa_obj <- prepare_mofa(mofa_obj, model_options = model_opts, training_options = train_opts)
model <- run_mofa(mofa_obj, use_basilisk = TRUE)

plot_factors(model, color_by = "group")  # coloring by treatment groups
plot_weights(model, view = "pathways", factor = 1)
plot_factor_cor(model)

p<-plot_variance_explained(model)
p<-p+theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("variance_explained.2.pdf", plot = p, width = 10, height = 8, units = "cm")

p<-plot_top_weights(model,
                 view = "pathways",
                 factor = 2,
                 nfeatures = 10
)
ggsave("pathwasy_factor2.pdf", plot = p, width = 24, height = 12, units = "cm")

p<-plot_top_weights(model,
                    view = "pathways",
                    factor = 1,
                    nfeatures = 10
)
ggsave("genes_factor2.pdf", plot = p, width = 24, height = 12, units = "cm")

p<-plot_top_weights(model,
                    view = "genes",
                    factor = 1,
                    nfeatures = 10
)
ggsave("genes_factor2.pdf", plot = p, width = 7, height = 8, units = "cm")

p<-plot_top_weights(model,
                    view = "descriptors",
                    factor = 4,
                    nfeatures = 10
)
ggsave("descriptor_factor5.pdf", plot = p, width = 7, height = 8, units = "cm")



p <- plot_factors(model, 
                  factors = c(1,5), 
                  color_by = "PSMD2",
                  dot_size = 2.5,
                  show_missing = T
)


plot_factors(model, 
            factors = c(2,5)
)

plot_variance_explained(model, plot_total = T)


plot_factors(model, factors = c(1, 5), color_by = "sample")


plot_factors_macpie<-function (object, factors = c(1, 2), groups = "all", show_missing = TRUE, 
                               scale = FALSE, color_by = NULL, shape_by = NULL, size_by = NULL,
                               color_name = NULL, shape_name = NULL, dot_size = 2, alpha = 1, 
                               legend = TRUE, stroke = NULL, return_data = FALSE) 
{
  if (!is(object, "MOFA")) 
    stop("'object' has to be an instance of MOFA")
  
  if (length(unique(factors)) == 1) {
    .args <- as.list(match.call()[-1])
    .args <- .args[names(.args) != "factors"]
    return(do.call(plot_factor, c(.args, list(factors = unique(factors)))))
  } else if (length(factors) > 2) {
    .args <- as.list(match.call()[-1])
    p <- do.call(.plot_multiple_factors, .args)
    return(p)
  }
  
  if (!is.null(color_by) && (length(color_by) == 1) && is.null(color_name)) 
    color_name <- color_by
  if (!is.null(shape_by) && (length(shape_by) == 1) && is.null(shape_name)) 
    shape_name <- shape_by
  
  factors <- .check_and_get_factors(object, factors)
  Z <- MOFA2::get_factors(object, factors = factors, groups = groups, as.data.frame = TRUE)
  
  #color_by <- .set_colorby(object, color_by)
  #shape_by <- .set_shapeby(object, shape_by)
  #size_by <- .set_colorby(object, size_by)
  
  Z <- Z[complete.cases(Z), ]
  df <- merge(Z, color_by, by = "sample")
  df <- merge(df, shape_by, by = "sample")
  df$shape_by <- as.character(df$shape_by)
  
  # Optional size_by merge
  if (!is.null(size_by)) {
    size_df <- data.frame(sample = size_by$sample, size_val = as.numeric(size_by$color_by))
    df <- merge(df, size_df, by = "sample", all.x = TRUE)
  } else {
    df$size_val <- dot_size  # fallback to default
  }
  
  if (isFALSE(show_missing)) 
    df <- filter(df, !is.na(color_by) & !is.na(shape_by))
  
  df <- spread(df, key = "factor", value = "value")
  df <- df[, c(colnames(df)[seq_len(4)], "size_val", factors)]
  df <- set_colnames(df, c(colnames(df)[seq_len(4)], "size", "x", "y"))
  
  if (scale) {
    df$x <- df$x / max(abs(df$x))
    df$y <- df$y / max(abs(df$y))
  }
  
  if (return_data) 
    return(df)
  
  if (is.null(stroke)) {
    stroke <- .select_stroke(N = length(unique(df$sample)))
  }
  
  p <- ggplot(df, aes(x = .data$x, y = .data$y, fill = .data$color_by, 
                      shape = .data$shape_by, size = .data$size)) +
    geom_point(alpha = alpha, stroke = stroke) +
    labs(x = factors[1], y = factors[2]) +
    theme_classic() +
    theme(
      axis.text = element_text(size = rel(0.8), color = "black"),
      axis.title = element_text(size = rel(1.1), color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5)
    )
  
  p <- .add_legend(p, df, legend, color_name, shape_name)
  if (!is.null(color_name)) p <- p + labs(fill = color_name)
  if (!is.null(shape_name)) p <- p + labs(shape = shape_name)
  p <- p + guides(size = guide_legend(title = "Dot Size"))
  
  return(p)
}

temp<-p



ggsave("MDEC.33_PSMD2.pdf", plot = p, width = 24, height = 12, units = "cm")



############ PROCEDURE TO MAKE A FUNCTION
#1. open terminal and pull from github
#git pull origin main

#2. create a branch in the format: function_author
#git checkout -b function_author

#3. in the console: create a function
#use_r(functionX)
#4. edit script
#load_all()
#5. check if the function exists in the global environment
#exists("plote_plot_plate_layout", where = globalenv(), inherits = FALSE)
#6. check and document
#check()
#7. click inside the function and
#Code > Insert roxygen skeleton.
#document()
#check()
#lint(filename="R/multi_enrich_pathways.R",linters = custom_linters)
#8. make the test
#use_test("functionX")
#test()
#9. update vignette and create html
#devtools::build_vignettes()
#10. terminal: git
#git add .
#git commit -m "message"
#git push origin function_author
#git checkout main
#git branch -d function_author


