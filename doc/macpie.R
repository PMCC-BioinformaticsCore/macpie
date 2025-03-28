## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----load_metadata------------------------------------------------------------

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage

# Load all functions
devtools::load_all()

library(macpie)
library(Seurat)
library(edgeR)
library(dplyr)
library(leidenbase)
library(gridExtra)
library(enrichR)
library(parallel)
library(mcprogress)
library(ggsci)
library(ggiraph)
library(pheatmap)
library(tidyr)
library(tibble)
library(enrichR)
library(variancePartition)
library(glmGamPoi)
library(PoiClaClu)
library(DT)

# Define project variables
project_name <- "PMMSq033"
project_metadata <- system.file("extdata/PMMSq033/PMMSq033_metadata_drugnames.csv", package = "macpie")

# Load metadata
metadata <- read_metadata(project_metadata)
metadata$Time <- as.factor(metadata$Time)
metadata$Concentration_1 <- as.factor(metadata$Concentration_1)
colnames(metadata)

# Validate metadata
validate_metadata(metadata)


## ----metadata_plot, fig.width = 8, fig.height = 6-----------------------------


plot_metadata_heatmap(metadata)


## ----load_data----------------------------------------------------------------

project_rawdata <- system.file("extdata/PMMSq033/raw_matrix", package = "macpie")
raw_counts_total <- Read10X(data.dir = project_rawdata)
keep <- rowSums(cpm(raw_counts_total) >= 10) >= 2
raw_counts <- raw_counts_total[keep, ]

# Create tidySeurat object
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 1,
                          min.features = 1)


## ----violin_plot, fig.width = 8, fig.height = 6-------------------------------
# Calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA")

# Join with metadata
mac <- mac %>%
  inner_join(metadata, by = c(".cell" = "Barcode"))

# Add unique identifier
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

# Example of a function from Seurat QC 
VlnPlot(mac, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 4, group.by = "Sample_type") & 
  scale_fill_manual(values = macpie_colours$discrete) &
  macpie_theme(show_x_title = F, show_y_title = F, legend_position_ = 'none', x_labels_angle = 45)


## ----subset_seurat, fig.width = 10, fig.height = 6----------------------------

unique(mac$Project)
mac <- mac %>%
  filter(Project == "Current")

# QC plot plate layout (all metadata columns can be used):
plot_plate_layout(mac, "nCount_RNA", "Sample_type")


## ----mds_plot, fig.width = 8, fig.height = 8----------------------------------

# Example of MDS function 
p <- plot_mds(mac, group_by = "Sample_type", label = "combined_id", n_labels = 30)
girafe(ggobj = p, fonts = list(sans = "sans"))


## ----plot_rle, fig.width = 8, fig.height = 7----------------------------------

# First we will subset the data to look at control, DMSO samples only
mac_dmso <- mac %>%
  filter(Treatment_1 == "DMSO")

# Run the RLE function
plot_rle(mac_dmso, label_column = "Row")
plot_rle(mac_dmso, label_column = "Row", normalisation = "SCT")
plot_rle(mac_dmso, label_column = "Row", normalisation = "edgeR")


## ----qc_stats, fig.width = 8, fig.height = 6----------------------------------

qc_stats <- compute_qc_metrics(mac, "combined_id")
qc_stats$stats_summary


## ----z_score_lollipop, fig.width = 8, fig.height = 10-------------------------

plot_qc_metrics(qc_stats, "combined_id", "z_score")

## ----distance_plot, fig.width = 8, fig.height = 6-----------------------------

plot_distance(mac, "combined_id", "Staurosporine_10")


## ----de_analysis, fig.width = 8, fig.height = 6-------------------------------

# First perform the differential expression analysis
treatment_samples <- "Staurosporine_0.1"
control_samples <- "DMSO_0"
top_table_edgeR <- compute_single_de(mac, treatment_samples, control_samples, method = "edgeR")

# Let's visualise the results with a volcano plot
plot_volcano(top_table_edgeR)


## ----plot_cpm, fig.width = 8, fig.height = 6----------------------------------
genes <- top_table_edgeR$gene[1:6]
group_by <- "combined_id"
plot_cpm(mac,genes, group_by, treatment_samples, control_samples)

## ----de_single_summary--------------------------------------------------------
datatable(summarise_de(top_table_edgeR, lfc_threshold = 1, padj_threshold = 0.05))

## ----pathway_analysis_single, fig.width = 8, fig.height = 15------------------

top_genes <- top_table_edgeR %>%
  filter(p_value_adj < 0.05) %>%
  select(gene) %>%
  pull()

enriched <- enrichR::enrichr(top_genes, c("MSigDB_Hallmark_2020","DisGeNET",
                                 "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"))
p1 <- enrichR::plotEnrich(enriched[[1]]) + 
  macpie_theme(legend_position_ = 'right') + scale_fill_gradientn(colors = macpie_colours$continuous)
p2 <- enrichR::plotEnrich(enriched[[2]]) + 
   macpie_theme(legend_position_ = 'right') + scale_fill_gradientn(colors = macpie_colours$continuous)
p3 <- enrichR::plotEnrich(enriched[[3]]) + 
   macpie_theme(legend_position_ = 'right') + scale_fill_gradientn(colors = macpie_colours$continuous)

gridExtra::grid.arrange(p1, p2, p3, ncol = 1)


## ----de_multi, fig.width = 8, fig.height = 5----------------------------------
mac$combined_id <- make.names(mac$combined_id)

treatments <- mac %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
mac <- compute_multi_de(mac, treatments, control_samples = "DMSO_0", method = "edgeR", num_cores = 1)


## ----plot_multi_de, fig.width=10, fig.height=6--------------------------------
plot_multi_de(mac, group_by = "combined_id", value = "log2FC", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")


## ----plot_multi_de_lcpm, fig.width=10, fig.height=6---------------------------
plot_multi_de(mac, group_by = "combined_id", value = "lcpm", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")

## ----de_multi_summary---------------------------------------------------------
datatable(summarise_de(mac, lfc_threshold = 1, padj_threshold = 0.01, multi=TRUE))

## ----enriched_pathways, fig.width = 8, fig.height = 12------------------------

# Load genesets from enrichr for a specific species or define your own
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
mac <- compute_multi_enrichr(mac, genesets = enrichr_genesets)

enriched_pathways_mat <- mac@tools$pathway_enrichment %>%
  select(combined_id, Term, Combined.Score) %>%
  pivot_wider(names_from = combined_id, values_from = Combined.Score) %>%
  column_to_rownames(var = "Term") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, log1p(.)))) %>%  # Replace NA with 0 across all columns
  as.matrix()

pheatmap(enriched_pathways_mat, color = macpie_colours$continuous) + macpie_theme()


## ----compute_multi_screen_profile, fig.width = 15, fig.height = 5-------------

mac <- compute_multi_screen_profile(mac, target = "Staurosporine_10", num_cores = 1)
mac_screen_profile <- mac@tools$screen_profile %>%
  mutate(logPadj = c(-log10(padj))) %>%
  arrange(desc(NES)) %>%
  mutate(target = factor(target, levels = unique(target))) 

ggplot(mac_screen_profile, aes(target, NES)) +
  #geom_point(aes(size = logPadj)) +
  geom_point() +
  facet_wrap(~pathway, scales = "free") +
  macpie_theme(x_labels_angle = 45, show_x_title = F)


