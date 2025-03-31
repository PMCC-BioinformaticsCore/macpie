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
library(Matrix)


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


## ----violin_plot, fig.width = 8, fig.height = 4-------------------------------
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
  scale_fill_manual(values = macpie_colours$discrete) 


## ----subset_seurat, fig.width = 8, fig.height = 3-----------------------------
unique(mac$Project)
mac <- mac %>%
  filter(Project == "Current")

# QC plot plate layout (all metadata columns can be used):
p <- plot_plate_layout(mac, "nCount_RNA", "combined_id")
girafe(ggobj = p, 
  fonts = list(sans = "sans"),
  options = list(
    opts_hover(css = "stroke:orange; stroke-width:1px;")  # <- slight darkening
  ))

## ----mds_plot, fig.width = 8, fig.height = 6----------------------------------
p <- plot_mds(mac, group_by = "Sample_type", label = "combined_id", n_labels = 30)
girafe(ggobj = p, fonts = list(sans = "sans"))


## ----qc_stats_plot, fig.width = 8, fig.height = 4-----------------------------
qc_stats <- compute_qc_metrics(mac, group_by = "combined_id", order_by = "median")
qc_stats$stats_summary

## ----qc_metrics, fig.width = 8, fig.height = 4--------------------------------

plot_qc_metrics_heatmap(qc_stats$stats_summary)


## ----fig.width = 8, fig.height = 6--------------------------------------------
plot_distance(mac, "combined_id", treatment = "DMSO_0")

## ----plot_rle, fig.width = 8, fig.height = 4----------------------------------
# First we will subset the data to look at control, DMSO samples only
mac_dmso <- mac %>%
  filter(Treatment_1 == "DMSO")

# Run the RLE function
plot_rle(mac_dmso, label_column = "Row", normalisation = "raw")
plot_rle(mac_dmso, label_column = "Row", normalisation = "edgeR")

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
summarise_de(top_table_edgeR, lfc_threshold = 1, padj_threshold = 0.05)

## ----pathway_analysis_single, fig.width = 8, fig.height = 15------------------

top_genes <- top_table_edgeR %>%
  filter(p_value_adj < 0.05) %>%
  select(gene) %>%
  pull()

enriched <- enrichR::enrichr(top_genes, c("MSigDB_Hallmark_2020","DisGeNET",
                                 "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"))
p1 <- enrichR::plotEnrich(enriched[[1]]) + 
  macpie_theme(legend_position_ = 'right') + 
  scale_fill_gradientn(colors = macpie_colours$continuous)
p2 <- enrichR::plotEnrich(enriched[[2]]) + 
  macpie_theme(legend_position_ = 'right') + 
  scale_fill_gradientn(colors = macpie_colours$continuous)
p3 <- enrichR::plotEnrich(enriched[[3]]) + 
  macpie_theme(legend_position_ = 'right') + 
  scale_fill_gradientn(colors = macpie_colours$continuous)

gridExtra::grid.arrange(p1, p2, p3, ncol = 1)


## ----de_multi, fig.width = 8, fig.height = 5----------------------------------
mac$combined_id <- make.names(mac$combined_id)

treatments <- mac %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
mac <- compute_multi_de(mac, treatments, control_samples = "DMSO_0", method = "edgeR", num_cores = 2)


## ----plot_multi_de, fig.width=10, fig.height=6--------------------------------
plot_multi_de(mac, group_by = "combined_id", value = "log2FC", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")


## ----plot_multi_de_lcpm, fig.width=10, fig.height=6---------------------------
plot_multi_de(mac, group_by = "combined_id", value = "lcpm", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")

## ----de_multi_summary---------------------------------------------------------
summarise_de(mac, lfc_threshold = 1, padj_threshold = 0.01, multi=TRUE)

## ----enriched_pathways, fig.width = 8, fig.height = 12------------------------

# Load genesets from enrichr for a specific species or define your own
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
mac <- compute_multi_enrichr(mac, genesets = enrichr_genesets)

enriched_pathways_mat <- mac@tools$pathway_enrichment %>%
  bind_rows() %>%
  select(combined_id, Term, Combined.Score) %>%
  pivot_wider(names_from = combined_id, values_from = Combined.Score) %>%
  column_to_rownames(var = "Term") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, log1p(.)))) %>%  # Replace NA with 0 across all columns
  as.matrix()

pheatmap(enriched_pathways_mat, color = macpie_colours$continuous_rev) + macpie_theme()


## ----compute_multi_screen_profile, fig.width = 8, fig.height = 5--------------

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


## ----load_two_plates----------------------------------------------------------
project_metadata <- system.file("extdata/PMMSq033/PMMSq033_metadata_drugnames.csv", package = "macpie")
metadata <- read_metadata(project_metadata)

pmm33_in <- system.file("extdata/PMMSq033/raw_matrix", package = "macpie")
pmm34_in <- system.file("extdata/PMMSq034/raw_matrix", package = "macpie")
raw_counts_total <- Read10X(data.dir = c(pmm33_in,pmm34_in))
keep <- rowSums(cpm(raw_counts_total) >= 10) >= 2
raw_counts <- raw_counts_total[keep, ]
combined <- CreateSeuratObject(counts=raw_counts,
                               min.cells = 1,
                          min.features = 1)

combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-|^MT-")
combined[["percent.ribo"]] <- PercentageFeatureSet(combined, pattern = "^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA")

#join with metadata
combined$Barcode <- str_replace_all(rownames(combined@meta.data),"[1|2]_","")

combined <- combined %>%
  inner_join(metadata, by = c("Barcode" = "Barcode"))

combined <- combined%>%
  filter(Project == "Current")

#add unique identifier
combined <- combined %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

combined$combined_id <- make.names(combined$combined_id)


## ----two_plates_plate_layout_plate1, fig.width = 8, fig.height = 6------------
p <- plot_plate_layout(combined%>%filter(orig.ident==1), "nCount_RNA", "combined_id")
girafe(ggobj = p, 
  fonts = list(sans = "sans"),
  options = list(
    opts_hover(css = "stroke:orange; stroke-width:1px;")  # <- slight darkening
  ))

## ----two_plates_plate_layout_plate2, fig.width = 8, fig.height = 6------------
p <- plot_plate_layout(combined%>%filter(orig.ident==2), "nCount_RNA", "combined_id")
girafe(ggobj = p, 
  fonts = list(sans = "sans"),
  options = list(
    opts_hover(css = "stroke:orange; stroke-width:1px;")  # <- slight darkening
  ))

## -----------------------------------------------------------------------------
combined_dmso <- combined %>%
  filter(Treatment_1 == "DMSO")

## ----multi_plates_plot_mds, fig.width = 8, fig.height = 4---------------------
plot_mds(combined, group_by = "orig.ident", label = "combined_id", n_labels = 30)
plot_mds(combined_dmso, group_by = "orig.ident", label = "combined_id", n_labels = 30)

## ----multi_plates_plot_rle, fig.width = 8, fig.height = 4---------------------
plot_rle(combined_dmso, label_column = "orig.ident", normalisation = "raw") + scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_rle(combined_dmso, label_column = "orig.ident", normalisation = "edgeR") + scale_x_discrete(drop = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----multi_plates_compute_single_de, fig.width = 8, fig.height = 6------------
treatment_samples <- "Staurosporine_0.1"
control_samples <- "DMSO_0"
subset <- combined[, grepl(paste0(treatment_samples, "|", control_samples), combined$combined_id)]
batch <- subset$orig.ident
combined_edgeR <- compute_single_de(combined, treatment_samples, control_samples, method = "edgeR", batch = batch)
plot_volcano(combined_edgeR)


## ----multiplate_plot_cpm, fig.width = 8, fig.height = 6-----------------------
genes <- combined_edgeR$gene[1:6]
group_by <- "combined_id"
plot_cpm(combined,genes, group_by, treatment_samples, control_samples)


## -----------------------------------------------------------------------------
summarise_de(combined_edgeR, lfc_threshold = 1, padj_threshold = 0.01, multi=FALSE)

## -----------------------------------------------------------------------------
treatments <- combined %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
combined <- compute_multi_de(combined, treatments, control_samples = "DMSO_0", method = "edgeR", num_cores = 2, batch = batch)

## -----------------------------------------------------------------------------
summarise_de(combined, lfc_threshold = 1, padj_threshold = 0.01, multi=TRUE)

## ----mutli_plates_plot_multi_de, fig.width=10, fig.height=6-------------------
plot_multi_de(combined, group_by = "combined_id", value = "log2FC", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")


