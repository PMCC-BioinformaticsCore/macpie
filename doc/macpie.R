## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  message = FALSE
)

## ----load_metadata------------------------------------------------------------
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

#load all functions
devtools::load_all()

# Define project variables
project_name <- "PMMSq033"
project_metadata <- system.file("extdata/PMMSq033/PMMSq033_metadata.csv", package = "macpie")

# Load metadata
metadata <- read_metadata(project_metadata)
colnames(metadata)

# Validate metadata
validate_metadata(metadata)

## ----metadata_plot, fig.width = 8, fig.height = 6-----------------------------
metadata_heatmap(metadata)


## ----load_data----------------------------------------------------------------
project_rawdata <- system.file("extdata/PMMSq033/raw_matrix", package = "macpie")
raw_counts_total <- Read10X(data.dir = project_rawdata)
keep <- rowSums(cpm(raw_counts_total) >= 10) >= 2
raw_counts <- raw_counts_total[keep, ]

#create tidySeurat object
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 1,
                          min.features = 1)


## ----violin_plot, fig.width = 8-----------------------------------------------
#calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA")

#join with metadata
mac <- mac %>%
  inner_join(metadata, by = c(".cell" = "Barcode"))

#add unique identifier
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

#example of a function from Seurat QC 
VlnPlot(mac, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), 
        ncol = 4, group.by="Sample_type")

## ----subset_seurat, fig.width = 8---------------------------------------------
unique(mac$Project)
mac <- mac %>%
  filter(Project == "Current")

#QC plot plate layout (all metadata columns can be used):
plate_layout(mac, "nCount_RNA", "Sample_type")

## ----mds_plot, fig.width = 8--------------------------------------------------
#example of MDS function 
p<-plot_mds(mac)
girafe(ggobj = p)

## ----rle_plot, fig.width = 8, fig.height=5------------------------------------

# First we will subset the data to look at control, DMSO samples only
mac_dmso <- mac %>%
  filter(Treatment_1 == "DMSO")

#RLE function
rle_plot(mac_dmso, label_column = "Row")
rle_plot(mac_dmso, label_column = "Row", normalisation = "SCT")
rle_plot(mac_dmso, label_column = "Row", normalisation = "edgeR")

## ----de_analysis, fig.width = 8, fig.height=5---------------------------------

# First perform the differential expression analysis
treatment_samples="Staurosporine_0.1"
control_samples<-"DMSO_0"

top_table_edgeR<-differential_expression(mac, treatment_samples, control_samples,method = "edgeR")

plot_volcano(top_table_edgeR)


## ----pathway_analysis_single, fig.width = 8, fig.height=15--------------------

top_genes<-top_table_edgeR %>%
  filter(p_value_adj<0.05) %>%
  select(gene) %>%
  pull()

enriched <- enrichR::enrichr(top_genes, c("MSigDB_Hallmark_2020","DisGeNET",
                                 "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO"))
p1<-enrichR::plotEnrich(enriched[[1]])
p2<-enrichR::plotEnrich(enriched[[2]])
p3<-enrichR::plotEnrich(enriched[[3]])

gridExtra::grid.arrange(p1, p2, p3, ncol = 1)

## ----de_multi, fig.width = 8, fig.height=5------------------------------------
de_list <- multi_DE(data = mac[50:150], 
                    treatment_samples = NULL, 
                    control_samples = "DMSO_0", 
                    method = "edgeR")

## ----enriched_pathways, fig.width = 8, fig.height=5---------------------------
de_df <- bind_rows(de_list)

#load genesets for human MSigDB_Hallmark_2020
file_path <- system.file("extdata", "PMMSq033/pathways.Rds", package = "macpie")
genesets <- readRDS(file_path)

enriched_pathways <- de_df %>%
  filter(p_value_adj<0.01) %>% #select DE genes based on FDR < 0.01
  group_by(combined_id) %>%
  filter(n_distinct(gene)>5) %>% #filter out samples with less than 5 DE genes
  reframe(enrichment=hyper_enrich_bg(gene, genesets=.env$genesets,background = "human")) %>%
  unnest(enrichment)

enriched_pathways_mat <- enriched_pathways %>%
  select(combined_id, Term, Combined.Score) %>%
  pivot_wider(
    names_from = combined_id,
    values_from = Combined.Score
  ) %>%
  column_to_rownames(var = "Term") %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, log1p(.)))) %>%  # Replace NA with 0 across all columns
  as.matrix()

pheatmap(enriched_pathways_mat)

