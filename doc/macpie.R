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
#example of MDS function, using limma
plot_mds(mac, "Sample_type")

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
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

treatment_samples="Staurosporine_0.1"
control_samples<-"DMSO_0"

top_table_edgeR<-differential_expression(mac, treatment_samples, control_samples,method = "edgeR")
top_table_Seurat<-differential_expression(mac, treatment_samples, control_samples,method = "Seurat_wilcox")

plot_volcano(top_table_edgeR)
plot_volcano(top_table_Seurat)


