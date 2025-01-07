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

#join with metadata
mac <- mac %>%
  inner_join(metadata, by = c(".cell" = "Barcode"))

#calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA")

## ----violin_plot, fig.width = 8-----------------------------------------------
#example of Seurat function being used
VlnPlot(mac, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4)

## ----plate_layout, fig.width = 8----------------------------------------------
#QC plot plate layout (all metadata columns can be used):
plate_layout(mac, "nCount_RNA", "Sample_type")

## ----mds_plot, fig.width = 8--------------------------------------------------
#example of MDS function, using limma
plot_mds(mac, "Sample_type")

## ----rle_plot, fig.width = 8, fig.height=5------------------------------------
#RLE function
rle_plot(mac, label_column = "Row")
rle_plot(mac, label_column = "Row", normalisation = "limma_voom")

