# to create small example data from a plate
# the small example data has 384 wells and 500 genes
library(devtools)
library(Matrix)

in_dir <- "inst/extdata/PMMSq033/raw_matrix/"
mat <- readMM(file.path(in_dir, "matrix.mtx.gz"))
bar <- readLines(gzfile(file.path(in_dir, "barcodes.tsv.gz")))
feat <- read.delim(gzfile(file.path(in_dir, "features.tsv.gz")),header = FALSE, stringsAsFactors = FALSE)

#barcodes
colnames(mat) <- bar
rownames(mat) <- feat$V2

#pick 500 genes
set.seed(1)
keep_genes<- sample(nrow(mat), 500)

#drop genes to 300 
mat_small <- mat[keep_genes, ,drop = FALSE]
bar_small <- bar
feat_small <- feat[keep_genes, ]

#create the dense small matrix 
dense_small <- as.matrix(mat_small)
rownames(dense_small) <- make.unique(rownames(dense_small))
rownames(dense_small) <- gsub("_", "-", rownames(dense_small), fixed = TRUE)
write.table(dense_small, 
            file = "inst/extdata/PMMSq033/PMMSq033_small.txt", 
            sep = "\t", quote = FALSE, col.names = NA)



# write three files for small mac
out_dir <- file.path("inst/extdata/PMMSq033/mini_mac/raw" )        # …/extdata/PMMSq033/mini_mac
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

writeMM(mat_small,        file.path(out_dir, "matrix.mtx.gz"))          # sparse counts
writeLines(bar_small,     file.path(out_dir, "barcodes.tsv.gz"))        # cells
write.table(feat_small[ , 1:2],  file = file.path(out_dir, "features.tsv.gz"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# read in the three input files for the small example data to save it as .rda 

# generate three slots for storing results from DE, pathway and UMAP from DE
# for the small example data
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
library(data.table)
devtools::load_all()


# Load metadata
project_name <- "PMMSq033_mini"
metadata <- read_metadata("inst/extdata/PMMSq033/PMMSq033_metadata_drugnames.csv")
metadata$Time <- as.factor(metadata$Time)
metadata$Concentration_1 <- as.factor(metadata$Concentration_1)
colnames(metadata)

# Validate metadata
validate_metadata(metadata)
mini_mac_in <-"inst/extdata/PMMSq033/mini_mac/raw/"
mini_mac_counts <- Read10X(data.dir = mini_mac_in)
mini_mac <- CreateSeuratObject(counts = mini_mac_counts,
                               project = project_name)
mini_mac <- NormalizeData(mini_mac)

# join with metadata
mini_mac <- mini_mac %>%
  inner_join(metadata, by = c(".cell" = "Barcode"))

# rownames(metadata) <- metadata$Barcode
# mini_mac <- Seurat::AddMetaData(mini_mac, metadata = metadata)


mini_mac <- mini_mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))

mini_mac@meta.data[["percent.mt"]] <- PercentageFeatureSet(mini_mac, pattern = "^mt-|^MT-")
mini_mac@meta.data[["percent.ribo"]] <- PercentageFeatureSet(mini_mac, pattern = "^Rp[slp][[:digit:]]|^Rpsa|^RP[SLP][[:digit:]]|^RPSA")


# Example of a function from Seurat QC 
VlnPlot(mini_mac, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 4, group.by = "Sample_type") & 
  scale_fill_manual(values = macpie_colours$discrete) 

mini_mac <- mini_mac %>%
  filter(Project == "Current")



#multiple DE for 10 treatments with highest concentration 
mini_mac$combined_id <- make.names(mini_mac$combined_id)

treatments <- mini_mac@meta.data %>%
  filter(Concentration_1 == 10) %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
# 10 treatments are 
#  [1] "Staurosporine_10"       "Paclitaxel_10"          "Chlorambucil_10"       
# [4] "Vinblastine_sulfate_10" "Etoposide_10"           "Cytarabine_10"         
# [7] "Camptothecin_10"        "Anastrozole_10"         "Sb590885_10"           
# [10] "Fluvastatin_sodium_10"  "Capivasertib_10"            "Nutlin.3a_10"              
# [13] "Ceralasertib_10"            "Erlotinib_hydrochloride_10" "Mk.2206_dihydrochloride_10"
mini_mac <- compute_multi_de(mini_mac, treatments[1:15], control_samples = "DMSO_0", method = "limma_voom", num_cores = 1)


#plot_volcano(mini_mac@tools[["diff_exprs"]][["Staurosporine_10"]])

# multi pathway 
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
mini_mac <- compute_multi_enrichr(mini_mac, genesets = enrichr_genesets)

usethis::use_data(mini_mac, overwrite = TRUE)


# generate pathway dataset for hypergeometric enrichment test
# convert the exiting pathways.Rds into data/
file_path <- system.file("extdata", "PMMSq033/pathways.Rds", package = "macpie")
genesets <- readRDS(file_path)
usethis::use_data(genesets)



