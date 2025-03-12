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
library(SingleCellExperiment)

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
project_metadata <- file.path(data_dir, project_name, paste0(project_name, "_metadata.csv"))
metadata <- read_metadata(project_metadata)

################## metadata ##################
# Mark's load metadata
project_metadata <- paste0(data_dir, project_name,"/",project_name,"_metadata.csv")
metadata<-read_metadata(project_metadata)

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
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_"))


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
girafe(ggobj = p)

# Compare normalisation methods using the RLE function
#To-do: verify that different plates would be plotted side-by-side
#To-do: add all Xin's plots
#To-do: QC summary

# RLE plot
mac_dmso <- mac %>%
  filter(Treatment_1 == "DMSO")
plot_rle(mac_dmso, label_column = "Row", normalisation = "edgeR")

# Differential expression ===================================================

## DE Single -----------------

treatment_samples <- "Staurosporine_0.1"
control_samples <- "DMSO_0"

# Perform differential expression
top_table <- run_differential_expression(mac, treatment_samples, control_samples, method = "edgeR")
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
enriched <- pathway_enrichment(top_genes, "MSigDB_Hallmark_2020", species = "human")
plotEnrich(enriched)

## Screen-level analysis --------------------------

### DGE analysis ------------------
#to-do: shorten
treatments <- mac %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()

##slow version
#de_list<-sapply(treatments,function(x){
#  differential_expression(mac, x, control_samples,method = "edgeR");
#  cat(".")
#})

# Load genesets from enrichr for a specific species or define your own
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")

# Update the mac object with differential expression
mac <- compute_multi_DE(mac, treatments, control_samples = "DMSO_0", method = "edgeR")
mac <- compute_enrichr_pathways(mac, genesets = enrichr_genesets)
mac <- compute_multi_screen_profile(mac, target = "Staurosporine_10")
mac <- prepare_de_umap(mac)
mac <- find_clusters_de_umap(mac, k = 3)
p <- multi_plot_umap(mac, group_by = "cluster", max_overlaps = 5)
girafe(ggobj = p)

# Get all the differential expression information in a tabular format
de_genes_per_comparison <- bind_rows(mac@tools$diff_exprs)
enriched_pathways_per_comparison <- mac@tools$pathway_enrichment
screen_preeofile_per_comparison <- mac@tools$screen_profile

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
#exists("plote_plate_layout", where = globalenv(), inherits = FALSE)
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


