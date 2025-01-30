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

#define longer length for description files
custom_linters <- lintr::linters_with_defaults(
  line_length_linter = lintr::line_length_linter(120) # Set max line length to 120
)

#load all previous libraries
load_all()

#select which sample dataset to use
project_name<-"PMMSq033"

#directory with data
data_dir<-"inst/extdata/"

################## metadata ##################
# Mark's load metadata
project_metadata<-paste0(data_dir,project_name,"/",project_name,"_metadata.csv")
metadata<-read_metadata(project_metadata)

# Mark's validate metadata
validate_metadata(metadata)

# Susi's function: heatmap visualisation of metadata
metadata_heatmap(metadata)

#only create an id if there are multiple plates
#to-do: test multiple plates
if(length(project_metadata)>1){
  metadata<-metadata %>%
    mutate(id=gsub("Plate","",Plate_ID)) %>%
    mutate(Barcode=paste0(id,"_",Barcode))
}

################## reads to Seurat object ##################
# TO-DO: check for multiple folders of data that should be imported at the same time
project_rawdata<-paste0(data_dir,project_name,"/raw_matrix")
raw_counts_total <- Read10X(data.dir = project_rawdata)
keep <- rowSums(cpm(raw_counts_total)>=10) >= 2
raw_counts <- raw_counts_total[keep,]

#create tidySeurat object
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 5,
                          min.features = 5)

#calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|
                                                ^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#TO-DO: verify that the mouse and human regexp for ribo proteins is correct (grepl("^RP",row.names(mac@assays$RNA$counts)))

#add metadata to the real data
mac<- mac %>%
  inner_join(metadata,by=c(".cell"="Barcode")) %>%
  filter(Project == "Current")

#create an ID that uniquely identifies samples based on the
#combination of treatment and treatment concentration
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_"))


################## QC ##################

#QC plot plate layout (all metadata columns can be used):
plate_layout(mac,"nCount_RNA","Sample_type")

#example of Seurat function being used
VlnPlot(mac,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"),
        group.by = "Sample_type",
        ncol = 4)

#example of MDS function, using limma
p<-plot_mds(mac)
girafe(ggobj = p)

#Compare normalisation methods using the RLE function
#To-do: verify that different plates would be plotted side-by-side
#To-do: add all Xin's plots
#To-do: QC summary

# Xin's rle plot
mac_dmso<- mac %>%
  filter(Treatment_1=="DMSO")
rle_plot(mac_dmso, label_column = "Row",normalisation="DESeq2")

################## Differential expression ##################

################## DE Single ##################

treatment_samples="Staurosporine_0.1"
control_samples<-"DMSO_0"

#perform differential expression
top_table<-differential_expression(mac, treatment_samples, control_samples,method = "limma_voom")
plot_volcano(top_table)

#perform pathway enrichment
top_genes<-top_table %>%
  filter(p_value_adj<0.01) %>%
  select(gene) %>%
  pull()

#basic pathway enrichment
enrichr_results <- enrichr(top_genes, "MSigDB_Hallmark_2020")
plotEnrich(enrichr_results[[1]])

enriched <- pathway_enrichment(top_genes, "MSigDB_Hallmark_2020", species = "human")
plotEnrich(enriched)

################## DE Multiple ##################

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

num_cores <- detectCores() - 2
de_list <- multi_DE(mac, treatments, control_samples, num_cores=num_cores, method = "edgeR")
de_df <- do.call("rbind", de_list)

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

########## Screen for similarity of profiles

fgsea_results <- screen_profile(de_list, target = "Staurosporine_10", n_genes_profile = 500)
fgsea_results %>%
  mutate(logPadj=c(-log10(padj))) %>%
  arrange(desc(NES)) %>%
  mutate(target = factor(target, levels = unique(target))) %>%
  ggplot(.,aes(target,NES))+
  #geom_point(aes(size = logPadj)) +
  geom_point() +
  facet_wrap(~pathway,scales = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




df<-do.call("rbind",de_list) %>%
  mutate(comparison=combined_id)
df_wide <- df %>%
  select(gene, comparison, t) %>%
  pivot_wider(names_from = comparison, values_from = t)

set.seed(10)
df_umap<-umap(t(df_wide[,-1]))


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
# lint(filename="R/functionX.R",linters = custom_linters)
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



