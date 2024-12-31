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
project_metadata<-paste0(data_dir,project_name,"/",project_name,"_metadata.csv")
project_rawdata<-paste0(data_dir,project_name,"/raw_matrix")

#load metadata
######## 1. Mark's function
#metadata<-read_metadata()
#TO-DO: fix the function so that it doesn't open a window and can be used from the command line
#TO-DO: check for importing cell count & FACS data
#TO-DO: fix metadata from the other project

metadata<-read.csv(project_metadata,header=T) %>%
  arrange("Barcode") %>%
  mutate(Treatment_1=gsub("-","_",Treatment_1)) %>%
  mutate(Treatment_conc=paste0(Treatment_1,"_",Concentration_1))

#only create an id if there are multiple plates
#TO-DO: test a project with multiple plates
if(length(project_metadata)>1){
  metadata<-metadata %>%
    mutate(id=gsub("Plate","",Plate_ID)) %>%
    mutate(Barcode=paste0(id,"_",Barcode))
}

######## 2. Mark's function
#validate_metadata(metadata)

#import reads
# TO-DO: check for multiple folders of data that should be imported at the same time
raw_counts_total <- Read10X(data.dir = project_rawdata)
keep <- rowSums(cpm(raw_counts_total)>=10) >= 2
raw_counts <- raw_counts_total[keep,]

#create tidySeurat object
mac <- CreateSeuratObject(counts = raw_counts,
                          project = project_name,
                          min.cells = 1,
                          min.features = 1)

#calculate percent of mitochondrial and ribosomal genes
mac[["percent.mt"]] <- PercentageFeatureSet(mac, pattern = "^mt-|^MT-")
mac[["percent.ribo"]] <- PercentageFeatureSet(mac, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa|
                                                ^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA|^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#TO-DO: verify that the mouse and human regexp for ribo proteins is correct (grepl("^RP",row.names(mac@assays$RNA$counts)))

#add metadata to the real data
mac<- mac %>%
  inner_join(metadata,by=c(".cell"="Barcode"))

#QC plot plate layout (all metadata columns can be used):
plate_layout(mac,"nCount_RNA","Sample_type")

#example of Seurat function being used
VlnPlot(mac, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 4)

#example of MDS function, using limma
plot_mds(mac,"Treatment_1")

#RLE function
rle_plot(mac, label_column = "Row")

#TO-DO: verify that different plates would be plotted side-by-side


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
#exists("plot_plate_layout", where = globalenv(), inherits = FALSE)
#6. check and document
#check()
#7. click inside the function and
#Code > Insert roxygen skeleton.
#document()
#check()
#lint(filename="R/functionX.R",linters = custom_linters)
#8. make the test
#use_test("functionX")
#test()
#9. terminal: git
#git add .
#git commit -m "message"
#git push origin function_author
#git checkout main
#git branch -d function_author







