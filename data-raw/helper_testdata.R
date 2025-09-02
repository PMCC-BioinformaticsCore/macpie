# to create small test example data from a plate
# the small example data has 384 wells and 200 genes from the mini_mac dataset
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(usethis)
  library(macpie)
})

# Load the existing packaged object
# (works when run from the package root after installing/loading deps)
load("data/mini_mac.rda")  # provides object 'mini_mac'

# Choose a small subset (edit sizes to taste)
set.seed(42)
keep_genes <- sample(rownames(mini_mac), size = min(200, nrow(mini_mac)))    # smaller gene set

testdata <- subset(mini_mac, features = keep_genes)

# Store a short description to help users
testdata@misc$description <- "Tiny testdata Seurat object derived from mini_mac for testing."

# multiple DE for 10 treatments with highest concentration 
testdata@meta.data$combined_id <- make.names(paste0(testdata$Treatment_1, "_",testdata$Concentration_1))

# change orig.ident
testdata$orig.ident <- "testdata"

treatments <- testdata@meta.data %>%
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
testdata <- compute_multi_de(testdata, treatments[1:15], control_samples = "DMSO_0", method = "limma_voom", num_cores = 1)


#plot_volcano(testdata@tools[["diff_exprs"]][["Staurosporine_10"]])

# multi pathway 
enrichr_genesets <- download_geneset("human", "MSigDB_Hallmark_2020")
testdata <- compute_multi_enrichr(testdata, genesets = enrichr_genesets)

# save it as rds for testdata
saveRDS(testdata, file = "tests/testthat/testdata/testdata.rds")


