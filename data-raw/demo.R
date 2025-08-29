usethis::use_data(demo, overwrite = TRUE)
# data-raw/demo.R
# Create a tiny demo dataset for Quick start in README.
# Source of truth: derived from the packaged 'mini_mac' to keep schema consistent.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(usethis)
})

# Load the existing packaged object
# (works when run from the package root after installing/loading deps)
load("data/mini_mac.rda")  # provides object 'mini_mac'

# Choose a small subset (edit sizes to taste)
set.seed(42)
keep_genes <- sample(rownames(mini_mac), size = min(200, nrow(mini_mac)))    # smaller gene set

demo <- subset(mini_mac, features = keep_genes)

# Optional: slim metadata to the columns you use in the README
cols_keep <- c(".cell", "Sample_type", "Treatment_1", "Concentration_1", "Project")
demo@meta.data <- demo@meta.data[, intersect(cols_keep, colnames(demo@meta.data)), drop = FALSE]

# Optional: store a short description to help users
demo@misc$description <- "Tiny demo Seurat object derived from mini_mac for README quick-start examples."

# Save to data/ as a compressed .rda
usethis::use_data(demo, overwrite = TRUE, compress = "xz")  # creates data/new_demo.rda
