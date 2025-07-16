
<!-- README.md is generated from README.Rmd. Please edit that file -->

# macpie

<img src="docs/logo.png" width="100" align="right" class="pkgdown-hide"/>

<!-- badges: start -->

[![R-CMD-check](https://github.com/PMCC-BioinformaticsCore/macpie/actions/workflows/r-cmd-check.yaml/badge.svg)](https://github.com/PMCC-BioinformaticsCore/macpie/actions/workflows/r-cmd-check.yaml)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15778812.svg)](https://doi.org/10.5281/zenodo.15778812)

<!-- badges: end -->

## Introduction

macpie is a R toolkit designed for researchers, originally with MAC-seq
data in mind, but validated for general High-Throughput Transcriptomics
(HTTr) data applications. Its primary aim is to deliver the latest tools
for quality control (QC), visualization, and analysis.

For processing raw sequencing data into count matrices, please refer to
our companion Nextflow workflow:  
[**dinoflow: Nextflow workflow for
MAC-seq**](https://github.com/PMCC-BioinformaticsCore/dinoflow)

## Documentation

Full documentation and step-by-step tutorials are available at:  
[**macpie documentation
site**](https://pmcc-bioinformaticscore.github.io/macpie/articles/macpie.html)

### Example Data

We provide both **full** and **subset** datasets for testing and
exploration:

- **Full example dataset** is hosted on Zenodo:  
  <https://doi.org/10.5281/zenodo.15778812>

- **Quick-start subset `mini_mac`** is bundled with the package and can
  be loaded directly in R:

``` r
data("mini_mac", package = "macpie")
```

## Requirements

- **R ≥ 4.3.3**

  The package is developed on R 4.3.3, and installation (including all
  dependencies) has been tested on the latest R 4.5.0.

### Installation and dependencies

All required R packages are automatically installed:

- via pak::pkg_install or devtools::install_github

- or use our pre-built Docker image for a ready-to-use environment.

#### Installing locally

To install the development version of `macpie`, we recommend using
`pak` - a fast package installer to install directly from GitHub. Make
sure you have the `pak` package installed first.

##### Configure your own GitHub PAT (optional)

`pak` uses the GitHub API to resolve versions and fetch metadata—and
unauthenticated requests are limited to 60 per hour. To avoid rate-limit
errors, you can add your own Personal Access Token (PAT):

``` r
# install gitcreds if you haven’t already
install.packages("gitcreds")

# this will prompt you to paste in your PAT
gitcreds::gitcreds_set()
```

``` r
# 1. Install pak (fast installer)
if (!requireNamespace("pak", quietly=TRUE)) {
  install.packages("pak")
}

# 2. Make sure Bioconductor repos are set
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
options(repos = BiocManager::repositories())

# 3. Re-install macpie + exactly its Depends & Imports via pak
pak::pkg_install(
  "PMCC-BioinformaticsCore/macpie",
  dependencies = c("Depends", "Imports")
)

# 4. Verify
# Should load without error:
library(macpie)

# Optional dependencies
# To unlock all macpie features, install the "Suggests" packages:
pak::pkg_install(c(
  "testthat", "devtools","leidenbase", "gridExtra", "variancePartition", 
  "rcdk", "SingleCellExperiment", "doParallel", "BiocParallel",
  "zinbwave", "EDASeq", "SummarizedExperiment", "mcprogress",
  "glmGamPoi", "fgsea", "webchem", "data.table", "PoiClaClu",
  "pheatmap", "purrr", "httr2", "readr", "patchwork",
  "scran", "Matrix", "umap", "ggrepel", "forcats",
  "gtools", "colorspace", "enrichR", "readxl", "MOFA2"
))

```

Optionally, you can also install the package by

``` r

# First, install devtools if not already installed
install.packages("devtools")

# Make sure BiocManager is installed, and point your repos at both CRAN + Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

options(repos = BiocManager::repositories())

# Install CRAN packages (Imports + Suggests) 
install.packages(c(
  # Imports
  "dplyr", "ggplot2", "rlang", "tidyseurat", "stringr", "tibble", "ggiraph",
  "tidyr", "scales", "drc", "tidyverse", "glue", "igraph", "unikn",
  # Suggests
  "leidenbase", "gridExtra", "variancePartition","rcdk", "doParallel", 
  "mcprogress", "fgsea", "webchem", "data.table","PoiClaClu", "pheatmap", 
  "purrr", "httr2", "readr", "patchwork", "Matrix", "umap", "ggrepel", 
  "forcats", "gtools", "colorspace", "enrichR", "readxl"
))

# Install Bioconductor packages (Imports + Suggests) 
BiocManager::install(c(
  # Imports
  "edgeR", "limma", "Biobase", "Seurat", "DESeq2", "RUVSeq",
  # Suggests
  "SingleCellExperiment", "BiocParallel", "zinbwave", "EDASeq",
  "SummarizedExperiment", "glmGamPoi", "scran"
))

# Install MOFA2 from its GitHub (it’s not on Bioconductor) - this is a Suggests package
devtools::install_github("bioFAM/MOFA2")

# Finally, install macpie itself
devtools::install_github("PMCC-BioinformaticsCore/macpie", dependencies = TRUE)
```

#### Using Docker image

Have your docker desktop running, open a terminal, paste the docker pull
command and install, depending on your platform.

1.  Pull the Docker image

``` bash
docker pull --platform linux/amd64 xliu81/macpie:v1.0.0
```

2.  Run the Docker container

``` bash
docker run --rm -ti \
  -e PASSWORD=password \
  -p 8787:8787 \
  --platform linux/amd64 \
  -v /path/to/your/macpie/:/home/rstudio/macpie:z\
  xliu81/macpie:v1.0.0
```

- Replace /path/to/your/macpie/ with the absolute path to your local
  repo.

- Copy and paste <http://localhost:8787> in your browser

  - Username: rstudio

  - Password: password (or the one you set in the docker run command)

After logging in, you’ll find your local directory mounted under:

``` bash
/home/rstudio/macpie/
```

## Quick start

In here we show using the `mini_mac` dataset for a quick start. From
each vignette on our website, we only include a couple of functions in
this quick start.

``` r
library(macpie)

# load mini_mac, 
# mini_mac is a tidySeurat object with matched metadata
data("mini_mac")

# Quality control
# Filter by counts per sample group
mini_mac <- filter_genes_by_expression(mini_mac,
                                  group_by = "combined_id",
                                  min_counts = 5,
                                  min_samples = 3)



# MDS plot
p <- plot_mds(mini_mac, group_by = "Sample_type", label = "combined_id", n_labels = 30)
girafe(ggobj = p, fonts = list(sans = "sans"))


# Correction of the batch effect
# First we will subset the data to look at control, DMSO samples only
mini_mac_dmso <- mini_mac %>%
  filter(Treatment_1 == "DMSO")

# Run the RLE function
plot_rle(mini_mac_dmso, label_column = "Row", normalisation = "limma_voom")


# Transcriptional analysis
# Single comparison
# First perform the differential expression analysis
treatment_samples <- "Staurosporine_10"
control_samples <- "DMSO_0"
top_table <- compute_single_de(mini_mac, treatment_samples, control_samples, method = "limma_voom")

top_genes <- top_table %>%
  filter(p_value_adj < 0.1) %>%
  select(gene) %>%
  pull()

# A volcano plot with very small number of genes, as it's a subset of the full dataset 
plot_volcano(top_table, max.overlaps = 16)

# Multiple comparisons
# Filter out lower concentrations of compounds and untreated samples
treatments <- mini_mac %>%
  filter(Concentration_1 == 10) %>%
  select(combined_id) %>%
  filter(!grepl("DMSO", combined_id)) %>%
  pull() %>%
  unique()
mini_mac <- compute_multi_de(mini_mac, treatments, control_samples = "DMSO_0", method = "limma_voom", num_cores = 1)

# plot shared differentially expressed genes
plot_multi_de(mini_mac, group_by = "combined_id", value = "log2FC", p_value_cutoff = 0.01, direction="up", n_genes = 5, control = "DMSO_0", by="fc")


```

## Citation

``` r
citation("macpie")
```

``` r
#>To cite the macpie package in publications, please use:

#>  Bartonicek N, Liu X, Twomey L, Meier M, Lupat R, Craig S, Yoannidis D, Li J, Semple T,
#>  Simpson K, Li M, Ramm S (2025). “macpie: a scalable workflow for high-throughput
#>  transcriptomic profiling.” Manuscript in preparation.

#>A BibTeX entry for LaTeX users is

#>  @Unpublished{,
#>    title = {macpie: a scalable workflow for high-throughput transcriptomic profiling},
#>    author = {Nenad Bartonicek and Xin Liu and Laura Twomey and Michelle Meier and Richard Lupat and Stuart Craig and David Yoannidis and #>Jason Li and Tim Semple and Kaylene J Simpson and Mark X Li and Susanne Ramm},
#>    organization = {Peter MacCallum Cancer Centre},
#>    year = {2025},
#>    note = {Manuscript in preparation},
#>  }
```
