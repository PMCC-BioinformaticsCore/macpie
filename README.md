
<!-- README.md is generated from README.Rmd. Please edit that file -->

# macpie

<!-- badges: start -->
<!-- badges: end -->

The goal of macpie package is to provide the users of Mac-seq data with
the most recent set of tools for QC, visualisation and analysis for
this high-throughput transcriptomic platform.

## Installation

You can install the development version of macpie like so:

``` r
remotes::install_github("https://github.com/PMCC-BioinformaticsCore/macpie")
```

## Dependencies

The simplest way is to use our docker container with all the R packages installed.
``` r
docker pull lauratwomey/macpie
```
Inside your docker desktop, open a terminal, paste the docker pull command and install, depending on your platform.
``` r
docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 -v /Users/username/macpie:/home/rstudio/macpie lauratwomey/macpie
docker run -d --platform linux/amd64 -p 8787:8787 -e PASSWORD=yourpassword -v /Users/username/macpie:/home/rstudio/macpie lauratwomey/macpie
```

## Example

This is a basic example which shows you how to import user-defined data and metadata with some basic QC.

``` r
library(macpie)

#load metadata
project_metadata <- system.file("extdata/PMMSq033/PMMSq033_metadata.csv", package = "macpie")

# Load metadata
metadata <- read_metadata(project_metadata)

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
  
#add unique identifier
mac <- mac %>%
  mutate(combined_id = str_c(Treatment_1, Concentration_1, sep = "_")) %>%
  mutate(combined_id = gsub(" ", "", .data$combined_id))  

#example of MDS function, using limma
plot_mds(mac, "Sample_type")

#RLE function
rle_plot(mac, label_column = "Row")

#Calculate qc metrics for each condition
qc_stats <- QC_metrics(mac, "combined_id")
qc_stats$stats_summary

#Plot qc metrics 
plot_QC_metrics(qc_stats, "combined_id", "z_score")

#Plot distance matrix in a heatmap for condition and control with high varibility from previous qc metrics step 
plot_distance(mac, "combined_id", "Staurosporine_10")
plot_distance(mac, "combined_id", "DMSO_0")





```
