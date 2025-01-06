
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
docker run --rm -ti -e PASSWORD=yourpassword -p 8787:8787 lauratwomey/macpie
docker run -d --platform linux/amd64 -p 8787:8787 -e PASSWORD=YOURD -v /Users/mark/ML_Macpie_dev:/home/rstudio/ML_Macpie_dev lauratwomey/macpie
```

## Example

This is a basic example which shows you how to import data and perform basic QC.

``` r
library(macpie)

## basic example code
mac<-read_macseq("directory_of_filtered_gene_matrix","path_to_metadata_file")

## summary of metadata
summarize_metadata(mac)

## plate layout
plate_layout(mac)


```

This is a basic example which shows you how to import user defined metadata and basic QC of metadata.

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

#example of MDS function, using limma
plot_mds(mac, "Sample_type")

#RLE function
rle_plot(mac, label_column = "Row")

```
