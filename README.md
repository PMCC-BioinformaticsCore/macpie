
<!-- README.md is generated from README.Rmd. Please edit that file -->

# macpie

<!-- badges: start -->
<!-- badges: end -->

The goal of macpie package is to provide the users of Mac-seq data with
the most recent set of tools for QC, visualisation and analysis for
this high-throughput transcriptomic platform.

## Installation guide

macpie can be installed from github using the install_github() function from the package remotes:

``` r
remotes::install_github("https://github.com/PMCC-BioinformaticsCore/macpie")
```

## Dependencies

The simplest way is to use our docker container with all the R packages installed.
To install the docker container:
``` r
docker pull xliu81/macpie::latest
```
To use the docker container:
Inside your docker desktop, open a terminal, and run:
``` r
docker run --rm -ti -d -e PASSWORD=yourpassword -p 8787:8787 -v /Users/username/macpie:/home/rstudio/macpie xliu81/macpie
```
You may need to specify the platform:
```
docker run --rm -ti -d --platform linux/amd64 -p 8787:8787 -e PASSWORD=yourpassword -v /Users/username/macpie:/home/rstudio/macpie xliu81/macpie

```
Note, the :z in the -v argument of the docker command is essential.

Then, to use macpie interactively, navigate to http://localhost:8787/
Username: rstudio
Password: yourpassword

You are now ready to use macpie!

## Example use

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

#example of MDS function, using limma
plot_mds(mac, "Sample_type")

#RLE function
rle_plot(mac, label_column = "Row")

```
