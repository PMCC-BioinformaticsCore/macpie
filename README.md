
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
remotes::install_github("https://github.com/qoiopipq/macpie")
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
2024 Nov findmetadata() function remains user interactive

``` r
library(macpie)
## basic example of using find_metadata()
## user should provide a dir path
metadata <- find_metadata("path_to_metadata_file")

##if no dir path is provided, it will trigger a file.choose() for user to find file
metadata <- find_metadata()

## summary of metadata QC
validate_metadata(metadata)

```
